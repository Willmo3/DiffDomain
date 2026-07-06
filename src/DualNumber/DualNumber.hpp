//
// Created by will on 3/12/26.
//

#ifndef DIFFDOMAIN_DUALNUMBER_H
#define DIFFDOMAIN_DUALNUMBER_H

// NOTE: synthesized zonotopes are an early feature for autodiff!
#define USE_SYNTHESIZED_ZONOS

#include <cstdint>
#include <ostream>
#include <ranges>
#include <vector>
#include <complex>
#include <tuple>

#include "cereal/cereal.hpp"
#include "Numeric.hpp"
#include "polylib/polylib.hpp"

/**
 * Forward-mode automatic differentiation via dual numbers.
 *
 * @tparam T A type satisfying the Numeric concept.
 * @author Will Morris
 */
template<typename T>
requires Numeric<T>
class DualNumber {
public:
    /*
     * Constructors
     */
    DualNumber() = default;
    /**
     * @param primal Primal (function) value. Derivative value is initialized to one, the derivative of a single variable.
     */
    explicit DualNumber(double primal) : _primal_value(T(primal)), _deriv_value(T(1)) {}
    /**
     * @param primal Primal (function) value.
     * @param deriv  Derivative value.
     */
    DualNumber(T primal, T deriv): _primal_value(primal), _deriv_value(deriv) {}

    /*
     * Accessors
     */

    /**
     * @return Radius of the primal value.
     */
    double radius() const {
        return _primal_value.radius();
    }
    /**
     * @return Minimum of the primal value.
     */
    double min() const {
        return _primal_value.min();
    }
    /**
     * @return Maximum of the primal value.
     */
    double max() const {
        return _primal_value.max();
    }

    /**
     * @return A reference to the primal value. Less expensive, but limited lifetime.
     */
    const T &primal_ref() const {
        return _primal_value;
    }
    /**
     * @return A reference to the derivative value. Less expensive, but limited lifetime.
     */
    const T &deriv_ref() const {
        return _deriv_value;
    }

    /**
     * @return A copy of the primal value.
     * Note that this is shallow, and heap resources -- i.e. vector contents -- will be shared.
     */
    [[nodiscard]] T primal_value() {
        return _primal_value;
    }
    /**
     * @return a deep copy of the derivative value.
     * Note that this is shallow, and heap resources -- i.e. vector contents -- will be shared.
     */
    [[nodiscard]] T deriv_value() {
        return _deriv_value;
    }

    /*
     * DualNumber–DualNumber arithmetic.
     * Derivative rules applied:
     *   Addition/subtraction : sum rule
     *   Multiplication       : product rule
     *   Division             : quotient rule
     */
    [[nodiscard]] DualNumber operator+(const DualNumber& rhs) const {
        return { _primal_value + rhs._primal_value, _deriv_value + rhs._deriv_value };
    }
    [[nodiscard]] DualNumber operator-(const DualNumber& rhs) const {
        return { _primal_value - rhs._primal_value, _deriv_value - rhs._deriv_value };
    }
    [[nodiscard]] DualNumber operator*(const DualNumber& rhs) const {
        return { _primal_value * rhs._primal_value,
                 _deriv_value * rhs._primal_value + _primal_value * rhs._deriv_value };
    }
    [[nodiscard]] DualNumber operator/(const DualNumber& rhs) const {
        return {
            _primal_value / rhs._primal_value,
            (_deriv_value * rhs._primal_value - _primal_value * rhs._deriv_value) / rhs._primal_value.pow(2u)
        };
    }

    /*
     * DualNumber–scalar arithmetic
     */
    [[nodiscard]] DualNumber operator+(double scalar) const {
        return { _primal_value + scalar, _deriv_value };
    }
    [[nodiscard]] DualNumber operator-(double scalar) const {
        return { _primal_value - scalar, _deriv_value };
    }
    [[nodiscard]] DualNumber operator*(double scalar) const {
        return { _primal_value * scalar, _deriv_value * scalar };
    }
    [[nodiscard]] DualNumber operator/(double scalar) const {
        return { _primal_value / scalar, _deriv_value / scalar };
    }

    /*
     * Unary operations
     */

    /**
     * Applies the power rule: d/dx(f^n) = n * f^(n-1) * f'.
     * @param n Non-negative integer exponent.
     * @return A new dual number representing the result of power rule application
     */
    [[nodiscard]] DualNumber pow(uint32_t n) const {
        if (n == 0) {
            // f^0 = 1 (constant), derivative = 0
            return DualNumber(_primal_value.pow(0), _deriv_value * 0.0);
        }

        return { _primal_value.pow(n), _primal_value.pow(n - 1) * static_cast<double>(n) * _deriv_value };
    }
    /**
     * Applies power rule with exponent of 1/2. d/dx (sqrt(f)) = 1/2 f^(-1/2) * f'
     * @return A new dual number representing the result of sqrt with chain rule applied.
     */
    [[nodiscard]] DualNumber sqrt() const {
        return { _primal_value.sqrt(), _deriv_value / (_primal_value.sqrt() * 2.0) };
    }
    /**
     * @return A new dual number representing the absolute value, d/dx(|f|) = sign(f) * f'.
     */
    [[nodiscard]] DualNumber abs() const {
        if (_primal_value >= 0.0) {
            return { _primal_value.abs(), _deriv_value };
        } else if (_primal_value < 0.0) {
            return { _primal_value.abs(), _deriv_value * -1.0 };
        }
        // Mixed sign (e.g. an interval spanning zero): conservative bound.
        return { _primal_value.abs(), _deriv_value.abs() };
    }
    /**
     * @return A new dual number representing the result of exp with chain rule appconstlied.
     */
    [[nodiscard]] DualNumber exp() const {
        auto exp_primal = _primal_value.exp();
        return { exp_primal, exp_primal * _deriv_value };
    }
    /**
     * @return A new dual number representing the result of tanh with chain rule applied.
     * Derivative: 1 - tanh^2(f) * f'
     */
    [[nodiscard]] DualNumber tanh() const {
        auto tanh_primal = _primal_value.tanh();
        auto tanh_derivative = tanh_primal.pow(2u) * -1 + 1;
        return { tanh_primal, tanh_derivative * _deriv_value };
    }
    /**
     * @return A new dual number representing the result of sigmoid with chain rule applied.
     * Derivative: sigmoid(f) * (1 - sigmoid(f)) * f'
     */
    [[nodiscard]] DualNumber sigmoid() const {
        auto sigmoid_primal = _primal_value.sigmoid();
        auto sigmoid_derivative = sigmoid_primal * (sigmoid_primal * -1 + 1);
        return { sigmoid_primal, sigmoid_derivative * _deriv_value };
    }
    /**
     * @return A new dual number representing the result of ReLU with chain rule applied.
     * Derivative: 0 if f < 0, 1 if f >= 0, or conservative bound for mixed sign intervals.
     */
    [[nodiscard]] DualNumber relu() const {
        if (_primal_value >= 0.0) {
            return { _primal_value.relu(), _deriv_value };
        } else if (_primal_value < 0.0) {
            return { _primal_value.relu(), _deriv_value * 0.0 };
        }
        // Mixed sign (e.g. an interval spanning zero): conservative bound.
        return { _primal_value.relu(), _deriv_value.abs() };
    }

    /*
     * Compositional operations
     */

    /**
     * @param rhs Other DualNumber to union with.
     * @return A new DualNumber whose components are the respective unions.
     */
    [[nodiscard]] DualNumber union_with(const DualNumber &rhs) const {
        return { _primal_value.union_with(rhs._primal_value),
                 _deriv_value.union_with(rhs._deriv_value) };
    }
    /**
     * @param n_splits Number of splits to produce.
     * @return Vector of n_splits DualNumbers.
     */
    [[nodiscard]] std::vector<DualNumber> split(uint32_t n_splits) const {
        auto primal_splits = _primal_value.split(n_splits);
        auto deriv_splits  = _deriv_value.split(n_splits);

        std::vector<DualNumber> result;
        result.reserve(n_splits);

        for (uint32_t i = 0; i < n_splits; ++i) {
            result.emplace_back(primal_splits[i], deriv_splits[i]);
        }
        return result;
    }

    /*
     * DualNumber–DualNumber comparison operators
     *
     * Equality (==, !=) compares both primal and derivative.
     * Ordering (<, <=, >, >=) is determined by the primal value alone.
     */
    bool operator==(const DualNumber &rhs) const {
        return _primal_value == rhs._primal_value && _deriv_value == rhs._deriv_value;
    }
    bool operator!=(const DualNumber &rhs) const {
        return _primal_value != rhs._primal_value || _deriv_value != rhs._deriv_value;
    }
    bool operator<(const DualNumber &rhs) const {
        return _primal_value < rhs._primal_value;
    }
    bool operator<=(const DualNumber &rhs) const {
        return _primal_value <= rhs._primal_value;
    }
    bool operator>(const DualNumber &rhs) const {
        return _primal_value > rhs._primal_value;
    }
    bool operator>=(const DualNumber &rhs) const {
        return _primal_value >= rhs._primal_value;
    }

    /*
     * DualNumber–scalar comparison operators
     *
     * Comparison is against the primal value only
     */
    bool operator<(double scalar) const {
        return _primal_value < scalar;
    }
    bool operator<=(double scalar) const {
        return _primal_value <= scalar;
    }
    bool operator>(double scalar) const {
        return _primal_value > scalar;
    }
    bool operator>=(double scalar) const {
        return _primal_value >= scalar;
    }

    /*
     * Utility functions
     */
    template<class Archive>
    void serialize(Archive &archive) {
        archive(cereal::make_nvp("primal", _primal_value),
                cereal::make_nvp("deriv",  _deriv_value));
    }

    friend std::ostream &operator<<(std::ostream &os, const DualNumber &dn) {
        os << "DualNumber(primal=" << dn._primal_value << ", deriv=" << dn._deriv_value << ")";
        return os;
    }

private:
    T _primal_value;
    T _deriv_value;
};

#ifdef USE_SYNTHESIZED_ZONOS
#include "MixedForm/MixedForm.hpp"
// Function-specific synthesis helpers
#include "synthesis/tanh.hpp"
#include "synthesis/quotient.hpp"
#include "synthesis/product.hpp"

inline Eigen::MatrixXd get_corners(const Winterval &bounds1, const Winterval &bounds2, const Winterval &bounds3, const Winterval &bounds4) {
    auto x1 = Eigen::VectorXd(2);
    auto x2 = Eigen::VectorXd(2);
    auto x3 = Eigen::VectorXd(2);
    auto x4 = Eigen::VectorXd(2);
    x1 << bounds1.min(), bounds1.max();
    x2 << bounds2.min(), bounds2.max();
    x3 << bounds3.min(), bounds3.max();
    x4 << bounds4.min(), bounds4.max();
    return cartesian_product(x1, x2, x3, x4);
}

// synthesized abstract transformers for autodiff derived from Pasado
// https://dl.acm.org/doi/pdf/10.1145/3622867


/**
 * Synthesized abstract transformer for autodiff of product rule.
 * We will use standard operations for the primal, but abstract the derivative.
 */
template<>
[[nodiscard]] inline DualNumber<MixedForm> DualNumber<MixedForm>::operator*(const DualNumber&rhs) const {
    auto corners = get_corners(_primal_value.interval_bounds(), rhs._primal_value.interval_bounds(), _deriv_value.interval_bounds(), rhs._deriv_value.interval_bounds());

    auto affine_rep = compute_opt_affine_product(*this, rhs, corners);
    auto interval_rep = compute_opt_interval_product(corners);
    auto derivative = MixedForm(affine_rep, interval_rep);

    return {
        _primal_value * rhs._primal_value,
        derivative
    };
}

/**
 * Synthesized abstract transformer for autodiff of quotient rule.
 */
template<>
[[nodiscard]] inline DualNumber<MixedForm> DualNumber<MixedForm>::operator/(const DualNumber &rhs) const {
    // Default to regular division if we have a range w/ negative values.
    if (rhs.min() <= 0 && rhs.max() >= 0) {
        return {
            _primal_value / rhs._primal_value,
            (_deriv_value * rhs._primal_value - _primal_value * rhs._deriv_value) / rhs._primal_value.pow(2u)
        };
    }

    auto corners = get_corners(_primal_value.interval_bounds(), rhs._primal_value.interval_bounds(), _deriv_value.interval_bounds(), rhs._deriv_value.interval_bounds());

    auto affine_result = compute_opt_affine_quotient(corners, *this, rhs);
    auto refined_interval = compute_opt_interval_quotient(corners, rhs._primal_value.min(), rhs._primal_value.max());
    auto derivative = MixedForm(affine_result, refined_interval);

    return {
        _primal_value / rhs._primal_value,
        derivative
    };
}

/**
 * Synthesized chain rule abstract transformer for the tanh function, from Pasado.
 * @return A new dual number representing the result of tanh with chain rule applied.
 *
 *
 */
template<>
inline DualNumber<MixedForm> DualNumber<MixedForm>::tanh() const {
    if (_deriv_value.min() == 0 && _deriv_value.max() == 0) {
        return { _primal_value.tanh(), MixedForm(0) };
    }

    auto affine_rep = compute_opt_affine_tanh(*this);
    auto interval_rep = compute_opt_interval_tanh(*this);

    return { _primal_value.tanh(), MixedForm(affine_rep, interval_rep) };
}

#endif

#endif //DIFFDOMAIN_DUALNUMBER_H
