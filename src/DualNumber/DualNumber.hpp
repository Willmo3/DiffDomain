//
// Created by will on 3/12/26.
//

#ifndef DIFFDOMAIN_DUALNUMBER_H
#define DIFFDOMAIN_DUALNUMBER_H

#include <cstdint>
#include <ostream>
#include <vector>

#include "cereal/cereal.hpp"
#include "Numeric.hpp"

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
    T primal_value() {
        return _primal_value;
    }
    /**
     * @return a deep copy of the derivative value.
     * Note that this is shallow, and heap resources -- i.e. vector contents -- will be shared.
     */
    T deriv_value() {
        return _deriv_value;
    }

    /*
     * DualNumber–DualNumber arithmetic.
     * Derivative rules applied:
     *   Addition/subtraction : sum rule
     *   Multiplication       : product rule
     *   Division             : quotient rule
     */
    DualNumber operator+(const DualNumber& rhs) const {
        return { _primal_value + rhs._primal_value, _deriv_value + rhs._deriv_value };
    }
    DualNumber operator-(const DualNumber& rhs) const {
        return { _primal_value - rhs._primal_value, _deriv_value - rhs._deriv_value };
    }
    DualNumber operator*(const DualNumber& rhs) const {
        return { _primal_value * rhs._primal_value,
                 _deriv_value * rhs._primal_value + _primal_value * rhs._deriv_value };
    }
    DualNumber operator/(const DualNumber& rhs) const {
        return {
            _primal_value / rhs._primal_value,
            (_deriv_value * rhs._primal_value - _primal_value * rhs._deriv_value) / rhs._primal_value.pow(2u)
        };
    }

    /*
     * DualNumber–scalar arithmetic
     */
    DualNumber operator+(double scalar) const {
        return { _primal_value + scalar, _deriv_value };
    }
    DualNumber operator-(double scalar) const {
        return { _primal_value - scalar, _deriv_value };
    }
    DualNumber operator*(double scalar) const {
        return { _primal_value * scalar, _deriv_value * scalar };
    }
    DualNumber operator/(double scalar) const {
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
    DualNumber pow(uint32_t n) const {
        if (n == 0) {
            // f^0 = 1 (constant), derivative = 0
            return DualNumber(_primal_value.pow(0), _deriv_value * 0.0);
        }

        return { _primal_value.pow(n), _primal_value.pow(n - 1) * static_cast<double>(n) * _deriv_value };
    }
    /**
     * @return A new dual number representing the absolute value, d/dx(|f|) = sign(f) * f'.
     */
    DualNumber abs() const {
        if (_primal_value >= 0.0) {
            return { _primal_value.abs(), _deriv_value };
        } else if (_primal_value < 0.0) {
            return { _primal_value.abs(), _deriv_value * -1.0 };
        }
        // Mixed sign (e.g. an interval spanning zero): conservative bound.
        return { _primal_value.abs(), _deriv_value.abs() };
    }
    /**
     * @return A new dual number representing the result of exp with chain rule applied.
     */
    DualNumber exp() const {
        auto exp_primal = _primal_value.exp();
        return { exp_primal, exp_primal * _deriv_value };
    }

    /*
     * Compositional operations
     */

    /**
     * @param rhs Other DualNumber to union with.
     * @return A new DualNumber whose components are the respective unions.
     */
    DualNumber union_with(const DualNumber &rhs) const {
        return { _primal_value.union_with(rhs._primal_value),
                 _deriv_value.union_with(rhs._deriv_value) };
    }
    /**
     * @param n_splits Number of splits to produce.
     * @return Vector of n_splits DualNumbers.
     */
    std::vector<DualNumber> split(uint32_t n_splits) const {
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

#endif //DIFFDOMAIN_DUALNUMBER_H

