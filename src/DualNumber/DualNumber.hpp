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

#include "cereal/cereal.hpp"
#include "Numeric.hpp"
#include "Eigen/Dense"
#include "../shared/polynomial/regression.hpp"

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
#include "shared/polynomial/rootfinding.hpp"

// synthesized abstract transformers for autodiff derived from Pasado
// https://dl.acm.org/doi/pdf/10.1145/3622867

inline Eigen::VectorXd quotient_rule(const Eigen::ArrayXd &f, const Eigen::ArrayXd &g,
                                          const Eigen::ArrayXd &f_prime, const Eigen::ArrayXd &g_prime) {
    return (f_prime * g - f * g_prime) / g.cwisePow(2);
}
inline double quotient_rule(double f, double g, double f_prime, double g_prime) {
    return (f_prime * g - f * g_prime) / (g * g);
}
inline Eigen::VectorXd product_rule(const Eigen::ArrayXd &f, const Eigen::ArrayXd &g,
                                          const Eigen::ArrayXd &f_prime, const Eigen::ArrayXd &g_prime) {
    return f_prime * g + f * g_prime;
}

/**
 * Compute the optimal interval bounds for the quotient rule.
 * @param corners Corner points of the hypercube to evaluate for bounds optimization
 * @param candidate_roots Internal critical roots
 * @parma
 */
inline Winterval compute_opt_interval_quotient(const Eigen::MatrixXd &corners,
                                               const std::vector<std::array<double, 4> > &candidate_roots,
                                               double g_min, double g_max) {

    // Evaluate quotient rule on all corner points
    std::vector<double> corner_evals;
    corner_evals.reserve(corners.rows());
    for (int i = 0; i < corners.rows(); ++i) {
        corner_evals.push_back(quotient_rule(corners(i, 0), corners(i, 1), corners(i, 2), corners(i, 3)));
    }

    // Evaluate quotient rule on valid internal critical points

    auto valid_root = [g_min, g_max](const std::array<double, 4> &a) -> bool {
        return (a[1] != 0) && ((a[0] * a[3]) != 0) && (((2 * a[0] * a[3]) / a[1]) >= g_min) && (
                   ((2 * a[0] * a[3]) / a[1]) <= g_max);
    };

    std::vector<double> root_evals;
    for (const auto &point: candidate_roots) {
        if (valid_root(point)) {
            auto f = point[0];
            auto f_prime = point[2];
            auto g_prime = point[3];
            double critical_g = (2 * f * g_prime) / f_prime;
            root_evals.push_back(quotient_rule(f, critical_g, f_prime, g_prime));
        }
    }

    // Combine all evaluations
    std::vector<double> all_evals;
    all_evals.reserve(corner_evals.size() + root_evals.size());
    all_evals.insert(all_evals.end(), corner_evals.begin(), corner_evals.end());
    all_evals.insert(all_evals.end(), root_evals.begin(), root_evals.end());

    // Find min and max
    double lower = *std::ranges::min_element(all_evals);
    double upper = *std::ranges::max_element(all_evals);

    return {lower, upper};
}

/**
 * Synthesized abstract transformer for autodiff of product rule.
 * We will use standard operations for the primal, but abstract the derivative.
 */
template<>
[[nodiscard]] inline DualNumber<MixedForm> DualNumber<MixedForm>::operator*(const DualNumber&rhs) const {
    constexpr uint32_t regression_number = 5u;

    auto primal_1_range = Eigen::VectorXd::LinSpaced(regression_number, _primal_value.min(), _primal_value.max());
    auto primal_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs._primal_value.min(), rhs._primal_value.max());
    auto deriv_1_range = Eigen::VectorXd::LinSpaced(regression_number, _deriv_value.min(), _deriv_value.max());
    auto deriv_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs._deriv_value.min(), rhs._deriv_value.max());
    auto xys = cartesian_product(primal_1_range, primal_2_range, deriv_1_range, deriv_2_range);

    auto f = xys.col(0).array();
    auto g = xys.col(1).array();
    auto f_prime = xys.col(2).array();
    auto g_prime = xys.col(3).array();
    auto sampled_derivatives = product_rule(f, g, f_prime, g_prime);
    // auto sampled_derivatives = product.col(0).array() * product.col(3).array() + product.col(1).array() * product.col(2).array();

    auto fit = regress(xys, sampled_derivatives);
    auto intercept = fit(0);
    auto c1 = fit(1);
    auto c2 = fit(2);
    auto c3 = fit(3);
    auto c4 = fit(4);

    // compute maximum possible deviation.
    auto corners = get_corners(_primal_value.interval_bounds(), rhs._primal_value.interval_bounds(), _deriv_value.interval_bounds(), rhs._deriv_value.interval_bounds());
    auto product_rule_on_corners = corners.col(0).array() * corners.col(3).array() + corners.col(1).array() * corners.col(2).array();

    auto evaluation = c1 * corners.col(0).array()
                                    + c2 * corners.col(1).array()
                                    + c3 * corners.col(2).array()
                                    + c4 * corners.col(3).array()
                                    + intercept;

    auto max_deviation = (product_rule_on_corners - evaluation).cwiseAbs().maxCoeff();
    auto affine_result = _primal_value.affine_rep() * c1 + rhs._primal_value.affine_rep() * c2 + _deriv_value.affine_rep() * c3 + rhs._deriv_value.affine_rep() * c4 + intercept;
    affine_result.add_noise_symbol(max_deviation);

    auto max_product_rule = product_rule_on_corners.maxCoeff();
    auto min_product_rule = product_rule_on_corners.minCoeff();

    auto derivative = MixedForm(affine_result, Winterval(min_product_rule, max_product_rule));

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
    constexpr uint32_t regression_number = 5u;

    // Default to regular division if we have a range w/ negative values.
    if (rhs.min() <= 0 && rhs.max() >= 0) {
        return {
            _primal_value / rhs._primal_value,
            (_deriv_value * rhs._primal_value - _primal_value * rhs._deriv_value) / rhs._primal_value.pow(2u)
        };
    }

    auto primal_1_range = Eigen::VectorXd::LinSpaced(regression_number, _primal_value.min(), _primal_value.max());
    auto primal_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs._primal_value.min(), rhs._primal_value.max());
    auto deriv_1_range = Eigen::VectorXd::LinSpaced(regression_number, _deriv_value.min(), _deriv_value.max());
    auto deriv_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs._deriv_value.min(), rhs._deriv_value.max());
    auto xys = cartesian_product(primal_1_range, primal_2_range, deriv_1_range, deriv_2_range);
    
    auto f = xys.col(0).array();
    auto g = xys.col(1).array();
    auto f_prime = xys.col(2).array();
    auto g_prime = xys.col(3).array();
    auto sampled_derivatives = quotient_rule(f, g, f_prime, g_prime);

    // Synthesize a well-formed linear approximation of the quotient rule output.
    auto fit = regress(xys, sampled_derivatives);
    auto intercept = fit(0);
    auto c1 = fit(1);
    auto c2 = fit(2);
    auto c3 = fit(3);
    auto c4 = fit(4);
    // To preserve soundness, remove any 0 coefficients.
    auto sig = 0.00001;
    if (c1 == 0) {
        c1 += sig;
    }
    if (c2 == 0) {
        c2 += sig;
    }
    if (c3 == 0) {
        c3 += sig;
    }
    if (c4 == 0) {
        c4 += sig;
    }

    // compute maximum possible deviation.
    auto corners = get_corners(_primal_value.interval_bounds(), rhs._primal_value.interval_bounds(), _deriv_value.interval_bounds(), rhs._deriv_value.interval_bounds());

    auto corner_f = corners.col(0).array();
    auto corner_g = corners.col(1).array();
    auto corner_f_prime = corners.col(2).array();
    auto corner_g_prime = corners.col(3).array();
    auto quotient_rule_on_corners = quotient_rule(corner_f, corner_g, corner_f_prime, corner_g_prime).array();

    auto evaluation = c1 * corner_f
                                    + c2 * corner_g
                                    + c3 * corner_f_prime
                                    + c4 * corner_g_prime
                                    + intercept;

    auto max_deviation = (quotient_rule_on_corners - evaluation).cwiseAbs().maxCoeff();

    // Now, to ensure soundness, we must iterate interior roots
    // Note that primal_1, deriv_1, and deriv_2 must be either lower or upper bounds
    // g^2 denominator is the remaining potential interior critical point
    auto candidate_root_points = std::vector<std::array<double, 4>>();

    auto polynomial_coefficients = std::vector<double>(4);
    // first coefficient always coefficient of G
    polynomial_coefficients[0] = c2;
    polynomial_coefficients[1] = 0;

    for (auto f: primal_1_range) {
        for (auto f_prime: deriv_1_range) {
            for (auto g_prime: deriv_2_range) {
                polynomial_coefficients[2] = f_prime;
                polynomial_coefficients[3] = -2 * f * g_prime;
                auto roots = poly_roots(polynomial_coefficients);
                std::vector<double> g_roots_constrained;

                // find values of x2 that solve the equation, then add all the permutations of points as candidates to be checked.
                for (auto root: roots) {
                    // Don't want complex roots or roots outside of the boundaries
                    if (root.imag() == 0 && root.real() >= rhs.primal_ref().min() && root.real() <= rhs.primal_ref().max()) {
                        g_roots_constrained.push_back(root.real());
                    }
                }

                for (auto g_root: g_roots_constrained) {
                    candidate_root_points.push_back({f, g_root, f_prime, g_prime});
                }
            }
        }
    }
    // Check if any candidate roots are the new max.
    for (auto points : candidate_root_points) {
        auto f = points[0];
        auto g = points[1];
        auto f_prime = points[2];
        auto g_prime = points[3];

        auto difference = std::abs((g * f_prime - f * g_prime) / std::pow(g, 2u));
        if (difference > max_deviation) {
            max_deviation = difference;
        }
    }

    auto affine_result = _primal_value.affine_rep() * c1 + rhs._primal_value.affine_rep() * c2 + _deriv_value.affine_rep() * c3 + rhs._deriv_value.affine_rep() * c4 + intercept;
    affine_result.add_noise_symbol(max_deviation);

    // Evaluating extrema for interval, we must be careful to remove any extrema that violate the equation mentioned in Pasado.
    auto lx2 = rhs._primal_value.min();
    auto ux2 = rhs._primal_value.max();

    auto refined_interval = compute_opt_interval_quotient(corners, candidate_root_points, lx2, ux2);

    auto derivative = MixedForm(affine_result, refined_interval);

    return {
        _primal_value / rhs._primal_value,
        derivative
    };
}

#endif

#endif //DIFFDOMAIN_DUALNUMBER_H
