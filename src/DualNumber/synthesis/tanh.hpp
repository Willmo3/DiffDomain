//
// Created by will on 7/6/26.
// Various routines for synthesizing derivatives of dual numbers.
//

#ifndef DIFFDOMAIN_TANH_HPP
#define DIFFDOMAIN_TANH_HPP
#include <cmath>
#include <complex>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

inline double arctanh(double y) {
    return 0.5 * log((1 + y) / (1 - y));
}

/**
 * Helper function to compute the complex cube root.
 * @param z Complex number.
 * @return Cube root of z.
 */
inline std::complex<double> cuberoot(const std::complex<double> &z) {
    return std::pow(z, 1.0 / 3.0);
}

/**
 * Vectorized chain rule.
 * @param primal Real value of the full function.
 * @param derivative Derivative of the inner function.
 * @return chain rule applied to each element of the operands.
 */
inline Eigen::VectorXd tanh_chain(const Eigen::VectorXd &primal, const Eigen::VectorXd &derivative) {
    auto tanh_deriv = Eigen::VectorXd::Ones(primal.size()).array() - primal.array().tanh() * primal.array().tanh();
    return tanh_deriv * derivative.array();
}
/**
 * Scalar chain rule.
 * @param primal Real value of the full function.
 * @param derivative Derivative of the inner function.
 * @return chain rule applied to operands
 */
inline double tanh_chain(double primal, double derivative) {
    auto tanh_deriv = 1.0 - std::tanh(primal) * std::tanh(primal);
    return tanh_deriv * derivative;
}

/**
 * Compute the roots of the synthesized transformer optimization problem using Cardano's formula.
 *
 * @param y real value to compute roots of.
 * @return Array containing the three roots (x1, x2, x3) as complex numbers.
 */
inline std::array<std::complex<double>, 3>
synthesized_transformer_roots_tanh(double y) {
    // Compute Cardano
    auto Q = std::complex(-1.0 / 3.0);
    auto R = std::complex(0.25 * y);
    auto sqrt_term = std::sqrt(std::complex(std::pow(Q, 3.0) + std::pow(R, 2.0)));
    auto S = cuberoot(R + sqrt_term);
    auto T = cuberoot(R - sqrt_term);

    // From substituted S and T retrive actual roots
    auto x1 = S + T;
    const auto i_sqrt3_2 = std::complex(0.0, std::sqrt(3.0) / 2.0);
    auto x2 = -0.5 * (S + T) + i_sqrt3_2 * (S - T);
    auto x3 = -0.5 * (S + T) - i_sqrt3_2 * (S - T);

    return { x1, x2, x3 };
}

/**
 * Compute the possible interior critical points of the tanh bounds synthesis problem defined in Pasado.
 * @param synthesized_coefficient Coefficient of the synthesized polynomial.
 * @param g_bound Bound for g(x) being considered -- lower or upper.
 * @param min_f Minimum value for f(x) -- will exclude exterior roots
 * @param max_f Maximum value for f(x) -- will exclude exterior roots
 * @return Vector of real roots from the inverse polynomial.
 */
inline std::vector<double> tanh_interior_critical_points(double synthesized_coefficient, double g_bound, double min_f, double max_f) {
    if (synthesized_coefficient == 0) {
        // 0-valued degenerate case: critical point is at exactly x=0
        // However, this may not be a valid critical point if 0 is outside of the boundaries of f.
        if (0 >= min_f && 0 <= max_f) {
            return {0};
        }
        return {};
    }

    auto val = synthesized_coefficient / g_bound;
    auto found_roots = synthesized_transformer_roots_tanh(val);
    std::vector<double> valid_roots;

    // imaginary
    auto tolerance = 1e-8;
    for (auto &root : found_roots) {
        // Roots are values of tanh(x): map back to x before checking against f's bounds.
        auto real = root.real();
        if (std::abs(root.imag()) >= tolerance || std::abs(real) >= 1) {
            continue;
        }
        auto x = arctanh(real);
        if (x >= min_f && x <= max_f) {
            valid_roots.push_back(x);
        }
    }

    return valid_roots;
}

inline AffineForm compute_opt_affine_tanh(const DualNumber<MixedForm> &value) {
    constexpr uint32_t sampled_points_per_dim = 8;

    auto primal_lower = value.primal_ref().min();
    auto primal_upper = value.primal_ref().max();
    auto deriv_lower = value.deriv_ref().min();
    auto deriv_upper = value.deriv_ref().max();

    auto primal_range = Eigen::VectorXd::LinSpaced(sampled_points_per_dim, primal_lower, primal_upper);
    auto deriv_range = Eigen::VectorXd::LinSpaced(sampled_points_per_dim, deriv_lower, deriv_upper);
    auto joint_range = cartesian_product(primal_range, deriv_range);

    auto primal_values = joint_range.col(0).array();
    auto deriv_values = joint_range.col(1).array();
    auto sampled_outputs = tanh_chain(primal_values, deriv_values);

    auto fit = regress_grid(joint_range, sampled_outputs);
    auto intercept = fit[0];
    auto c1 = fit[1];
    auto c2 = fit[2];
    // Ensure first coefficient is non-zero.
    constexpr double sigma = 1e-8;
    if (c1 == 0.0) {
        c1 += sigma;
    }

    // Ensure soundness by evaluating function at all possible extrema.
    // These are found by applying the cubic formula to the bounds optimization problem
    // The corners are additional candidate critical points
    auto candidate_extrema_lower = tanh_interior_critical_points(c1, deriv_lower, primal_lower, primal_upper);
    auto candidate_extrema_upper = tanh_interior_critical_points(c1, deriv_upper, primal_lower, primal_upper);

    std::vector primal_eval_points = {primal_lower, primal_upper};
    primal_eval_points.insert(primal_eval_points.end(), candidate_extrema_lower.begin(), candidate_extrema_lower.end());
    primal_eval_points.insert(primal_eval_points.end(), candidate_extrema_upper.begin(), candidate_extrema_upper.end());
    std::vector deriv_eval_points = {deriv_lower, deriv_upper};

    double max_deviation = 0;
    for (auto primal_val : primal_eval_points) {
        for (auto deriv_val : deriv_eval_points) {
            auto sampled_output = tanh_chain(primal_val, deriv_val);
            auto linear_approx = c1 * primal_val + c2 * deriv_val + intercept;
            auto deviation = std::abs(sampled_output - linear_approx);

            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }
    }

    auto result = value.primal_ref().affine_rep() * c1 + value.deriv_ref().affine_rep() * c2 + intercept;
    result.add_noise_symbol(max_deviation);
    return result;
}

inline Winterval compute_opt_interval_tanh(const DualNumber<MixedForm> &value) {
    const auto &primal = value.primal_ref();
    const auto &deriv = value.deriv_ref();

    std::vector corner_evals = {
        tanh_chain(primal.min(), deriv.min()),
        tanh_chain(primal.min(), deriv.max()),
        tanh_chain(primal.max(), deriv.min()),
        tanh_chain(primal.max(), deriv.max())
    };

    auto lower_extrema = tanh_interior_critical_points(0, deriv.min(), primal.min(), primal.max());
    auto upper_extrema = tanh_interior_critical_points(0, deriv.max(), primal.min(), primal.max());

    for (auto candidate: lower_extrema) {
        corner_evals.push_back(tanh_chain(candidate, deriv.min()));
    }
    for (auto candidate : upper_extrema) {
        corner_evals.push_back(tanh_chain(candidate, deriv.max()));
    }

    return {
        *std::ranges::min_element(corner_evals),
        *std::ranges::max_element(corner_evals)
    };
}


#endif //DIFFDOMAIN_TANH_HPP
