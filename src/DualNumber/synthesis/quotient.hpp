//
// Created by will on 7/6/26.
// Helpers for computing the quotient rule.
//

#ifndef DIFFDOMAIN_QUOTIENT_HPP
#define DIFFDOMAIN_QUOTIENT_HPP

inline Eigen::VectorXd quotient_rule(const Eigen::ArrayXd &f, const Eigen::ArrayXd &g,
                                          const Eigen::ArrayXd &f_prime, const Eigen::ArrayXd &g_prime) {
    return (f_prime * g - f * g_prime) / g.cwisePow(2);
}
inline double quotient_rule(double f, double g, double f_prime, double g_prime) {
    return (f_prime * g - f * g_prime) / (g * g);
}
// Using Pasado formula, check if the root is valid
inline bool valid_root(const Eigen::ArrayXd &root, double g_min, double g_max) {
    return (root[1] != 0)
        && (root[0] * root[3] != 0)
        && (2 * root[0] * root[3] / root[1] >= g_min)
        && (2 * root[0] * root[3] / root[1] <= g_max);
}

/**
 * Determine valid internal roots for the quotient rule synthesized linear abstract transformer.
 * @param primal_1_bounds bounds on f, (min, max)
 * @param deriv_1_bounds bounds on f'
 * @param primal_2_bounds bounds on g
 * @param deriv_2_bounds bounds on g'
 * @param g_linear_coeff synthesized coefficient for g
 * @return A list of candidate root points, which should be evaluated for soundness.
 */
inline std::vector<std::array<double, 4>> find_candidate_root_points(
    const std::array<double, 2> &primal_1_bounds,
    const std::array<double, 2> &deriv_1_bounds,
    const std::array<double, 2> &primal_2_bounds,
    const std::array<double, 2> &deriv_2_bounds,
    const double g_linear_coeff) {

    auto candidate_root_points = std::vector<std::array<double, 4>>();

    auto polynomial_coefficients = std::vector<double>(4);
    // first coefficient always coefficient of G
    polynomial_coefficients[0] = g_linear_coeff;
    polynomial_coefficients[1] = 0;

    for (auto f: primal_1_bounds) {
        for (auto f_prime: deriv_1_bounds) {
            for (auto g_prime: deriv_2_bounds) {
                polynomial_coefficients[2] = f_prime;
                polynomial_coefficients[3] = -2 * f * g_prime;
                auto roots = poly_roots(polynomial_coefficients);
                std::vector<double> g_roots_constrained;

                // find values of x2 that solve the equation, then add all the permutations of points as candidates to be checked.
                for (auto root: roots) {
                    // Don't want complex roots or roots outside of the boundaries
                    if (root.imag() == 0 && root.real() >= primal_2_bounds[0] && root.real() <= primal_2_bounds[1]) {
                        g_roots_constrained.push_back(root.real());
                    }
                }

                for (auto g_root: g_roots_constrained) {
                    candidate_root_points.push_back({f, g_root, f_prime, g_prime});
                }
            }
        }
    }

    return candidate_root_points;
}

/**
 * Compute the optimal affine bounds for the quotient rule.
 * @param corners Cartesian product of corners of primal, derivative values.
 * @param lhs Left operand, i.e. f
 * @param rhs Right operand, i.e. g
 * @return Optimal bounds for the affine quotient
 */
inline AffineForm compute_opt_affine_quotient(const Eigen::MatrixXd &corners, const DualNumber<MixedForm> &lhs, const DualNumber<MixedForm> &rhs) {
    constexpr uint32_t regression_number = 5u;

    auto primal_1_range = Eigen::VectorXd::LinSpaced(regression_number, lhs.primal_ref().min(), lhs.primal_ref().max());
    auto primal_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs.primal_ref().min(), rhs.primal_ref().max());
    auto deriv_1_range = Eigen::VectorXd::LinSpaced(regression_number, lhs.deriv_ref().min(), lhs.deriv_ref().max());
    auto deriv_2_range = Eigen::VectorXd::LinSpaced(regression_number, rhs.deriv_ref().min(), rhs.deriv_ref().max());
    auto xys = cartesian_product(primal_1_range, primal_2_range, deriv_1_range, deriv_2_range);

    auto f = xys.col(0).array();
    auto g = xys.col(1).array();
    auto f_prime = xys.col(2).array();
    auto g_prime = xys.col(3).array();
    auto sampled_derivatives = quotient_rule(f, g, f_prime, g_prime);

    // Synthesize a well-formed linear approximation of the quotient rule output.
    auto fit = regress_svd(xys, sampled_derivatives);
    auto intercept = fit(0);
    auto c1 = fit(1);
    auto c2 = fit(2);
    auto c3 = fit(3);
    auto c4 = fit(4);
    // To preserve soundness, remove any 0 coefficients.
    // Specifically, this ensures that f''(x) - A = 0.
    // If this does not hold, then the procedure will evaluate non-critical points.
    auto sig = 1e-8;
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

    auto max_deviation = (quotient_rule_on_corners - evaluation).abs().maxCoeff();

    // Now, to ensure soundness, we must iterate interior roots
    // Note that primal_1, deriv_1, and deriv_2 must be either lower or upper bounds
    // g^2 denominator is the remaining potential interior critical point
    auto candidate_root_points = find_candidate_root_points(
        {lhs.primal_ref().min(), lhs.primal_ref().max()},
        {lhs.deriv_ref().min(), lhs.deriv_ref().max()},
        {rhs.primal_ref().min(), rhs.primal_ref().max()},
        {rhs.deriv_ref().min(), rhs.deriv_ref().max()},
        c2);

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

    auto affine_result = rhs.primal_ref().affine_rep() * c1 + rhs.primal_ref().affine_rep() * c2 + lhs.deriv_ref().affine_rep() * c3 + rhs.deriv_ref().affine_rep() * c4 + intercept;
    affine_result.add_noise_symbol(max_deviation);
    return affine_result;
}

/**
 * Compute the optimal interval bounds for the quotient rule.
 * @param corners Corner points of the hypercube to evaluate for bounds optimization
 * @param g_min Minimum value of g
 * @param g_max Maximum value of g
 */
inline Winterval compute_opt_interval_quotient(const Eigen::MatrixXd &corners,
                                               const double g_min, const double g_max) {
    // Evaluate quotient rule on all corner points
    std::vector<double> corner_evals;
    corner_evals.reserve(corners.rows());
    for (int i = 0; i < corners.rows(); ++i) {
        corner_evals.push_back(quotient_rule(corners(i, 0), corners(i, 1), corners(i, 2), corners(i, 3)));
    }

    // Evaluate quotient rule on valid internal critical points
    std::vector<double> root_evals;
    for (int i = 0; i < corners.rows(); ++i) {
        auto row = corners.row(i).array();
        if (valid_root(row, g_min, g_max)) {
            auto f = row[0];
            auto f_prime = row[2];
            auto g_prime = row[3];
            double critical_g = (2 * f * g_prime) / f_prime;
            root_evals.push_back(quotient_rule(f, critical_g, f_prime, g_prime));
        }
    }

    // Combine all evaluations
    if (root_evals.empty()) {
        return { *std::ranges::min_element(corner_evals), *std::ranges::max_element(corner_evals) };
    }
    // corner evaluations cannot be empty
    auto lower = *std::min(std::ranges::min_element(corner_evals), std::ranges::min_element(root_evals));
    auto upper = *std::max(std::ranges::max_element(corner_evals), std::ranges::max_element(root_evals));
    return {lower, upper};
}

#endif //DIFFDOMAIN_QUOTIENT_HPP