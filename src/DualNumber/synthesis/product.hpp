//
// Created by will on 7/6/26.
//

#ifndef DIFFDOMAIN_PRODUCT_HPP
#define DIFFDOMAIN_PRODUCT_HPP

inline Eigen::VectorXd product_rule(const Eigen::ArrayXd &f, const Eigen::ArrayXd &g,
                                          const Eigen::ArrayXd &f_prime, const Eigen::ArrayXd &g_prime) {
    return f_prime * g + f * g_prime;
}

inline Winterval compute_opt_interval_product(const Eigen::MatrixXd &corners) {
    auto product_rule_on_corners = corners.col(0).array() * corners.col(3).array() + corners.col(1).array() * corners.col(2).array();
    return { product_rule_on_corners.minCoeff(), product_rule_on_corners.maxCoeff() };
}

inline AffineForm compute_opt_affine_product(const DualNumber<MixedForm> &lhs, const DualNumber<MixedForm> &rhs, const Eigen::MatrixXd &corners) {
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
    auto sampled_derivatives = product_rule(f, g, f_prime, g_prime);
    // auto sampled_derivatives = product.col(0).array() * product.col(3).array() + product.col(1).array() * product.col(2).array();

    auto fit = regress_grid(xys, sampled_derivatives);
    auto intercept = fit(0);
    auto c1 = fit(1);
    auto c2 = fit(2);
    auto c3 = fit(3);
    auto c4 = fit(4);

    // compute maximum possible deviation.
    auto product_rule_on_corners = corners.col(0).array() * corners.col(3).array() + corners.col(1).array() * corners.col(2).array();

    auto evaluation = c1 * corners.col(0).array()
                                    + c2 * corners.col(1).array()
                                    + c3 * corners.col(2).array()
                                    + c4 * corners.col(3).array()
                                    + intercept;

    auto max_deviation = (product_rule_on_corners - evaluation).cwiseAbs().maxCoeff();
    auto affine_result = lhs.primal_ref().affine_rep() * c1 + rhs.primal_ref().affine_rep() * c2 + lhs.deriv_ref().affine_rep() * c3 + rhs.deriv_ref().affine_rep() * c4 + intercept;
    affine_result.add_noise_symbol(max_deviation);

    return affine_result;
}

#endif //DIFFDOMAIN_PRODUCT_HPP
