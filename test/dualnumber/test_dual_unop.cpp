//
// Created by will on 3/12/26.
//

#include <cmath>
#include <gtest/gtest.h>

#include "DualNumber/DualNumber.hpp"
#include "Real/Real.hpp"
#include "Winterval/Winterval.hpp"

// x = 3, dx/dx = 1  (the independent variable)
static DualNumber<Winterval> makeX() { return { Winterval(3, 3), Winterval(1, 1) }; }

// ── pow ───────────────────────────────────────────────────────────────────────

TEST(dual_unop, pow_zero) {
    // f^0 = 1 (constant), derivative = 0
    EXPECT_EQ(makeX().pow(0), DualNumber(Winterval(1, 1), Winterval(0, 0)));
}
TEST(dual_unop, pow_zero_scalar) {
    EXPECT_EQ(3, DualNumber(Real(1), Real(1)).pow(3).deriv_value());
}

TEST(dual_unop, pow_one) {
    // d/dx(x^1) = 1
    EXPECT_EQ(makeX().pow(1), DualNumber(Winterval(3, 3), Winterval(1, 1)));
}

TEST(dual_unop, pow_two) {
    // d/dx(x^2) = 2x = 6 at x = 3
    EXPECT_EQ(makeX().pow(2), DualNumber(Winterval(9, 9), Winterval(6, 6)));
}

TEST(dual_unop, pow_three) {
    // d/dx(x^3) = 3x^2 = 27 at x = 3
    EXPECT_EQ(makeX().pow(3), DualNumber(Winterval(27, 27), Winterval(27, 27)));
}

// ── abs ───────────────────────────────────────────────────────────────────────

TEST(dual_unop, abs_positive) {
    // Positive primal: sign = +1, derivative unchanged.
    DualNumber x(Winterval(3, 3), Winterval(2, 2));
    EXPECT_EQ(x.abs(), DualNumber(Winterval(3, 3), Winterval(2, 2)));
}

TEST(dual_unop, abs_negative) {
    // Negative primal: sign = -1, derivative negated.
    DualNumber x(Winterval(-3, -3), Winterval(2, 2));
    EXPECT_EQ(x.abs(), DualNumber(Winterval(3, 3), Winterval(-2, -2)));
}

TEST(dual_unop, abs_zero) {
    // Zero is non-negative: derivative unchanged (valid subgradient).
    DualNumber x(Winterval(0, 0), Winterval(2, 2));
    EXPECT_EQ(x.abs(), DualNumber(Winterval(0, 0), Winterval(2, 2)));
}

TEST(dual_unop, abs_mixed_sign_interval) {
    // Interval spans zero: conservative bound – abs of derivative returned.
    DualNumber x(Winterval(-2, 3), Winterval(1, 1));
    auto result = x.abs();
    EXPECT_EQ(result, DualNumber(Winterval(0, 3), Winterval(1, 1)));
}

TEST(dual_unop, exp) {
    DualNumber x(Winterval(0, 0), Winterval(1, 1));
    auto result = x.exp();
    EXPECT_NEAR(result.primal_ref().min(), 1.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 1.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().min(), 1.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 1.0, 1e-6);
}

TEST(dual_unop, tanh) {
    DualNumber x(Winterval(1, 3), Winterval(1, 1));
    auto result = x.tanh().pow(2);
    ASSERT_NEAR(result.primal_ref().min(), 0.580026, 1e-6);
    ASSERT_NEAR(result.primal_ref().max(), 0.990134, 1e-6);
    ASSERT_NEAR(result.deriv_ref().min(), 0.01502783249528329, 1e-6);
    ASSERT_NEAR(result.deriv_ref().max(), 0.83579493009898309, 1e-6);
}

TEST(dual_unop, sigmoid) {
    DualNumber x(Winterval(1, 3), Winterval(1, 1));
    auto result = x.sigmoid();
    ASSERT_NEAR(result.primal_ref().min(), 0.731059, 1e-6);
    ASSERT_NEAR(result.primal_ref().max(), 0.952574, 1e-6);
    ASSERT_NEAR(result.deriv_ref().min(), 0.034671, 1e-6);
    ASSERT_NEAR(result.deriv_ref().max(), 0.256187, 1e-6);
}

TEST(dual_unop, relu) {
    DualNumber x(Winterval(2, 3), Winterval(1, 1));
    auto result = x.relu();
    EXPECT_EQ(result.primal_ref().min(), 2.0);
    EXPECT_EQ(result.primal_ref().max(), 3.0);
    EXPECT_EQ(result.deriv_ref().min(), 1.0);
    EXPECT_EQ(result.deriv_ref().max(), 1.0);

    x = DualNumber(Winterval(-3, -1), Winterval(2, 2));
    result = x.relu();
    EXPECT_EQ(result.primal_ref().min(), 0.0);
    EXPECT_EQ(result.primal_ref().max(), 0.0);
    EXPECT_EQ(result.deriv_ref().min(), 0.0);
    EXPECT_EQ(result.deriv_ref().max(), 0.0);

    x = DualNumber(Winterval(-2, 3), Winterval(1, 2));
    result = x.relu();
    EXPECT_EQ(result.primal_ref().min(), 0.0);
    EXPECT_EQ(result.primal_ref().max(), 3.0);
    EXPECT_EQ(result.deriv_ref().min(), 1.0);
    EXPECT_EQ(result.deriv_ref().max(), 2.0);
}

TEST(dual_unop, sqrt_point_value) {
    // d/dx(sqrt(x)) = 1/(2*sqrt(x))
    // At x = 4: sqrt(4) = 2, derivative = 1/(2*2) = 0.25
    DualNumber x(Winterval(4, 4), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_NEAR(result.primal_ref().min(), 2.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 2.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().min(), 0.25, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 0.25, 1e-6);
}

TEST(dual_unop, sqrt_zero) {
    // At x = 0: sqrt(0) = 0, derivative = 1/(2*sqrt(0)) → infinity
    DualNumber x(Winterval(0, 0), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_EQ(result.primal_ref().min(), 0.0);
    EXPECT_EQ(result.primal_ref().max(), 0.0);
    EXPECT_TRUE(std::isnan(result.deriv_ref().min()));
    EXPECT_TRUE(std::isnan(result.deriv_ref().max()));
}

TEST(dual_unop, sqrt_one) {
    // At x = 1: sqrt(1) = 1, derivative = 1/(2*1) = 0.5
    DualNumber x(Winterval(1, 1), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_NEAR(result.primal_ref().min(), 1.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 1.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().min(), 0.5, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 0.5, 1e-6);
}

TEST(dual_unop, sqrt_interval_1_to_4) {
    // sqrt([1, 4]) = [1, 2]
    // d/dx(sqrt(x)) = 1/(2*sqrt(x))
    // At x=1: deriv = 0.5, at x=4: deriv = 0.25
    DualNumber x(Winterval(1, 4), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_NEAR(result.primal_ref().min(), 1.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 2.0, 1e-6);
    // Derivative is decreasing: 1/(2*sqrt(x)) is smaller for larger x
    EXPECT_NEAR(result.deriv_ref().min(), 0.25, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 0.5, 1e-6);
}

TEST(dual_unop, sqrt_interval_9_to_25) {
    // sqrt([9, 25]) = [3, 5]
    // At x=9: deriv = 1/(2*3) = 1/6 ≈ 0.1667
    // At x=25: deriv = 1/(2*5) = 0.1
    DualNumber x(Winterval(9, 25), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_NEAR(result.primal_ref().min(), 3.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 5.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().min(), 0.1, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 1.0/6.0, 1e-6);
}

TEST(dual_unop, sqrt_with_scaled_derivative) {
    // Test with non-unit derivative: dx/dy where y is inner function
    // sqrt(x) with dx/dy = 2
    DualNumber x(Winterval(4, 4), Winterval(2, 2));
    auto result = x.sqrt();
    // sqrt(4) = 2
    // d/dx(sqrt(x)) = 1/(2*sqrt(x)) = 0.25
    // Chain rule: 0.25 * 2 = 0.5
    EXPECT_NEAR(result.primal_ref().min(), 2.0, 1e-6);
    EXPECT_NEAR(result.primal_ref().max(), 2.0, 1e-6);
    EXPECT_NEAR(result.deriv_ref().min(), 0.5, 1e-6);
    EXPECT_NEAR(result.deriv_ref().max(), 0.5, 1e-6);
}

TEST(dual_unop, sqrt_primal_bounds_are_correct) {
    // Verify that primal bounds satisfy sqrt property
    DualNumber x(Winterval(0, 16), Winterval(1, 1));
    auto result = x.sqrt();
    
    // sqrt(0) = 0, sqrt(16) = 4
    EXPECT_LE(result.primal_ref().min(), 0.0);
    EXPECT_GE(result.primal_ref().max(), 4.0);
    
    // Verify intermediate values are bounded
    for (double val = 0; val <= 16; val += 1.0) {
        double sqrt_val = std::sqrt(val);
        ASSERT_LE(result.primal_ref().min(), sqrt_val);
        ASSERT_GE(result.primal_ref().max(), sqrt_val);
    }
}

TEST(dual_unop, sqrt_negative_interval) {
    // sqrt of negative interval produces NaN
    DualNumber x(Winterval(-4, -1), Winterval(1, 1));
    auto result = x.sqrt();
    EXPECT_TRUE(std::isnan(result.primal_ref().min()));
    EXPECT_TRUE(std::isnan(result.primal_ref().max()));
}