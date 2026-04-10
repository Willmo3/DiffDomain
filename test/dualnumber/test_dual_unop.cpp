//
// Created by will on 3/12/26.
//

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
