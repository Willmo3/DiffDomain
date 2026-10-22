//
// Created by will on 9/30/25.
//

// NOTE: gtest keeps static state in a given test suite (i.e. affine).
// Therefore, we switch noise symbols.

#include <cmath>
#include <gtest/gtest.h>

#include "Caffeine/AffineForm.hpp"

TEST(affine_unop, negation) {
    auto value = -AffineForm(Winterval(-2, 3));
    ASSERT_NEAR(value.center(), -0.5, 0.001);
    ASSERT_NEAR(value.radius(), 2.5, 0.001);
}

TEST(affine_unop, abs) {
    auto value = AffineForm(Winterval(1, 2)).abs();
    ASSERT_NEAR(value.center(), 1.5, 0.001);
    ASSERT_NEAR(value.radius(), 0.5, 0.001);

    value = AffineForm(Winterval(-2, -1)).abs();
    ASSERT_NEAR(value.center(), 1.5, 0.001);
    ASSERT_NEAR(value.radius(), 0.5, 0.001);

    // Straddles zero: must contain [0, 2]. Computed as 2 * relu(x) - x.
    value = AffineForm(Winterval(-1, 3)).abs();
    ASSERT_LE(value.min(), 0.0);
    ASSERT_GE(value.max(), 3.0);
    ASSERT_NEAR(value.center(), 1.25, 0.001);
    ASSERT_NEAR(value.radius(), 1.75, 0.001);
}

TEST(affine_unop, pow) {
    auto value = AffineForm(Winterval(1, 2)).pow(3);
    // Range [-0.875, 8] (was [-1.25, 8] before pow squared with a non-negative quadratic term).
    ASSERT_NEAR(value.center(), 3.5625, 0.001);
    ASSERT_NEAR(value.radius(), 4.4375, 0.001);
}

TEST(affine_unop, pow_zero) {
    auto value = AffineForm(Winterval(-5, 3)).pow(0);
    ASSERT_NEAR(value.center(), 1.0, 0.001);
    ASSERT_NEAR(value.radius(), 0.0, 0.001);
}

TEST(affine_unop, pow_one_identity) {
    auto a = AffineForm(Winterval(1, 4));
    ASSERT_EQ(a.pow(1).to_interval(), a.to_interval());
}

TEST(affine_unop, pow_two_soundness) {
    auto result = AffineForm(Winterval(1, 3)).pow(2);
    ASSERT_LE(result.min(), 1.0);
    ASSERT_GE(result.max(), 9.0);
}

TEST(affine_unop, pow_two_straddling_zero_non_negative) {
    // x * x would give [-1, 1]; squaring knows the result is non-negative.
    auto result = AffineForm(Winterval(-1, 1)).pow(2);
    ASSERT_DOUBLE_EQ(result.min(), 0.0);
    ASSERT_DOUBLE_EQ(result.max(), 1.0);
    ASSERT_FALSE(std::isnan(result.sqrt().center()));
}

TEST(affine_unop, pow_two_positive_interval_soundness) {
    auto result = AffineForm(Winterval(2, 4)).pow(2);
    ASSERT_LE(result.min(), 4.0);
    ASSERT_GE(result.max(), 16.0);
}

TEST(affine_unop, exp) {
    auto result = AffineForm(Winterval(1, 2)).exp();
    ASSERT_NEAR(result.min(), 2.142368, 1e-5);
    ASSERT_NEAR(result.max(), std::exp(2.0), 1e-5);
}

TEST(affine_unop, exp_soundness) {
    auto result = AffineForm(Winterval(0, 2)).exp();
    ASSERT_LE(result.min(), 1.0);
    ASSERT_GE(result.max(), std::exp(2.0));
}

TEST(affine_unop, exp_exact) {
    auto result = AffineForm(1.0).exp();
    ASSERT_DOUBLE_EQ(result.min(), std::exp(1.0));
    ASSERT_DOUBLE_EQ(result.max(), std::exp(1.0));
}

TEST(affine_unop, exp_unbounded) {
    auto result = (AffineForm(1.0) / AffineForm(Winterval(-1, 1))).exp();
    ASSERT_EQ(result.min(), -INFINITY);
    ASSERT_EQ(result.max(), INFINITY);
}

TEST(affine_unop, tanh) {
    auto result = AffineForm(Winterval(0.1, 0.2));
    result = result.tanh();
    ASSERT_NEAR(result.min(), 0.099668, 1e-5);
    ASSERT_NEAR(result.max(), 0.197375, 1e-5);
}

TEST(affine_unop, sigmoid) {
    auto result = AffineForm(Winterval(0.1, 0.2)).sigmoid();
    ASSERT_NEAR(result.min(), 0.524979, 1e-5);
    ASSERT_NEAR(result.max(), 0.549834, 1e-5);

    auto unit = AffineForm(Winterval(0.1, 0.1)).sigmoid();
    ASSERT_NEAR(unit.min(), 0.524979, 1e-5);
    ASSERT_NEAR(unit.max(), 0.524979, 1e-5);
}

TEST(affine_unop, relu) {
    auto result = AffineForm(Winterval(-100, 2)).relu();
    ASSERT_NEAR(result.min(), -1.960784, 1e-5);
    ASSERT_EQ(result.max(), 2);
}

TEST(affine_unop, relu_unbounded) {
    auto result = (AffineForm(1.0) / AffineForm(Winterval(-1, 1))).relu();
    ASSERT_EQ(result.min(), -INFINITY);
    ASSERT_EQ(result.max(), INFINITY);
}

TEST(affine_unop, relu_huge_bounds_no_overflow) {
    // u - l and u * l overflow here; the DeepZ terms used to come out NaN.
    // Radius 1.6e308 is finite, but u - l = 3.2e308 overflows.
    auto x = AffineForm(Winterval(-0.8e308, 0.8e308)) + AffineForm(Winterval(-0.8e308, 0.8e308));
    auto result = x.relu();
    ASSERT_LE(result.min(), 0.0);
    ASSERT_GE(result.max(), 1.6e308);
    ASSERT_TRUE(std::isfinite(result.min()) && std::isfinite(result.max()));
}

TEST(affine_unop, abs_huge_bounds_no_overflow) {
    auto x = AffineForm(Winterval(-0.8e308, 0.8e308)) + AffineForm(Winterval(-0.8e308, 0.8e308));
    auto result = x.abs();
    ASSERT_LE(result.min(), 0.0);
    ASSERT_GE(result.max(), 1.6e308);
    ASSERT_TRUE(std::isfinite(result.min()) && std::isfinite(result.max()));
}

TEST(affine_unop, sqrt_point_interval) {
    auto result = AffineForm(Winterval(4.0)).sqrt();
    ASSERT_NEAR(result.min(), 2, 1e-5);
    ASSERT_NEAR(result.max(), 2, 1e-5);
}

TEST(affine_unop, sqrt_zero) {
    auto result = AffineForm(Winterval(0.0)).sqrt();
    ASSERT_NEAR(result.min(), 0.0, 1e-5);
    ASSERT_NEAR(result.max(), 0.0, 1e-5);
}

TEST(affine_unop, sqrt_small_positive_interval) {
    auto result = AffineForm(Winterval(1.0, 4.0)).sqrt();
    ASSERT_NEAR(result.center(), 1.5, 0.1);
    ASSERT_GT(result.radius(), 0.0);
}

TEST(affine_unop, sqrt_soundness_1_to_4) {
    auto result = AffineForm(Winterval(1.0, 4.0)).sqrt();
    ASSERT_LE(result.min(), std::sqrt(1.0));
    ASSERT_GE(result.max(), std::sqrt(4.0));
}

TEST(affine_unop, sqrt_soundness_9_to_25) {
    auto result = AffineForm(Winterval(9.0, 25.0)).sqrt();
    ASSERT_LE(result.min(), 3.0);
    ASSERT_GE(result.max(), 5.0);
}

TEST(affine_unop, sqrt_soundness_0_to_1) {
    auto result = AffineForm(Winterval(0, 1.0)).sqrt();
    ASSERT_LE(result.min(), 0.0);
    ASSERT_GE(result.max(), 1.0);
}

TEST(affine_unop, sqrt_negative_interval) {
    auto result = AffineForm(Winterval(-4.0, -1.0)).sqrt();
    ASSERT_TRUE(std::isnan(result.min()));
    ASSERT_TRUE(std::isnan(result.max()));
}

TEST(affine_unop, sqrt_unbounded) {
    // Center overflows to inf (bounds [inf, inf]); the Chebyshev terms used to yield NaN.
    auto result = (AffineForm(1e200) * AffineForm(1e200)).sqrt();
    ASSERT_EQ(result.min(), -INFINITY);
    ASSERT_EQ(result.max(), INFINITY);
}

TEST(affine_unop, pow_non_power_of_two) {
    ASSERT_DOUBLE_EQ(AffineForm(2.0).pow(5).center(), 32.0);
    ASSERT_DOUBLE_EQ(AffineForm(2.0).pow(6).center(), 64.0);
    ASSERT_DOUBLE_EQ(AffineForm(2.0).pow(7).center(), 128.0);

    auto result = AffineForm(Winterval(1, 2)).pow(6);
    ASSERT_LE(result.min(), 1.0);
    ASSERT_GE(result.max(), 64.0);
}

TEST(affine_unop, sqrt_soundness_below_one) {
    // sqrt(min) + sqrt(max) < 1 used to floor to 0 and divide by it.
    auto result = AffineForm(Winterval(0.01, 0.04)).sqrt();
    ASSERT_LE(result.min(), 0.1);
    ASSERT_GE(result.max(), 0.2);
    ASSERT_TRUE(std::isfinite(result.min()) && std::isfinite(result.max()));
}

TEST(affine_unop, sigmoid_large_input) {
    auto result = AffineForm(Winterval(1000, 1001)).sigmoid();
    ASSERT_NEAR(result.min(), 1.0, 1e-9);
    ASSERT_NEAR(result.max(), 1.0, 1e-9);
}

TEST(affine_unop, exp_overflow) {
    // exp(710) overflows: the chord slope was inf and the center inf - inf = NaN.
    auto result = AffineForm(Winterval(700, 710)).exp();
    ASSERT_LE(result.min(), std::exp(700.0));
    ASSERT_EQ(result.max(), INFINITY);
}

TEST(affine_unop, exp_underflow) {
    // exp underflows to 0 at both ends: the chord slope was 0 and log(0) = -inf gave NaN.
    auto result = AffineForm(Winterval(-800, -790)).exp();
    ASSERT_LE(result.min(), 0.0);
    ASSERT_GE(result.max(), 0.0);
    ASSERT_TRUE(std::isfinite(result.min()) && std::isfinite(result.max()));
}

TEST(affine_unop, tanh_sigmoid_huge_bounds) {
    // Saturated endpoints: min_deviation = 0 times an overflowing width gave NaN.
    auto x = AffineForm(Winterval(1e308, 1.7e308));
    ASSERT_DOUBLE_EQ(x.tanh().min(), 1.0);
    ASSERT_DOUBLE_EQ(x.tanh().max(), 1.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().min(), 1.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().max(), 1.0);
}

TEST(affine_unop, tanh_sigmoid_unbounded) {
    auto x = AffineForm(1.0) / AffineForm(Winterval(-1, 1));
    ASSERT_DOUBLE_EQ(x.tanh().min(), -1.0);
    ASSERT_DOUBLE_EQ(x.tanh().max(), 1.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().min(), 0.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().max(), 1.0);
}

TEST(affine_unop, tanh_sigmoid_nan_interval) {
    // Center overflow plus an infinite term concretizes to NaN; tanh and sigmoid are still bounded.
    auto x = AffineForm(1e200) * AffineForm(1e200) + AffineForm(1.0) / AffineForm(Winterval(-1, 1));
    ASSERT_DOUBLE_EQ(x.tanh().min(), -1.0);
    ASSERT_DOUBLE_EQ(x.tanh().max(), 1.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().min(), 0.0);
    ASSERT_DOUBLE_EQ(x.sigmoid().max(), 1.0);
}
