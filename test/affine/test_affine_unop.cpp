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

    value = AffineForm(Winterval(-1, 2)).abs();
    ASSERT_NEAR(value.center(), 0.25, 0.001);
    ASSERT_NEAR(value.radius(), 0.75, 0.001);
}

TEST(affine_unop, pow) {
    auto value = AffineForm(Winterval(1, 2)).pow(3);
    ASSERT_NEAR(value.center(), 3.375, 0.001);
    ASSERT_NEAR(value.radius(), 4.625, 0.001);
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

TEST(affine_unop, pow_two_positive_interval_soundness) {
    auto result = AffineForm(Winterval(2, 4)).pow(2);
    ASSERT_LE(result.min(), 4.0);
    ASSERT_GE(result.max(), 16.0);
}

TEST(affine_unop, exp) {
    auto result = AffineForm(Winterval(1, 2)).exp();
    ASSERT_NEAR(result.min(), -2.816363, 1e-5);
    ASSERT_NEAR(result.max(), 11.771873, 1e-5);
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