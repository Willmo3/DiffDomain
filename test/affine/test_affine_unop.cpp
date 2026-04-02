//
// Created by will on 9/30/25.
//

// NOTE: gtest keeps static state in a given test suite (i.e. affine).
// Therefore, we switch noise symbols.

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
    // Any form to the power 0 is the constant 1.
    auto value = AffineForm(Winterval(-5, 3)).pow(0);
    ASSERT_NEAR(value.center(), 1.0, 0.001);
    ASSERT_NEAR(value.radius(), 0.0, 0.001);
}

TEST(affine_unop, pow_one_identity) {
    // Power of 1 must preserve the interval bounds.
    auto a = AffineForm(Winterval(1, 4));
    ASSERT_EQ(a.pow(1).to_interval(), a.to_interval());
}

TEST(affine_unop, pow_two_soundness) {
    // [1,3]^2: exact range is [1,9]; the affine approximation must be a sound over-approximation.
    auto result = AffineForm(Winterval(1, 3)).pow(2);
    ASSERT_LE(result.min(), 1.0);
    ASSERT_GE(result.max(), 9.0);
}

TEST(affine_unop, pow_two_positive_interval_soundness) {
    // [2,4]^2: exact range [4,16]; result must contain it.
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
    ASSERT_NEAR(result.min(), -0.046326, 1e-5);
    ASSERT_NEAR(result.max(), 0.353398, 1e-5);
}
