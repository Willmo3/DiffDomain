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