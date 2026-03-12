//
// Created by will on 2/17/26.
//

#include "gtest/gtest.h"
#include "Caffeine/AffineForm.hpp"

TEST(affine_misc, min_max) {
    // Just a sanity test to make sure the test suite is working.
    auto a = AffineForm(Winterval(-2, 3));
    ASSERT_NEAR(a.min(), -2, 0.001);
    ASSERT_NEAR(a.max(), 3, 0.001);

    a = a - AffineForm(Winterval(4, 5));
    ASSERT_NEAR(a.min(), -7, 0.001);
    ASSERT_NEAR(a.max(), -1, 0.001);
}

TEST(affine_misc, add_noise_symbol) {
    auto a = AffineForm(Winterval(-2, 3));
    a.add_noise_symbol(5);
    ASSERT_NEAR(a.min(), -7, 0.001);
    ASSERT_NEAR(a.max(), 8, 0.001);
}

TEST(affine_misc, center_and_radius) {
    auto a = AffineForm(Winterval(1, 5));
    ASSERT_NEAR(a.center(), 3.0, 0.001);
    ASSERT_NEAR(a.radius(), 2.0, 0.001);
}

TEST(affine_misc, to_interval_roundtrip) {
    // An AffineForm constructed from an interval should concretize back to that same interval.
    ASSERT_EQ(AffineForm(Winterval(-3, 7)).to_interval(), Winterval(-3, 7));
}

TEST(affine_misc, default_constructor) {
    // Default-constructed form should be the zero point.
    AffineForm a;
    ASSERT_NEAR(a.center(), 0.0, 0.001);
    ASSERT_NEAR(a.radius(), 0.0, 0.001);
}
