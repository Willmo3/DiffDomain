//
// Created by will on 9/30/25.
//

#include <math.h>
#include <gtest/gtest.h>

#include "Caffeine/AffineForm.hpp"

TEST(affine_scalar, mult) {
    auto base = AffineForm(Winterval(-2, 3)) * -2;
    auto min = base.to_interval().min();
    auto max = base.to_interval().max();

    ASSERT_EQ(-6, min);
    ASSERT_EQ(4, max);
}
TEST(affine_scalar, div) {
    auto base = AffineForm(Winterval(-2, 3)) / -2;
    auto min = base.to_interval().min();
    auto max = base.to_interval().max();

    ASSERT_EQ(-1.5, min);
    ASSERT_EQ(1, max);

    // TODO: fix scalar division divide by 0 to be consistent with inverse multiplication.
    base = base / 0;
    // ASSERT_EQ(-INFINITY, base.to_interval().min());
    // ASSERT_EQ(INFINITY, base.to_interval().max());
}
TEST(affine_scalar, add) {
    auto base = AffineForm(Winterval(-2, 3)) + 2;
    auto min = base.to_interval().min();
    auto max = base.to_interval().max();
    ASSERT_EQ(0, min);
    ASSERT_EQ(5, max);

    base = base + -1;
    ASSERT_EQ(-1, base.to_interval().min());
    ASSERT_EQ(4, base.to_interval().max());
}
TEST(affine_scalar, sub) {
    auto base = AffineForm(Winterval(-2, 3)) - 2;
    auto min = base.to_interval().min();
    auto max = base.to_interval().max();

    ASSERT_EQ(-4, min);
    ASSERT_EQ(1, max);
}

TEST(affine_scalar, mult_by_zero) {
    // Multiplying by zero must collapse the form to a point at zero.
    auto result = AffineForm(Winterval(-3, 7)) * 0.0;
    ASSERT_NEAR(result.center(), 0.0, 0.001);
    ASSERT_NEAR(result.radius(), 0.0, 0.001);
}

TEST(affine_scalar, mult_by_one_identity) {
    // Multiplying by one must preserve the interval.
    auto a = AffineForm(Winterval(-3, 7));
    ASSERT_EQ((a * 1.0).to_interval(), a.to_interval());
}

TEST(affine_scalar, mult_by_negative) {
    // [-2,3] * -1 = [-3, 2]
    auto result = AffineForm(Winterval(-2, 3)) * -1.0;
    ASSERT_NEAR(result.min(), -3.0, 0.001);
    ASSERT_NEAR(result.max(),  2.0, 0.001);
}
