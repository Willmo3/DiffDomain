//
// Created by will on 10/9/25.
//

#include <cmath>
#include <gtest/gtest.h>
#include "MixedForm/MixedForm.hpp"

TEST(mixed_scalar, addition) {
    auto form = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto result = form + 2.0;
    EXPECT_EQ(result.interval_bounds(), Winterval(4, 4));
}
TEST(mixed_scalar, subtraction) {
    auto form = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto result = form - 2.0;
    EXPECT_EQ(result.interval_bounds(), Winterval(0, 0));
}
TEST(mixed_scalar, multiplication) {
    auto form = MixedForm(AffineForm(Winterval(0, 3)), Winterval(-1, 2));
    auto result = form * 2.0;
    EXPECT_EQ(result.interval_bounds(), Winterval(0, 4));
}
TEST(mixed_scalar, division) {
    auto form = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    auto result = form / 2.0;
    EXPECT_EQ(result.interval_bounds(), Winterval(0.5, 1));

    form = MixedForm(AffineForm(Winterval(-4, 2)), Winterval(-2, 1));
    result = form / 0;
    EXPECT_EQ(result.interval_bounds(), Winterval(-INFINITY, INFINITY));
}
TEST(mixed_scalar, mult_by_zero) {
    // Multiplying by 0 must collapse bounds to [0, 0].
    auto form = MixedForm(Winterval(1, 5));
    EXPECT_EQ((form * 0.0).interval_bounds(), Winterval(0, 0));
}

TEST(mixed_scalar, mult_by_negative) {
    // [1,3] * -2 = [-6, -2]
    auto form = MixedForm(Winterval(1, 3));
    EXPECT_EQ((form * -2.0).interval_bounds(), Winterval(-6, -2));
}

TEST(mixed_scalar, add_negative_scalar) {
    // Adding a negative scalar shifts bounds left.
    auto form = MixedForm(Winterval(2, 5));
    EXPECT_EQ((form + -1.0).interval_bounds(), Winterval(1, 4));
}
