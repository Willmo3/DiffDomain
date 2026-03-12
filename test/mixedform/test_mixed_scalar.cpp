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
