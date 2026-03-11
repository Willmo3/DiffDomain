//
// Created by will on 12/7/25.
//

#include "gtest/gtest.h"
#include "MixedForm/MixedForm.hpp"

TEST(mixed_comparison, binary_equality) {
    auto a = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(a == a);
}
TEST(mixed_comparison, binary_comparison) {
    auto a = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(3, 6)), Winterval(0, 3));
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(a <= b);
    EXPECT_FALSE(a > b);
    EXPECT_FALSE(a >= b);
}

TEST(mixed_comparison, scalar_gt) {
    auto form = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    EXPECT_TRUE(form > 0);
    EXPECT_FALSE(form > 1);
}
TEST(mixed_comparison, scalar_lt) {
    auto form = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    EXPECT_TRUE(form < 5);
    EXPECT_FALSE(form < 2);
}
TEST(mixed_comparison, scalar_lte) {
    auto form = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    EXPECT_TRUE(form <= 2);
}
TEST(mixed_comparison, scalar_gte) {
    auto form = MixedForm(AffineForm(Winterval(1, 4)), Winterval(-1, 2));
    EXPECT_TRUE(form >= 1);
}