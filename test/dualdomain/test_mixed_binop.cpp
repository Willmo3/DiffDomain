//
// Created by will on 10/17/25.
//

#include <cmath>
#include <gtest/gtest.h>

#include "DualDomain/MixedForm.hpp"

TEST(mixed_binop, test_mixed_addition) {
    auto a = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(1, 4)), Winterval(0, 3));

    auto result = a + b;
    EXPECT_EQ(result.interval_bounds(), Winterval(3, 5));
}
TEST(mixed_binop, test_mixed_subtraction) {
    auto a = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(1, 4)), Winterval(0, 3));

    auto result = a - b;
    EXPECT_EQ(result.interval_bounds(), Winterval(-1, 1));
}
TEST(mixed_binop, test_mixed_multiplication) {
    auto a = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(1, 4)), Winterval(0, 3));

    auto result = a * b;
    EXPECT_EQ(result.interval_bounds(), Winterval(2, 6));
}
TEST(mixed_binop, test_mixed_division) {
    auto a = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(-1, 4)), Winterval(0, 2));
    auto result = a / b;
    // NOTE: weird behavior when dividing by interval that includes 0, because the affine reps produce a tighter, but still sound, bound.
    EXPECT_EQ(result.interval_bounds(), Winterval(1, 12));
}
/*
 * Compositional
 */
TEST(mixed_binop, test_mixed_union) {
    auto a = MixedForm(AffineForm(Winterval(2, 3)), Winterval(-1, 2));
    auto b = MixedForm(AffineForm(Winterval(1, 4)), Winterval(0, 3));
    auto result = a.union_with(b);
    EXPECT_EQ(result.interval_bounds(), Winterval(1, 3));
}
TEST(mixed_binop, test_mixed_split) {
    auto form = MixedForm(AffineForm(Winterval(0, 8)), Winterval(-1, 7));
    auto splits = form.split(4);

    EXPECT_EQ(splits[0].interval_bounds(), Winterval(0, 1.75));
    EXPECT_EQ(splits[1].interval_bounds(), Winterval(1.75, 3.5));
    EXPECT_EQ(splits[2].interval_bounds(), Winterval(3.5, 5.25));
    EXPECT_EQ(splits[3].interval_bounds(), Winterval(5.25, 7.0));
}