//
// Created by will on 2/17/26.
//

#include "gtest/gtest.h"
#include "MixedForm/MixedForm.hpp"

TEST(mixed_misc, test_mixed_min) {
    auto form = MixedForm(Winterval(-3, 2));
    EXPECT_EQ(form.min(), -3);
    EXPECT_EQ(form.max(), 2);

    form = form * MixedForm(Winterval( -13, 2));
    EXPECT_EQ(form.min(), -26);
    EXPECT_EQ(form.max(), 39);
}

TEST(mixed_misc, construction_from_winterval_only) {
    // Constructing from a Winterval alone should use the interval as both representations.
    auto form = MixedForm(Winterval(2, 6));
    EXPECT_EQ(form.interval_bounds(), Winterval(2, 6));
}

TEST(mixed_misc, construction_from_affineform_only) {
    // Constructing from an AffineForm alone should concretise via to_interval().
    auto form = MixedForm(AffineForm(Winterval(2, 6)));
    EXPECT_EQ(form.interval_bounds(), Winterval(2, 6));
}

TEST(mixed_misc, radius) {
    EXPECT_EQ(MixedForm(Winterval(2, 6)).radius(), 2.0);
}

TEST(mixed_misc, interval_tightening_wide_affine) {
    // A wide affine form clipped by a tight Winterval gives the tighter bounds.
    auto form = MixedForm(AffineForm(Winterval(0, 8)), Winterval(2, 5));
    EXPECT_EQ(form.interval_bounds(), Winterval(2, 5));
}

TEST(mixed_misc, interval_tightening_wide_winterval) {
    // A tight affine form is not widened by a broad Winterval.
    auto form = MixedForm(AffineForm(Winterval(2, 5)), Winterval(0, 8));
    EXPECT_EQ(form.interval_bounds(), Winterval(2, 5));
}

TEST(mixed_misc, split_into_two) {
    auto splits = MixedForm(Winterval(0, 6)).split(2);
    ASSERT_EQ(splits.size(), 2u);
    EXPECT_EQ(splits[0].interval_bounds(), Winterval(0, 3));
    EXPECT_EQ(splits[1].interval_bounds(), Winterval(3, 6));
}
