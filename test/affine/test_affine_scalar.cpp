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