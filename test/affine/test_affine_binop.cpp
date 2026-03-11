//
// Created by will on 9/30/25.
//
#include <gtest/gtest.h>

#include "Caffeine/AffineForm.hpp"

TEST(affine_binop, add) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base + AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 5, 0.001);
    ASSERT_NEAR(next.radius(), 3, 0.001);
}

TEST(affine_binop, sub) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base - AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), -4, 0.001);
    ASSERT_NEAR(next.radius(), 3, 0.001);

    // test correspondence of variables.
    next = next - base;
    ASSERT_NEAR(next.coeff_of(2), 0, 0.001);
    ASSERT_NEAR(next.coeff_of(3), 0.5, 0.001);
}

TEST(affine_binop, mult) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base * AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 2.25, 0.001);
    ASSERT_NEAR(next.radius(), 12.75, 0.001);
    ASSERT_NEAR(next.coeff_of(4), -11.25, 0.001);
    ASSERT_NEAR(next.coeff_of(5), -0.25, 0.001);
    ASSERT_NEAR(next.coeff_of(6), 1.25, 0.001);
}

TEST(affine_binop, div) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base / AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 0.112500, 0.001);
    ASSERT_NEAR(next.radius(), 0.637500, 0.001);
    ASSERT_NEAR(next.coeff_of(7), -0.562500, 0.001);
    ASSERT_NEAR(next.coeff_of(8), 0.010000, 0.001);
    ASSERT_NEAR(next.coeff_of(9), 0.002500, 0.001);
    ASSERT_NEAR(next.coeff_of(10), 0.062500, 0.001);
}

TEST(affine_binop, union_with) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base.union_with(AffineForm(Winterval(4, 5)));
    ASSERT_NEAR(next.center(), 1.5, 0.001);
    ASSERT_NEAR(next.radius(), 3.5, 0.001);
}

TEST(affine_binop, split) {
    auto base = AffineForm(Winterval(0, 6));
    auto splits = base.split(3);

    ASSERT_EQ(Winterval(0, 2), splits[0].to_interval());
    ASSERT_EQ(Winterval(2, 4), splits[1].to_interval());
    ASSERT_EQ(Winterval(4, 6), splits[2].to_interval());
}