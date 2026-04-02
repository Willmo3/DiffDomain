//
// Created by will on 10/17/25.
//

#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include "MixedForm/MixedForm.hpp"

TEST(mixed_unop, test_mixed_abs) {
    auto form = MixedForm(Winterval(-10, 3));
    auto result = form.abs();
    EXPECT_EQ(result.interval_bounds(), Winterval(0, 10));
}
TEST(mixed_unop, test_mixed_pow) {
    auto form = MixedForm(Winterval(-3, 2));
    auto result = form.pow(2);
    EXPECT_EQ(result.interval_bounds(), Winterval(-6, 9));
}
TEST(mixed_unop, test_mixed_pow_zero) {
    // x^0 = 1 for any x ≠ 0.
    auto form = MixedForm(Winterval(2, 4));
    EXPECT_EQ(form.pow(0).interval_bounds(), Winterval(1, 1));
}

TEST(mixed_unop, test_mixed_pow_one_identity) {
    // x^1 = x.
    auto form = MixedForm(Winterval(2, 4));
    EXPECT_EQ(form.pow(1).interval_bounds(), form.interval_bounds());
}

TEST(mixed_unop, test_mixed_abs_positive_interval) {
    auto form = MixedForm(Winterval(1, 5));
    EXPECT_EQ(form.abs().interval_bounds(), Winterval(1, 5));
}

TEST(mixed_unop, test_mixed_exp) {
    auto form = MixedForm(Winterval(1, 2));
    auto result = form.exp();
    ASSERT_NEAR(2.718282, result.min(), 1e-6);
    ASSERT_NEAR(7.389056, result.max(), 1e-6);
}

TEST(mixed_unop, test_mixed_tanh) {
    auto form = MixedForm(Winterval(1, 2));
    auto result = form.tanh();
    ASSERT_NEAR(0.761594155956, result.min(), 1e-6);
    ASSERT_NEAR(0.964027580076, result.max(), 1e-6);
}