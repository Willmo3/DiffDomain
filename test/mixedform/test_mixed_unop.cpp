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

TEST(mixed_unop, test_mixed_sigmoid) {
    auto form = MixedForm(Winterval(1, 2));
    auto result = form.sigmoid();
    ASSERT_NEAR(0.731058578630, result.min(), 1e-6);
    ASSERT_NEAR(0.880797077977, result.max(), 1e-6);
}

TEST(mixed_unop, test_mixed_relu) {
    auto form = MixedForm(Winterval(-2, 3));
    auto result = form.relu();
    EXPECT_EQ(result.min(), 0.0);
    EXPECT_EQ(result.max(), 3.0);
}

TEST(mixed_unop, test_mixed_sqrt) {
    auto form = MixedForm(Winterval(1, 4));
    auto result = form.sqrt();
    ASSERT_NEAR(1.0, result.min(), 1e-6);
    ASSERT_NEAR(2.0, result.max(), 1e-6);
}

TEST(mixed_unop, test_mixed_sqrt_point_interval) {
    auto form = MixedForm(Winterval(4.0));
    auto result = form.sqrt();
    ASSERT_NEAR(result.min(), 2.0, 1e-5);
    ASSERT_NEAR(result.max(), 2.0, 1e-5);
}

TEST(mixed_unop, test_mixed_sqrt_zero) {
    auto form = MixedForm(Winterval(0.0));
    auto result = form.sqrt();
    ASSERT_NEAR(result.min(), 0.0, 1e-5);
    ASSERT_NEAR(result.max(), 0.0, 1e-5);
}

TEST(mixed_unop, test_mixed_sqrt_1_to_4) {
    auto form = MixedForm(Winterval(1.0, 4.0));
    auto result = form.sqrt();
    
    // sqrt(1) = 1, sqrt(4) = 2
    ASSERT_LE(result.min(), 1.0);
    ASSERT_GE(result.max(), 2.0);
}

TEST(mixed_unop, test_mixed_sqrt_9_to_25) {
    auto form = MixedForm(Winterval(9.0, 25.0));
    auto result = form.sqrt();
    
    // sqrt(9) = 3, sqrt(25) = 5
    ASSERT_LE(result.min(), 3.0);
    ASSERT_GE(result.max(), 5.0);
}

TEST(mixed_unop, test_mixed_sqrt_negative_interval) {
    auto form = MixedForm(Winterval(-4.0, -1.0));
    auto result = form.sqrt();
    
    // Both min and max should be NaN for negative intervals
    ASSERT_TRUE(std::isnan(result.min()));
    ASSERT_TRUE(std::isnan(result.max()));
}

TEST(mixed_unop, test_mixed_sqrt_preserves_relationships) {
    // sqrt should preserve order relationships for positive intervals
    auto form1 = MixedForm(Winterval(1.0, 4.0));
    auto form2 = MixedForm(Winterval(4.0, 9.0));
    
    auto result1 = form1.sqrt();
    auto result2 = form2.sqrt();
    
    // result1 should still be less than result2
    ASSERT_LT(result1.max(), result2.min());
}