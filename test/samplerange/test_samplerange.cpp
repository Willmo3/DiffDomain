//
// Test scalar operations for SampleRange
// Created by will on 3/26/26.
//
#include <gtest/gtest.h>
#include <cmath>
#include "SampleRange/SampleRange.hpp"

TEST(samplerange_scalar, multiply_by_zero) {
    SampleRange a(2.0, 5.0);
    SampleRange result = a * 0.0;

    EXPECT_NEAR(result.min(), 0.0, 1e-10);
    EXPECT_NEAR(result.max(), 0.0, 1e-10);
}

TEST(samplerange_comparison, samplerange_inequality) {
    SampleRange a(1.0, 5.0);
    SampleRange b(1.0, 5.0);

    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

TEST(samplerange_comparison, scalar_less_than) {
    SampleRange a(2.0, 5.0);

    // max < scalar
    EXPECT_TRUE(a < 6.0);
    EXPECT_FALSE(a < 4.0);
}

TEST(samplerange_comparison, scalar_greater_than) {
    SampleRange a(2.0, 5.0);

    // min > scalar
    EXPECT_TRUE(a > 1.0);
    EXPECT_FALSE(a > 3.0);
}

TEST(samplerange_misc, default_constructor) {
    SampleRange a;
}

TEST(samplerange_misc, scalar_constructor) {
    SampleRange a(5.0);
}

TEST(samplerange_misc, range_constructor) {
    SampleRange a(2.0, 8.0);
}

TEST(samplerange_misc, linspace_properties) {
    SampleRange a(0.0, 10.0);
}

TEST(samplerange_misc, radius_calculation) {
    SampleRange a(3.0, 11.0);
}

TEST(samplerange_misc, zero_width_range) {
    SampleRange a(5.0, 5.0);

    EXPECT_NEAR(a.radius(), 0.0, 1e-10);
    EXPECT_EQ(a.min(), 5.0);
    EXPECT_EQ(a.max(), 5.0);
}

TEST(samplerange_misc, negative_range) {
    SampleRange a(-10.0, -2.0);
}

TEST(samplerange_misc, large_split) {
    SampleRange a(0.0, 100.0);
    auto splits = a.split(10);

    ASSERT_EQ(splits.size(), 10);
}

TEST(samplerange_unary, tanh) {
    SampleRange a(0.0, 1.0);
    SampleRange result = a.tanh();

    EXPECT_GT(result.min(), std::tanh(0.0));
    EXPECT_LT(result.max(), std::tanh(1.0));
}
TEST(samplerange_unary, sigmoid) {
    SampleRange a(0.0, 1.0);
    SampleRange result = a.sigmoid();

    EXPECT_GT(result.min(), 0.5);
    EXPECT_LT(result.max(), 0.74);
}

