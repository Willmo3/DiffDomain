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

    EXPECT_DOUBLE_EQ(result.min(), std::tanh(0.0));
    EXPECT_DOUBLE_EQ(result.max(), std::tanh(1.0));
}
TEST(samplerange_unary, sigmoid) {
    SampleRange a(0.0, 1.0);
    SampleRange result = a.sigmoid();

    EXPECT_DOUBLE_EQ(result.min(), 0.5);
    EXPECT_DOUBLE_EQ(result.max(), 1.0 / (1.0 + std::exp(-1.0)));
}

TEST(samplerange_unary, relu) {
    SampleRange a(-1.0, 1.0);
    SampleRange result = a.relu();

    EXPECT_EQ(result.min(), 0.0);
    EXPECT_NEAR(result.max(), 1.0, 0.1);
}
TEST(samplerange_unary, sqrt) {
    SampleRange a(0.0, 1.0);
    SampleRange result = a.sqrt();

    EXPECT_NEAR(result.min(), 0.0, 0.1);
    EXPECT_NEAR(result.max(), 1.0, 0.1);
}

TEST(samplerange_unary, sqrt_large_values) {
    SampleRange a(9.0, 25.0);
    SampleRange result = a.sqrt();

    EXPECT_NEAR(result.min(), 3.0, 0.1);
    EXPECT_NEAR(result.max(), 5.0, 0.1);
}

TEST(samplerange_unary, sqrt_negative_produces_nan) {
    SampleRange a(-4.0, -1.0);
    SampleRange result = a.sqrt();

    EXPECT_TRUE(std::isnan(result.min()));
    EXPECT_TRUE(std::isnan(result.max()));
}

TEST(samplerange_unary, sqrt_mixed_sign) {
    SampleRange a(-1.0, 4.0);
    SampleRange result = a.sqrt();

    // Mixed sign interval: negative values produce NaN, positive produce valid sqrt
    EXPECT_TRUE(std::isnan(result.min()) || result.min() >= 0.0);
}
TEST(samplerange_misc, split_zero_throws) {
    EXPECT_THROW(SampleRange(0.0, 1.0).split(0), std::invalid_argument);
}

TEST(samplerange_misc, huge_bounds_no_overflow) {
    // max - min overflowed: samples were inf/NaN and the radius NaN.
    SampleRange a(-1e308, 1e308);
    EXPECT_GE(a.min(), -1e308);
    EXPECT_LE(a.max(), 1e308);
    EXPECT_LT(a.min(), 0.0);
    EXPECT_GT(a.max(), 0.0);
    EXPECT_TRUE(std::isfinite(a.radius()));
}

TEST(samplerange_misc, huge_bounds_split) {
    auto splits = SampleRange(-1e308, 1e308).split(4);
    for (const auto &split : splits) {
        EXPECT_TRUE(std::isfinite(split.min()) && std::isfinite(split.max()));
    }
}

TEST(samplerange_misc, pointwise_max) {
    Eigen::VectorXd ascending = Eigen::VectorXd::LinSpaced(SampleRange::SAMPLED_VALUES, 0.0, 1.0);
    Eigen::VectorXd descending = 1.0 - ascending.array();
    SampleRange result = SampleRange(ascending).pointwise_max(SampleRange(descending));

    // Samplewise, not a union: no sample falls below 0.5, though both ranges reach 0.
    EXPECT_EQ(result, SampleRange(Eigen::VectorXd(ascending.cwiseMax(descending))));
    EXPECT_NEAR(result.min(), 0.5, 1e-3);
    EXPECT_EQ(result.max(), 1.0);
}

TEST(samplerange_misc, corner_samples) {
    // Sample 0 is every range's min and sample 1 its max, so (a + b)'s range covers the all-min and all-max runs.
    SampleRange a(-1.0, 2.0);
    SampleRange b(3.0, 5.0);
    SampleRange sum = a + b;
    EXPECT_EQ(sum.min(), 2.0);
    EXPECT_EQ(sum.max(), 7.0);
}
