//
// Created by will on 12/7/25.
//

#include "gtest/gtest.h"
#include "Caffeine/AffineForm.hpp"

TEST(affine_comparison, scalar_comparisons) {
    auto base = AffineForm(Winterval(-2, 3));
    ASSERT_FALSE(base > -2);
    ASSERT_TRUE(base >= -2);
    ASSERT_FALSE(base < 3);
    ASSERT_TRUE(base <= 3);
}

TEST(affine_comparison, binary_comparisons) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = AffineForm(Winterval(4, 5));
    ASSERT_FALSE(base > next);
    ASSERT_FALSE(base >= next);
    ASSERT_TRUE(base < next);
    ASSERT_TRUE(base <= next);
}

TEST(affine_comparison, equality) {
    auto first = AffineForm(Winterval(-2, 3));
    ASSERT_TRUE(first == first);

    auto third = first + AffineForm(Winterval(1, 1));
    ASSERT_FALSE(first == third);
}

TEST(affine_comparison, inequality) {
    auto a = AffineForm(Winterval(-2, 3));
    auto b = AffineForm(Winterval(4, 5));
    ASSERT_TRUE(a != b);
    ASSERT_FALSE(a != a);
}

TEST(affine_comparison, ordering_with_overlap) {
    // Overlapping intervals: neither is strictly greater than the other.
    auto a = AffineForm(Winterval(0, 4));
    auto b = AffineForm(Winterval(2, 6));
    ASSERT_FALSE(a < b);
    ASSERT_FALSE(a > b);
    ASSERT_FALSE(a <= b);
    ASSERT_FALSE(a >= b);
}

TEST(affine_comparison, scalar_gte_at_boundary) {
    // Value is exactly at the lower bound – >= should be true, > should be false.
    auto a = AffineForm(Winterval(2, 5));
    ASSERT_TRUE(a >= 2);
    ASSERT_FALSE(a > 2);
}

TEST(affine_comparison, scalar_lte_at_boundary) {
    // Value is exactly at the upper bound – <= should be true, < should be false.
    auto a = AffineForm(Winterval(2, 5));
    ASSERT_TRUE(a <= 5);
    ASSERT_FALSE(a < 5);
}
