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
