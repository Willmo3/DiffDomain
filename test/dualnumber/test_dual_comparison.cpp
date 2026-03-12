//
// Created by will on 3/12/26.
//

#include <gtest/gtest.h>

#include "DualNumber/DualNumber.hpp"
#include "Winterval/Winterval.hpp"

// ── Equality / inequality ────────────────────────────────────────────────────

TEST(dual_comparison, equal_to_self) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(x == x);
    EXPECT_FALSE(x != x);
}

TEST(dual_comparison, equal_same_values) {
    // Two independently constructed objects with identical primal and deriv.
    DualNumber a(Winterval(3, 3), Winterval(1, 1));
    DualNumber b(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);
}

TEST(dual_comparison, not_equal_different_primal) {
    DualNumber a(Winterval(3, 3), Winterval(1, 1));
    DualNumber b(Winterval(5, 5), Winterval(1, 1));
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

TEST(dual_comparison, not_equal_different_deriv) {
    // Same primal, different derivative: objects are not equal.
    DualNumber a(Winterval(3, 3), Winterval(1, 1));
    DualNumber b(Winterval(3, 3), Winterval(2, 2));
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

// ── Ordering (based on primal only) ──────────────────────────────────────────

TEST(dual_comparison, ordering_strict) {
    // Derivative values differ to confirm ordering ignores them.
    DualNumber a(Winterval(3, 3), Winterval(1, 1));
    DualNumber b(Winterval(5, 5), Winterval(0, 0));
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(a <= b);
    EXPECT_FALSE(a > b);
    EXPECT_FALSE(a >= b);
    EXPECT_TRUE(b > a);
    EXPECT_TRUE(b >= a);
}

TEST(dual_comparison, ordering_equal_primals) {
    // Equal primals → a <= b and a >= b both true, strict inequalities false.
    DualNumber a(Winterval(3, 3), Winterval(1, 1));
    DualNumber b(Winterval(3, 3), Winterval(2, 2));
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a >= b);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(a > b);
}

// ── Scalar comparisons ────────────────────────────────────────────────────────

TEST(dual_comparison, scalar_gt) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(x > 2.0);
    EXPECT_FALSE(x > 3.0);
}

TEST(dual_comparison, scalar_lt) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(x < 4.0);
    EXPECT_FALSE(x < 3.0);
}

TEST(dual_comparison, scalar_gte) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(x >= 3.0);
    EXPECT_FALSE(x >= 4.0);
}

TEST(dual_comparison, scalar_lte) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    EXPECT_TRUE(x <= 3.0);
    EXPECT_FALSE(x <= 2.0);
}

