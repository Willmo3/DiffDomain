//
// Created by will on 3/12/26.
//

#include <gtest/gtest.h>

#include "DualNumber/DualNumber.hpp"
#include "Winterval/Winterval.hpp"

// Helpers:
//   x  – the independent variable at x = 3, dx/dx = 1
//   y  – a constant                 y = 2, dy/dx = 0
static DualNumber<Winterval> makeX() { return { Winterval(3, 3), Winterval(1, 1) }; }
static DualNumber<Winterval> makeY() { return { Winterval(2, 2), Winterval(0, 0) }; }

// ── Addition ──────────────────────────────────────────────────────────────────

TEST(dual_binop, addition_primal) {
    // primal: 3 + 2 = 5
    EXPECT_EQ((makeX() + makeY()).min(), 5.0);
}

TEST(dual_binop, addition_sum_rule) {
    // (f + g)' = f' + g' = 1 + 0 = 1
    EXPECT_EQ(makeX() + makeY(), DualNumber(Winterval(5, 5), Winterval(1, 1)));
}

// ── Subtraction ───────────────────────────────────────────────────────────────

TEST(dual_binop, subtraction_sum_rule) {
    // (f - g)' = f' - g' = 1 - 0 = 1
    EXPECT_EQ(makeX() - makeY(), DualNumber(Winterval(1, 1), Winterval(1, 1)));
}

// ── Multiplication ────────────────────────────────────────────────────────────

TEST(dual_binop, multiplication_product_rule_constant) {
    // x * y where y is constant:  (f * g)' = f'g + fg' = 1*2 + 3*0 = 2
    EXPECT_EQ(makeX() * makeY(), DualNumber(Winterval(6, 6), Winterval(2, 2)));
}

TEST(dual_binop, multiplication_product_rule_both_variable) {
    // x * x:  (x * x)' = 1*3 + 3*1 = 6  (i.e. d/dx(x^2) = 2x at x=3)
    auto x = makeX();
    EXPECT_EQ(x * x, DualNumber(Winterval(9, 9), Winterval(6, 6)));
}

// ── Division ──────────────────────────────────────────────────────────────────

TEST(dual_binop, division_quotient_rule) {
    // x=4, y=2 (constant):  (f/g)' = (f'g - fg') / g^2 = (1*2 - 4*0) / 4 = 0.5
    DualNumber x(Winterval(4, 4), Winterval(1, 1));
    auto result = x / makeY();
    EXPECT_EQ(result, DualNumber(Winterval(2, 2), Winterval(0.5, 0.5)));
}

// ── Compositional ─────────────────────────────────────────────────────────────

TEST(dual_binop, union_with) {
    DualNumber a(Winterval(1, 3), Winterval(0, 1));
    DualNumber b(Winterval(2, 4), Winterval(1, 2));
    // Union of each component independently.
    EXPECT_EQ(a.union_with(b), DualNumber(Winterval(1, 4), Winterval(0, 2)));
}

