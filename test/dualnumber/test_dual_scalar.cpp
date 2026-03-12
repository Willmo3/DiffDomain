//
// Created by will on 3/12/26.
//

#include <gtest/gtest.h>

#include "DualNumber/DualNumber.hpp"
#include "Winterval/Winterval.hpp"

// x = 3, dx = 2  (arbitrary non-unit derivative to make scaling visible)
static DualNumber<Winterval> makeX() { return { Winterval(3, 3), Winterval(2, 2) }; }

TEST(dual_scalar, addition) {
    // d/dx(f + c) = f'   primal shifts, derivative unchanged.
    EXPECT_EQ(makeX() + 5.0, DualNumber(Winterval(8, 8), Winterval(2, 2)));
}

TEST(dual_scalar, subtraction) {
    // d/dx(f - c) = f'
    EXPECT_EQ(makeX() - 1.0, DualNumber(Winterval(2, 2), Winterval(2, 2)));
}

TEST(dual_scalar, multiplication) {
    // d/dx(f * c) = f' * c = 2 * 4 = 8
    EXPECT_EQ(makeX() * 4.0, DualNumber(Winterval(12, 12), Winterval(8, 8)));
}

TEST(dual_scalar, division) {
    // d/dx(f / c) = f' / c = 2 / 2 = 1
    EXPECT_EQ(makeX() / 2.0, DualNumber(Winterval(1.5, 1.5), Winterval(1, 1)));
}

