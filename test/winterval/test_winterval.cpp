//
// Created by will on 9/30/25.
//

#include <cmath>
#include <gtest/gtest.h>

// Note: cereal root must be in the build path
#include "cereal/archives/json.hpp"
#include "Winterval/Winterval.hpp"

// TODO: split tests

TEST(winterval, contains) {
    auto winter = Winterval(0, 2);
    ASSERT_TRUE(winter.contains(0));
    ASSERT_FALSE(winter.contains(-1));
    ASSERT_TRUE(winter.contains(1));
}

TEST(winterval, add) {
    auto left = Winterval(0, 2);
    auto right = Winterval(-2, 1);
    auto sum = left + right;

    ASSERT_EQ(-2, sum.min());
    ASSERT_EQ(3, sum.max());
}

TEST(winterval, sub) {
    auto left = Winterval(0, 2);
    auto right = Winterval(-1, 2);
    auto difference = left - right;

    ASSERT_EQ(-2, difference.min());
    ASSERT_EQ(3, difference.max());
}

TEST(winterval, mul) {
    auto left = Winterval(1, 2);
    auto right = Winterval(-2, 1);
    auto mul = left * right;

    ASSERT_EQ(-4, mul.min());
    ASSERT_EQ(2, mul.max());
}

TEST(winterval, div) {
    auto left = Winterval(-1, 2);
    auto right = Winterval(-2, 1);
    auto div = left / right;

    ASSERT_EQ(-INFINITY, div.min());
    ASSERT_EQ(INFINITY, div.max());

    left = Winterval(0, 2);
    div = left / right;

    ASSERT_EQ(-1, div.min());
    ASSERT_EQ(INFINITY, div.max());

    left = Winterval(-2, 0);
    div = left / right;

    ASSERT_EQ(-INFINITY, div.min());
    ASSERT_EQ(1, div.max());

    right = Winterval(0, 0);
    div = left / right;

    ASSERT_EQ(0, div.min());
    ASSERT_EQ(0, div.max());
}

TEST(winterval, tanh) {
    auto left = Winterval(-1, 4).tanh().tanh();
    ASSERT_NEAR(-0.642015, left.min(), 1e-6);
    ASSERT_NEAR(0.761312, left.max(), 1e-6);
}

TEST(winterval, pow) {
    auto base = Winterval(-1, 4);

    auto pow_0 = base.pow(0);
    ASSERT_EQ(1, pow_0.min());
    ASSERT_EQ(1, pow_0.max());

    auto pow_1 = base.pow(1);
    ASSERT_EQ(-1, pow_1.min());
    ASSERT_EQ(4, pow_1.max());

    auto pow_2 = base.pow(2);
    ASSERT_EQ(-4, pow_2.min());
    ASSERT_EQ(16, pow_2.max());

    auto pow_3 = base.pow(3);
    ASSERT_EQ(-1, pow_3.min());
    ASSERT_EQ(64, pow_3.max());
}

TEST(winterval, scalar) {
    auto base = Winterval(-1, 4);

    auto prod = base * 2;
    ASSERT_EQ(-2, prod.min());
    ASSERT_EQ(8, prod.max());

    auto sum = base + 1;
    ASSERT_EQ(0, sum.min());
    ASSERT_EQ(5, sum.max());

    auto difference = base - 3;
    ASSERT_EQ(-4, difference.min());
    ASSERT_EQ(1, difference.max());

    auto quotient = base / 2;
    ASSERT_EQ(-0.5, quotient.min());
    ASSERT_EQ(2, quotient.max());

    auto div = base / 0;
    ASSERT_EQ(-INFINITY, div.min());
    ASSERT_EQ(INFINITY, div.max());
}

TEST(winterval, abs) {
    auto base = Winterval(-6, 4);
    auto abs = base.abs();

    ASSERT_EQ(0, abs.min());
    ASSERT_EQ(6, abs.max());
}

TEST(winterval, mid) {
    auto base = Winterval(-2, 2);
    ASSERT_EQ(0, base.mid());
}

TEST(winterval, radius) {
    ASSERT_EQ(1, Winterval(0, 2).radius());
    ASSERT_EQ(2, Winterval(-2, 2).radius());
    ASSERT_EQ(1.5, Winterval(0, 3).radius());
    ASSERT_EQ(0, Winterval(0, 0).radius());
}

TEST(winterval, serialize) {
    std::stringstream ss;
    // Scoping necessary to ensure data leaves scope.
    {
        cereal::JSONOutputArchive archive(ss);
        auto w = Winterval(-1, 2);
        archive( w );
    }

    cereal::JSONInputArchive archive(ss);
    Winterval w2 = Winterval(0, 0);
    archive( w2 );
    ASSERT_EQ(-1, w2.min());
    ASSERT_EQ(2, w2.max());
}

TEST(winterval, union_with) {
    auto a = Winterval(-2, 1);
    auto b = Winterval(0, 3);
    auto u = a.union_with(b);

    ASSERT_EQ(-2, u.min());
    ASSERT_EQ(3, u.max());
}

TEST(winterval, split) {
    auto a = Winterval(0, 8);
    auto splits = a.split(4);

    ASSERT_EQ(Winterval(0, 2), splits[0]);
    ASSERT_EQ(Winterval(2, 4), splits[1]);
    ASSERT_EQ(Winterval(4, 6), splits[2]);
    ASSERT_EQ(Winterval(6, 8), splits[3]);
}

TEST(winterval, comparisons) {
    auto a = Winterval(0, 2);
    auto b = Winterval(3, 5);
    auto c = Winterval(1, 4);

    ASSERT_TRUE(a < b);
    ASSERT_TRUE(a <= b);
    ASSERT_FALSE(a > b);
    ASSERT_FALSE(a >= b);

    ASSERT_FALSE(c < b);
    ASSERT_FALSE(c <= b);
    ASSERT_FALSE(c > b);
    ASSERT_FALSE(c >= b);

    ASSERT_FALSE(a <= 1);
    ASSERT_TRUE(a < 3);
    ASSERT_FALSE(a >= 3);
    ASSERT_FALSE(a > 2);
}

// ── Additional coverage ───────────────────────────────────────────────────────

TEST(winterval, reversed_bounds_constructor) {
    // Constructor must swap min/max when given in descending order.
    auto w = Winterval(5, 2);
    ASSERT_EQ(2, w.min());
    ASSERT_EQ(5, w.max());
}

TEST(winterval, equality_and_inequality) {
    ASSERT_TRUE(Winterval(1, 3) == Winterval(1, 3));
    ASSERT_FALSE(Winterval(1, 3) == Winterval(1, 4));
    ASSERT_TRUE(Winterval(1, 3) != Winterval(1, 4));
    ASSERT_FALSE(Winterval(1, 3) != Winterval(1, 3));
}

TEST(winterval, degenerate_point_interval) {
    auto w = Winterval(3, 3);
    ASSERT_EQ(3, w.min());
    ASSERT_EQ(3, w.max());
    ASSERT_EQ(3, w.mid());
    ASSERT_EQ(0, w.radius());
    ASSERT_TRUE(w.contains(3));
    ASSERT_FALSE(w.contains(2));
}

TEST(winterval, mul_both_mixed_sign) {
    // [-1,2] * [-2,3]: corner products are 2,-3,-4,6 → [-4, 6]
    auto result = Winterval(-1, 2) * Winterval(-2, 3);
    ASSERT_EQ(-4, result.min());
    ASSERT_EQ(6, result.max());
}

TEST(winterval, div_strictly_positive_divisor) {
    // [2,4] / [1,2]: candidates 2,4,1,2 → [1, 4]
    auto result = Winterval(2, 4) / Winterval(1, 2);
    ASSERT_EQ(1, result.min());
    ASSERT_EQ(4, result.max());
}

TEST(winterval, pow_even_positive_interval) {
    // [1,3]^2 = [1, 9]
    auto result = Winterval(1, 3).pow(2);
    ASSERT_EQ(1, result.min());
    ASSERT_EQ(9, result.max());
}

TEST(winterval, pow_odd_negative_interval) {
    // [-3,-1]^3 = [-27, -1] (monotone for odd powers)
    auto result = Winterval(-3, -1).pow(3);
    ASSERT_EQ(-27, result.min());
    ASSERT_EQ(-1, result.max());
}

TEST(winterval, abs_all_positive) {
    // abs is the identity for positive intervals.
    auto result = Winterval(2, 5).abs();
    ASSERT_EQ(2, result.min());
    ASSERT_EQ(5, result.max());
}

TEST(winterval, abs_all_negative) {
    // [-5,-2] → [2, 5]
    auto result = Winterval(-5, -2).abs();
    ASSERT_EQ(2, result.min());
    ASSERT_EQ(5, result.max());
}

TEST(winterval, scalar_negative_multiply) {
    // [1,3] * -2 = [-6, -2]
    auto result = Winterval(1, 3) * -2.0;
    ASSERT_EQ(-6, result.min());
    ASSERT_EQ(-2, result.max());
}

TEST(winterval, split_into_one) {
    // Splitting into one part must return the original interval.
    auto splits = Winterval(2, 6).split(1);
    ASSERT_EQ(1u, splits.size());
    ASSERT_EQ(Winterval(2, 6), splits[0]);
}

TEST(winterval, union_with_subsumes) {
    // When one interval contains the other, the union equals the larger one.
    ASSERT_EQ(Winterval(1, 5), Winterval(1, 5).union_with(Winterval(2, 3)));
    // Commutative.
    ASSERT_EQ(Winterval(1, 5), Winterval(2, 3).union_with(Winterval(1, 5)));
}

TEST(winterval, contains_exact_boundaries) {
    auto w = Winterval(1, 5);
    ASSERT_TRUE(w.contains(1));   // inclusive lower bound
    ASSERT_TRUE(w.contains(5));   // inclusive upper bound
    ASSERT_FALSE(w.contains(0));
    ASSERT_FALSE(w.contains(6));
}
