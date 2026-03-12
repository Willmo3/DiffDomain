//
// Created by will on 3/12/26.
//

#include <gtest/gtest.h>
#include <sstream>

#include "DualNumber/DualNumber.hpp"
#include "Numeric.hpp"
#include "Winterval/Winterval.hpp"

// Verify at compile time that DualNumber<Winterval> satisfies the Numeric concept.
static_assert(Numeric<DualNumber<Winterval>>, "DualNumber<Winterval> must satisfy Numeric");

// ── Accessors delegate to the primal ─────────────────────────────────────────

TEST(dual_misc, min_max_radius) {
    // Primal: [2, 6], radius = (6 - 2) / 2 = 2.
    DualNumber x(Winterval(2, 6), Winterval(0, 1));
    EXPECT_EQ(x.min(), 2.0);
    EXPECT_EQ(x.max(), 6.0);
    EXPECT_EQ(x.radius(), 2.0);
}

// ── split ─────────────────────────────────────────────────────────────────────

TEST(dual_misc, split_produces_correct_count) {
    DualNumber x(Winterval(0, 4), Winterval(0, 2));
    EXPECT_EQ(x.split(4).size(), 4u);
}

TEST(dual_misc, split_components_in_lockstep) {
    // Primal [0,4] and deriv [0,2] each split into 2 equal pieces.
    DualNumber x(Winterval(0, 4), Winterval(0, 2));
    auto splits = x.split(2);
    ASSERT_EQ(splits.size(), 2u);
    EXPECT_EQ(splits[0], DualNumber(Winterval(0, 2), Winterval(0, 1)));
    EXPECT_EQ(splits[1], DualNumber(Winterval(2, 4), Winterval(1, 2)));
}

// ── union_with ────────────────────────────────────────────────────────────────

TEST(dual_misc, union_with_merges_both_components) {
    DualNumber a(Winterval(1, 3), Winterval(0, 1));
    DualNumber b(Winterval(2, 4), Winterval(1, 2));
    EXPECT_EQ(a.union_with(b), DualNumber(Winterval(1, 4), Winterval(0, 2)));
}

// ── stream output ─────────────────────────────────────────────────────────────

TEST(dual_misc, stream_output_is_non_empty) {
    DualNumber x(Winterval(3, 3), Winterval(1, 1));
    std::ostringstream oss;
    oss << x;
    EXPECT_FALSE(oss.str().empty());
}

