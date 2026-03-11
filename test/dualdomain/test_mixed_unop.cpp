//
// Created by will on 10/17/25.
//

#include <gtest/gtest.h>
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