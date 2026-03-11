//
// Created by will on 2/17/26.
//

#include "gtest/gtest.h"
#include "DualDomain/MixedForm.hpp"

TEST(mixed_misc, test_mixed_min) {
    auto form = MixedForm(Winterval(-3, 2));
    EXPECT_EQ(form.min(), -3);
    EXPECT_EQ(form.max(), 2);

    form = form * MixedForm(Winterval( -13, 2));
    EXPECT_EQ(form.min(), -26);
    EXPECT_EQ(form.max(), 39);
}