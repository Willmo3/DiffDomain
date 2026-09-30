//
// Created by will on 9/30/25.
//
#include <gtest/gtest.h>
#include <random>

#include "Caffeine/AffineForm.hpp"

TEST(affine_binop, add) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base + AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 5, 0.001);
    ASSERT_NEAR(next.radius(), 3, 0.001);
}

TEST(affine_binop, sub) {
    // Symbols are allocated in order: s = base's, s + 1 = the rhs's.
    auto s = AffineForm::next_noise_symbol();
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base - AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), -4, 0.001);
    ASSERT_NEAR(next.radius(), 3, 0.001);

    // test correspondence of variables.
    next = next - base;
    ASSERT_NEAR(next.coeff_of(s), 0, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 1), 0.5, 0.001);
}

TEST(affine_binop, mult) {
    // s = base's, s + 1 = the rhs's, s + 2 = the multiplication error term.
    auto s = AffineForm::next_noise_symbol();
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base * AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 2.25, 0.001);
    ASSERT_NEAR(next.radius(), 12.75, 0.001);
    ASSERT_NEAR(next.coeff_of(s), -11.25, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 1), -0.25, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 2), 1.25, 0.001);
}

TEST(affine_binop, div) {
    // s = base's, s + 1 = the rhs's, s + 2 = the inverse's error term, s + 3 = the multiplication error term.
    auto s = AffineForm::next_noise_symbol();
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base / AffineForm(Winterval(4, 5));
    ASSERT_NEAR(next.center(), 0.112500, 0.001);
    ASSERT_NEAR(next.radius(), 0.637500, 0.001);
    ASSERT_NEAR(next.coeff_of(s), -0.562500, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 1), 0.010000, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 2), 0.002500, 0.001);
    ASSERT_NEAR(next.coeff_of(s + 3), 0.062500, 0.001);
}

TEST(affine_binop, union_with) {
    auto base = AffineForm(Winterval(-2, 3));
    auto next = base.union_with(AffineForm(Winterval(4, 5)));
    ASSERT_NEAR(next.center(), 1.5, 0.001);
    ASSERT_NEAR(next.radius(), 3.5, 0.001);
}

TEST(affine_binop, union_with_sound_under_shared_symbols) {
    // The join must contain each operand for every assignment of the shared symbols, not just its interval.
    // Operands are linear combinations of three shared symbols (scalar ops allocate no new symbols).
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(-2, 2);
    for (int trial = 0; trial < 500; ++trial) {
        auto s = AffineForm::next_noise_symbol();
        AffineForm x(Winterval(-1, 1)), y(Winterval(-1, 1)), z(Winterval(-1, 1));   // symbols s, s + 1, s + 2
        auto a = x * u(rng) + y * u(rng) + u(rng);
        auto b = (trial % 2 ? -a : x * u(rng) + z * u(rng)) + u(rng);           // include the x vs -x case
        auto joined = a.union_with(b);

        auto coeff = [&](const AffineForm &f, AffineForm::noise_symbol_t sym) { double c = f.coeff_of(sym); return c == c ? c : 0.0; };
        for (int e0 = -1; e0 <= 1; ++e0) for (int e1 = -1; e1 <= 1; ++e1) for (int e2 = -1; e2 <= 1; ++e2) {
            double eps[3] = {double(e0), double(e1), double(e2)};
            auto at = [&](const AffineForm &f) { double v = f.center(); for (int i = 0; i < 3; ++i) v += coeff(f, s + i) * eps[i]; return v; };
            // Symbols other than the shared three are free in [-1, 1].
            double free = joined.radius();
            for (int i = 0; i < 3; ++i) free -= std::abs(coeff(joined, s + i));
            ASSERT_LE(std::abs(at(a) - at(joined)), free + 1e-9) << "trial " << trial;
            ASSERT_LE(std::abs(at(b) - at(joined)), free + 1e-9) << "trial " << trial;
        }
    }
}

TEST(affine_binop, split) {
    auto base = AffineForm(Winterval(0, 6));
    auto splits = base.split(3);

    ASSERT_EQ(Winterval(0, 2), splits[0].to_interval());
    ASSERT_EQ(Winterval(2, 4), splits[1].to_interval());
    ASSERT_EQ(Winterval(4, 6), splits[2].to_interval());
}