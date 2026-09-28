//
// Created by will on 2/17/26.
//

#include "gtest/gtest.h"
#include <thread>
#include <unordered_set>

#include "Caffeine/AffineForm.hpp"

TEST(affine_misc, min_max) {
    // Just a sanity test to make sure the test suite is working.
    auto a = AffineForm(Winterval(-2, 3));
    ASSERT_NEAR(a.min(), -2, 0.001);
    ASSERT_NEAR(a.max(), 3, 0.001);

    a = a - AffineForm(Winterval(4, 5));
    ASSERT_NEAR(a.min(), -7, 0.001);
    ASSERT_NEAR(a.max(), -1, 0.001);
}

TEST(affine_misc, add_noise_symbol) {
    auto a = AffineForm(Winterval(-2, 3));
    a.add_noise_symbol(5);
    ASSERT_NEAR(a.min(), -7, 0.001);
    ASSERT_NEAR(a.max(), 8, 0.001);
}

TEST(affine_misc, center_and_radius) {
    auto a = AffineForm(Winterval(1, 5));
    ASSERT_NEAR(a.center(), 3.0, 0.001);
    ASSERT_NEAR(a.radius(), 2.0, 0.001);
}

TEST(affine_misc, to_interval_roundtrip) {
    // An AffineForm constructed from an interval should concretize back to that same interval.
    ASSERT_EQ(AffineForm(Winterval(-3, 7)).to_interval(), Winterval(-3, 7));
}

TEST(affine_misc, default_constructor) {
    // Default-constructed form should be the zero point.
    AffineForm a;
    ASSERT_NEAR(a.center(), 0.0, 0.001);
    ASSERT_NEAR(a.radius(), 0.0, 0.001);
}

TEST(affine_misc, noise_symbols_unique_across_threads) {
    // Built without OpenMP, like PDEnclose's copy: symbol allocation must still be thread-safe.
    constexpr int n_threads = 8, per_thread = 20000;
    std::vector<std::vector<AffineForm>> forms(n_threads);
    std::vector<std::thread> threads;
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&forms, t] {
            for (int i = 0; i < per_thread; ++i) forms[t].emplace_back(Winterval(0, 1));
        });
    }
    for (auto &th : threads) th.join();

    // Each form holds exactly one symbol; to_string prints it as "Noise symbols: (<symbol>: <coeff>),".
    std::unordered_set<std::string> symbols;
    for (auto &per : forms) {
        for (auto &form : per) {
            auto s = form.to_string();
            auto start = s.find("Noise symbols: (") + 16;
            symbols.insert(s.substr(start, s.find(':', start) - start));
        }
    }
    ASSERT_EQ(symbols.size(), static_cast<size_t>(n_threads * per_thread));
}
