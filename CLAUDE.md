# DiffDomain

C++20 numeric abstract domains for abstract interpretation of PDE solvers (co-designed with PDEnclose).

## Build & test
- Deps: `cd scripts && ./update_deps` (must run from `scripts/`; fetches cereal, Eigen, Polylib into `lib/`, which is gitignored).
- Build: `mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make`
- Test binaries land in `out/tests/`: `test_caffeine`, `test_winterval`, `test_mixed`, `test_dualnumber`, `test_samplerange` (GoogleTest, fetched via FetchContent).
- `cd scripts && ./run_unit_tests` runs all of them.

## Layout
- `src/Numeric.hpp` — the `Numeric` concept every domain must satisfy (arithmetic, scalar ops, pow/sqrt/abs/exp/tanh/sigmoid/relu, comparisons, radius/min/max, union_with/split, `T(double)`, `operator<<`).
- `src/Real/` — double wrapper.
- `src/Winterval/` — interval arithmetic (base for others).
- `src/Caffeine/` — affine forms (`AffineForm`); depends on Winterval. Noise symbols capped by `MAX_NOISE_SYMBOLS`; timing macros `AFFINE_TIME_*`.
- `src/MixedForm/` — reduced product of AffineForm + Winterval.
- `src/SampleRange/` — unsound sampling domain (Eigen, `SAMPLED_VALUES` samples).
- `src/DualNumber/` — header-only template `DualNumber<T: Numeric>` for forward-mode AD. `synthesis/` holds synthesized derivative transformers (product, quotient, tanh) gated by `USE_SYNTHESIZED_ZONOS`; uses polylib + Eigen.
- `src/CMakeLists.txt` — per-domain libs plus a unified `domains` lib linking all of them and polylib.
- `test/<domain>/` — gtest files split into binop/unop/scalar/comparison/misc.
- `release/diffdomain-release.tar.gz` — packaged sources for external CMake consumers.
- `dockerfile` — Fedora 40 build env.

## Gotchas
- AffineForm has global noise-symbol state (`max_noise_symbol`), so ids depend on everything allocated earlier. Tests must refer to symbols relative to `AffineForm::next_noise_symbol()`, never by absolute id — test run order is not stable (LTO reverses it).
- Includes are rooted at `src/` and `lib/` (e.g. `#include "Winterval/Winterval.hpp"`, `"cereal/..."`).
- New domain ops must be added to every domain to keep satisfying `Numeric`.
