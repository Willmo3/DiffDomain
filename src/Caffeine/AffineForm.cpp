//
// Created by will on 9/23/25.
//

#include "AffineForm.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>

/*
 * Noise symbol management
 */
std::atomic<AffineForm::noise_symbol_t> AffineForm::max_noise_symbol = 0;

AffineForm::noise_symbol_t AffineForm::new_noise_symbol() {
    // Relaxed: only uniqueness matters, not ordering relative to other memory.
    return max_noise_symbol.fetch_add(1, std::memory_order_relaxed);
}
AffineForm::noise_symbol_t AffineForm::next_noise_symbol() {
    return max_noise_symbol.load(std::memory_order_relaxed);
}

void print_debug_info(const std::string &op_name, const std::chrono::high_resolution_clock::time_point &start_time) {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    std::cout << "[DEBUG] Operation " << op_name << " took " << duration << " microseconds." << std::endl;
}

/*
 * Constructors
 */
AffineForm::AffineForm(double center, std::unordered_map<noise_symbol_t, double> starting_coeffs):
    _center(center), _coefficients(std::move(starting_coeffs)) {

    // Reset noise symbols if needed.
    // NOTE: this branch will almost never be taken!
    if (_coefficients.size() > MAX_NOISE_SYMBOLS) {
        collapse();
    }
}
AffineForm::AffineForm(const Winterval &interval): _center((interval.min() + interval.max()) / 2),
    _coefficients(std::unordered_map<noise_symbol_t, double>()) {
    _coefficients.insert(std::pair(new_noise_symbol(), (interval.min() - interval.max()) / 2));
}
AffineForm::AffineForm(double value): _center(value), _coefficients(std::unordered_map<noise_symbol_t, double>()) {}
AffineForm::AffineForm(): _center(0), _coefficients(std::unordered_map<noise_symbol_t, double>()) {}

/*
 * Unary operators
 */
AffineForm AffineForm::operator-() const {
    auto value = clone();
    value._center = -_center;
    for (auto &coeff : value._coefficients | std::views::values) {
        coeff *= -1;
    }

    return value;
}
AffineForm AffineForm::abs() const {
#   ifdef AFFINE_TIME_ABS
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    // Strictly negative
    if (this->operator<(0)) {
        return this->operator*(-1);
    }
    // Strictly positive
    if (this->operator>(0)) {
        return clone();
    }
    // Straddles
    // Return abs(this.center / 2) + sum (this.noise / 2)
    auto value = clone();
    value._center = std::abs(value._center / 2);
    for (auto &coeff : value._coefficients | std::views::values) {
        coeff /= 2;
    }

#   ifdef AFFINE_TIME_ABS
    print_debug_info("absolute value", time);
#   endif

    return value;
}

/*
 * Affine arithmetic operators.
 */
AffineForm AffineForm::operator+(const AffineForm &other) const {
#   ifdef AFFINE_TIME_ADD
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    auto value = clone();
    value._center += other._center;

    for (auto [symbol, coeff] : other._coefficients) {
        // Outer product is union of both fields' error symbols. Common error symbols are added.
        if (auto [it, inserted] = value._coefficients.try_emplace(symbol, coeff); !inserted) {
            it->second += coeff;
        }
    }

#   ifdef AFFINE_TIME_ADD
    print_debug_info("addition", time);
#   endif

    if (value._coefficients.size() > MAX_NOISE_SYMBOLS) {
        value.collapse();
    }
    // Since affine addition introduces no new error, we don't need to add a new value!
    return value;
}
AffineForm AffineForm::operator-(const AffineForm &other) const {
#   ifdef AFFINE_TIME_SUB
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    auto value = clone();
    value._center -= other._center;

    for (auto [symbol, coeff] : other._coefficients) {
        if (auto [it, inserted] = value._coefficients.try_emplace(symbol, -coeff); !inserted) {
            it->second -= coeff;
        }
    }

#   ifdef AFFINE_TIME_SUB
    print_debug_info("subtraction", time);
#   endif

    if (value._coefficients.size() > MAX_NOISE_SYMBOLS) {
        value.collapse();
    }
    return value;
}
AffineForm AffineForm::operator*(const AffineForm &right) const {
#   ifdef AFFINE_TIME_MULT
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    auto result = AffineForm(this->_center * right._center, std::unordered_map<noise_symbol_t, double>());
    // Affine form multiplication is an outer product.

    // Perform product for all error symbols in rhs.
    for (auto [symbol, coeff] : right._coefficients) {
        auto a = this->_center * coeff;
        if (auto it = this->_coefficients.find(symbol); it == this->_coefficients.end()) {
            // Add missing error terms scaled by left's center.
            result._coefficients.emplace(symbol, a);
        } else {
            result._coefficients.emplace(symbol, a + right._center * it->second);
        }
    }

    // Now go and perform similar calculation for error symbols in lhs that weren't caught earlier.
    // result already holds every rhs symbol, so try_emplace only inserts lhs-only symbols.
    for (auto [symbol, coeff] : this->_coefficients) {
        result._coefficients.try_emplace(symbol, right._center * coeff);
    }

    // Affine multiplication adds a noise symbol.
    // For now, we add an error w/ coeff rad * rad, following Affapy impl.
    result._coefficients[new_noise_symbol()] = this->radius() * right.radius();
    if (result._coefficients.size() > MAX_NOISE_SYMBOLS) {
        result.collapse();
    }

#   ifdef AFFINE_TIME_MULT
    print_debug_info("multiplication", time);
#   endif

    return result;
}
AffineForm AffineForm::operator/(const AffineForm &right) const {
    return operator*(right.inv());
}
AffineForm AffineForm::pow(uint32_t power) const {
#   ifdef AFFINE_TIME_POW
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    // TODO: as we adapt the numeric API, we could switch this to use negative numbers w/ the inverse strategy.
    if (power == 0) {
        // Our implementation always returns affine forms, even if the power is 0 -- in this case, an exact affine form.
        return { 1, std::unordered_map<noise_symbol_t, double>() };
    }

    auto odd_power = power > 1 && power % 2 == 1;
    if (odd_power) {
        // Can descend logarithmically given an even power. We will do the extra mult later.
        power -= 1;
    }

    auto result = clone();
    while (power > 1) {
        // Perform multiply and half power each time until down to pow 1, unit operation.
        // Insight: squaring intermediate results allows our quick logarithmic descent.
        result = result * result;
        power /= 2;
    }

    if (odd_power) {
        // Now perform that last standard multiplication we saved.
        return result * *this;
    }

#   ifdef AFFINE_TIME_POW
    print_debug_info("multiplication", time);
#   endif

    return result;
}
AffineForm AffineForm::sqrt() const {
    auto interval = this->to_interval();
    auto min = interval.min();
    auto max = interval.max();

    if (min < 0) {
        // return NAN if negative value possible -- not defined for real semantics
        return AffineForm(NAN);
    }
    if (min == 0 && max == 0) {
        // special case: min, max both 0 -- we know value is 0.
        return AffineForm(0);
    }

    auto t = std::floor(std::sqrt(min) + std::sqrt(max));
    auto alpha = 1 / t;
    auto dzeta = t / 8 + 0.5 * std::sqrt(min * max) / t;

    auto rdelta = std::sqrt(max) - std::sqrt(min);
    auto delta = rdelta * rdelta / std::floor(8 * t);

    return approximate_affine_form(alpha, dzeta, delta);
}
AffineForm AffineForm::exp() const {
    auto rad = radius();
    if (rad == INFINITY || std::isnan(rad)) {
        return AffineForm(Winterval(0, INFINITY));
    }

    auto interval = Winterval(_center - rad, _center + rad);
    auto min = interval.min();
    auto max = interval.max();
    auto exp_min = std::exp(min);
    auto exp_max = std::exp(max);

    // Linear interpolation of (min, exp(min)) and (max, exp(max)) fns
    auto alpha = (exp_max - exp_min) / (max - min);
    auto log_alpha = std::log(alpha); // max point

    // zeta = offset of center
    auto zeta = alpha * (1 - log_alpha);

    // error of new noise symbol
    auto max_delta = alpha * (log_alpha - (1 - min)) + exp_min;
    auto delta = max_delta / 2; // Convert to radius, coeff of new noise symbol.

    return approximate_affine_form(alpha, zeta, delta);
}
AffineForm AffineForm::tanh() const {
    // use approximation from DeepZ
    auto interval = this->to_interval();
    auto tanh_min = std::tanh(interval.min());
    auto tanh_max = std::tanh(interval.max());

    auto derivative_min = 1 - std::pow(tanh_min, 2);
    auto derivative_max = 1 - std::pow(tanh_max, 2);

    auto min_deviation = std::min(derivative_min, derivative_max);

    auto offset = 0.5 * (tanh_max + tanh_min - min_deviation * (interval.max() + interval.min()));
    auto noise_coeff = 0.5 * (tanh_max - tanh_min - min_deviation * (interval.max() - interval.min()));

    auto result = *this * min_deviation;
    result._center += offset;
    result.add_noise_symbol(noise_coeff);
    return result;
}
AffineForm AffineForm::sigmoid() const {
    auto interval = this->to_interval();

    auto exp_min = std::exp(interval.min());
    auto exp_max = std::exp(interval.max());

    auto sigmoid_min = exp_min / (1 + exp_min);
    auto sigmoid_max = exp_max / (1 + exp_max);

    auto derivative_min = exp_min / std::pow(1 + exp_min, 2);
    auto derivative_max = exp_max / std::pow(1 + exp_max, 2);
    auto min_deviation = std::min(derivative_min, derivative_max);

    auto offset = 0.5 * (sigmoid_max + sigmoid_min - min_deviation * (interval.max() + interval.min()));
    auto noise_coeff = 0.5 * (sigmoid_max - sigmoid_min - min_deviation * (interval.max() - interval.min()));

    auto result = *this * min_deviation;
    result._center += offset;
    result.add_noise_symbol(noise_coeff);
    return result;
}
AffineForm AffineForm::relu() const {
    auto interval = this->to_interval();
    if (interval.min() > 0) {
        return clone();
    }
    if (interval.max() <= 0) {
        return AffineForm(0);
    }

    // Intersections derived form DeepZ
    auto upper = interval.max();
    auto lower = interval.min();

    auto intersection_point = upper / (upper - lower);
    auto error_rad = -0.5 * upper * lower / (upper - lower);

    auto result = *this * intersection_point;
    result._center += error_rad;
    result.add_noise_symbol(error_rad);
    return result;
}

/*
 * Compositional operations
 */
AffineForm AffineForm::union_with(const AffineForm &other) const {
    // using join formula from Taylor1+ (https://link.springer.com/chapter/10.1007/978-3-642-02658-4_47)
    auto result = clone();
    auto interval_union = to_interval().union_with(other.to_interval());
    result._center = interval_union.mid();

    for (auto [symbol, coeff] : other._coefficients) {
        // Add missing error terms with other's magnitude; otherwise, union takes the min of the two.
        // result starts as a clone of this, so an existing entry still holds this's coefficient.
        if (auto [it, inserted] = result._coefficients.try_emplace(symbol, coeff); !inserted) {
            it->second = std::min(std::abs(it->second), std::abs(coeff));
        }
    }

    auto error = interval_union.max() - result.center() - result.radius();
    if (error > 0) {
        result._coefficients[new_noise_symbol()] = error;
    }
    if (result._coefficients.size() > MAX_NOISE_SYMBOLS) {
        result.collapse();
    }

    return result;
}
std::vector<AffineForm> AffineForm::split(uint32_t n_splits) const {
    if (n_splits == 0) {
        throw std::invalid_argument("Cannot split into zero intervals.");
    }

    auto interval = to_interval();
    auto split_intervals = interval.split(n_splits);
    std::vector<AffineForm> result_forms;
    result_forms.reserve(n_splits);

    for (auto &split_interval : split_intervals) {
        result_forms.emplace_back(split_interval);
    }

    return result_forms;
}

/*
 * Binary affine comparison operators.
 */
bool AffineForm::operator==(const AffineForm &other) const {
    return _center == other._center && _coefficients == other._coefficients;
}
bool AffineForm::operator!=(const AffineForm &other) const {
    return !operator==(other);
}
bool AffineForm::operator<(const AffineForm &other) const {
    return to_interval() < other.to_interval();
}
bool AffineForm::operator<=(const AffineForm &other) const {
    return to_interval() <= other.to_interval();
}
bool AffineForm::operator>(const AffineForm &other) const {
    return to_interval() > other.to_interval();
}
bool AffineForm::operator>=(const AffineForm &other) const {
    return to_interval() >= other.to_interval();
}

/*
 * Scalar comparison operators.
 */
bool AffineForm::operator<(double other) const {
    return to_interval() < other;
}
bool AffineForm::operator>(double other) const {
    return to_interval() > other;
}
bool AffineForm::operator<=(double other) const {
    return to_interval() <= other;
}
bool AffineForm::operator>=(double other) const {
    return to_interval() >= other;
}

/*
 * Scalar arithmetic operators
 */
AffineForm AffineForm::operator*(double other) const {
    auto value = clone();
    value._center *= other;
    for (auto &coeff : value._coefficients | std::views::values) {
        coeff *= other;
    }
    return value;
}
AffineForm AffineForm::operator+(double other) const {
    auto value = clone();
    value._center += other;
    // Notice: addition does not affect error symbols.
    // effectively, the polytope is simply being translated.
    return value;
}
AffineForm AffineForm::operator-(double other) const {
    auto value = clone();
    value._center -= other;
    return value;
}
AffineForm AffineForm::operator/(double other) const {
    // Special case: 0. In affine forms, this will set all terms to 0, leading to a unit form.
    if (other == 0) {
        return { 0, std::unordered_map<noise_symbol_t, double>() };
    }
    return operator*(1 / other);
}

/*
 * Explicit modifiers
 */
void AffineForm::add_noise_symbol(double coeff) {
    _coefficients.insert({new_noise_symbol(), coeff});
    if (_coefficients.size() > MAX_NOISE_SYMBOLS) {
        collapse();
    }
}

/*
 * Accessors
 */
std::string AffineForm::to_string() const {
    std::string retval = std::string();
    retval += "Interval concretization: ";
    auto rad = radius();
    auto interval = Winterval(_center - rad, _center + rad);
    retval += "[" + std::to_string(interval.min());
    retval += ", ";
    retval += std::to_string(interval.max()) + "]\n";
    retval += "Center: " + std::to_string(_center) + "\n";
    retval += "Radius: " + std::to_string(rad) + "\n";
    retval += "Noise symbols:";
    for (auto [symbol, coeff] : _coefficients) {
        retval += " (" + std::to_string(symbol) + ": " + std::to_string(coeff) + "),";
    }
    retval += "\n";
    return retval;
}

double AffineForm::center() const {
    return _center;
}
double AffineForm::radius() const {
    return std::accumulate(_coefficients.begin(), _coefficients.end(), 0.0,
        [](auto sum, auto pair) { return sum + std::abs(pair.second); });
}

double AffineForm::min() const {
    return _center - radius();
}
double AffineForm::max() const {
    return _center + radius();
}
Winterval AffineForm::to_interval() const {
    // One pass over the coefficients rather than one each for min() and max().
    auto rad = radius();
    return {_center - rad, _center + rad};
}

double AffineForm::coeff_of(noise_symbol_t symbol) const {
    auto it = _coefficients.find(symbol);
    return it == _coefficients.end() ? NAN : it->second;
}

/*
 * Non-affine approximators
 */

/*

Approximating a non-affine form follows a general pattern:
- Create a new center with new center alpha * old_center + zeta
- Scale each error term by alpha
- Add a new error term with coeff delta.

*/
AffineForm AffineForm::approximate_affine_form(double alpha, double zeta, double delta) const {
#   ifdef AFFINE_TIME_OTHER
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    auto center = alpha * _center + zeta;

    auto map = std::unordered_map<noise_symbol_t, double>();
    for (auto [symbol, coeff] : _coefficients) {
        auto new_value = alpha * coeff;
        map.insert({symbol, new_value});
    }

    map.insert({new_noise_symbol(), delta});

#   ifdef AFFINE_TIME_OTHER
    print_debug_info("approximation constructions", time);
#   endif

    return { center, std::move(map) };
}

/*
 * Use mini-range approximation rather than standard Chebyshev.
 * Credit: libaffa
 */
AffineForm AffineForm::inv() const {
#   ifdef AFFINE_TIME_INV
    auto time = std::chrono::high_resolution_clock::now();
#   endif

    auto interval = this->to_interval();
    if (interval.contains(0)) {
        // If interval contains 0, infinity will be in this interval, or the interval will just be [0, 0].
        // The blowup to infinity wipes away the dependence of the variables and adds a term of unlimited magnitude.
        auto retval = clone();
        retval.add_noise_symbol(INFINITY);
        return retval;
    }

    auto a = interval.abs().min();
    auto b = interval.abs().max();
    // Derivative of 1/x: -1/x^2. Notice: series defined over first derivative.
    auto alpha = -1 / std::pow(b, 2);

    auto range = Winterval(1/a - alpha * a, 2 / b);
    auto zeta = range.mid();

    // if negative value included in interval, flip result.
    if (interval.min() < 0) {
        zeta = -zeta;
    }

#   ifdef AFFINE_TIME_INV
    print_debug_info("multiplication", time);
#   endif

    // New noise term will be radius of mini-range interval.
    return approximate_affine_form(alpha, zeta, range.radius());
}

/*
 * Internal helpers
 */
// Internal clone constructor
AffineForm AffineForm::clone() const {
    return { this->_center, this->_coefficients };
}

void AffineForm::collapse() {
    // Approximate by collapsing all noise symbols into one.
    double error_magnitude = 0;
    for (auto coeff: _coefficients | std::views::values) {
        error_magnitude += std::abs(coeff);
    }

    _coefficients = std::unordered_map<noise_symbol_t, double>();
    _coefficients.insert({new_noise_symbol(), error_magnitude});
}


/*
 * Associated operators.
 */
std::ostream& operator<<(std::ostream &os, const AffineForm &rhs) {
    os << "a" << rhs.to_interval();
    return os;
}