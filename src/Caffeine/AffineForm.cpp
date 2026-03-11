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
#include <vector>

/*
 * Noise symbol management
 */
AffineForm::noise_symbol_t AffineForm::max_noise_symbol = 0;

AffineForm::noise_symbol_t AffineForm::new_noise_symbol() {
    noise_symbol_t symbol;
#   pragma omp atomic capture
    symbol = max_noise_symbol++;
    return symbol;
}

void print_debug_info(const std::string &op_name, const std::chrono::high_resolution_clock::time_point &start_time) {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    std::cout << "[DEBUG] Operation " << op_name << " took " << duration << " microseconds." << std::endl;
}

/*
 * Constructors
 */
AffineForm::AffineForm(double center, const std::unordered_map<noise_symbol_t, double> &starting_coeffs):
    _center(center), _coefficients(starting_coeffs) {

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
AffineForm::AffineForm(): _center(0), _coefficients(std::unordered_map<noise_symbol_t, double>()) {}

/*
 * Unary operators
 */
AffineForm AffineForm::operator-() const {
    auto value = clone();
    value._center = -_center;
    for (const auto symbol: _coefficients | std::views::keys) {
        value._coefficients[symbol] *= -1;
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
    for (auto symbol: value._coefficients | std::views::keys) {
        value._coefficients[symbol] /= 2;
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
        if (!value._coefficients.contains(symbol)) {
            value._coefficients[symbol] = coeff;
        } else {
            value._coefficients[symbol] += coeff;
        }
    }

#   ifdef AFFINE_TIME_ADD
    print_debug_info("addition", time);
#   endif

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
        if (!value._coefficients.contains(symbol)) {
            value._coefficients[symbol] = -coeff;
        } else {
            value._coefficients[symbol] -= coeff;
        }
    }

#   ifdef AFFINE_TIME_SUB
    print_debug_info("subtraction", time);
#   endif

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
        if (!this->_coefficients.contains(symbol)) {
            // Add missing error terms scaled by left's center.
            result._coefficients[symbol] = this->_center * coeff;
        } else if (this->_coefficients.contains(symbol)) {
            auto a = this->_center * coeff;
            auto b = right._center * this->_coefficients.at(symbol);
            result._coefficients[symbol] = a + b;
        }
    }

    // Now go and perform similar calculation for error symbols in lhs that weren't caught earlier.
    for (auto [symbol, coeff] : this->_coefficients) {
        if (!right._coefficients.contains(symbol)) {
            result._coefficients[symbol] = right._center * coeff;
        }
    }

    // Affine multiplication adds a noise symbol.
    // For now, we add an error w/ coeff rad * rad, following Affapy impl.
    result._coefficients[new_noise_symbol()] = this->radius() * right.radius();

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
AffineForm AffineForm::union_with(const AffineForm &other) const {
    auto interval1 = this->to_interval();
    auto interval2 = other.to_interval();
    auto union_interval = interval1.union_with(interval2);
    // Notice the loss of dependence information here.
    return AffineForm(union_interval);
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
    for (const auto symbol: _coefficients | std::views::keys) {
        value._coefficients[symbol] *= other;
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
}

/*
 * Accessors
 */
std::string AffineForm::to_string() const {
    std::string retval = std::string();
    retval += "Interval concretization: ";
    retval += "[" + std::to_string(to_interval().min());
    retval += ", ";
    retval += std::to_string(to_interval().max()) + "]\n";
    retval += "Center: " + std::to_string(_center) + "\n";
    retval += "Radius: " + std::to_string(radius()) + "\n";
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
    return to_interval().min();
}
double AffineForm::max() const {
    return to_interval().max();
}

double AffineForm::coeff_of(noise_symbol_t symbol) const {
    if (_coefficients.contains(symbol)) {
        return _coefficients.at(symbol);
    }
    return NAN;
}

Winterval AffineForm::to_interval() const {
    // Note: unable to do accumulation because of behavior with unordered maps.
    const double rad = radius();

    auto min = _center - rad;
    auto max = _center + rad;
    if (std::isnan(min)) {
        min = INFINITY * -1;
    }
    if (std::isnan(max)) {
        max = INFINITY;
    }

    return {min, max};
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

    return { center, map };
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
        // We will hence construct a new affine form from the corresponding interval, since no dependencence can be preserved.
        return AffineForm(this->to_interval());
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
std::ostream& operator<<(std::ostream &os, AffineForm rhs) {
    os << "a" << rhs.to_interval();
    return os;
}