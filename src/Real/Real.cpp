//
// Created by will on 9/17/25.
//

#include <cmath>

#include "Real.hpp"

#include <iosfwd>
#include <string>

Real::Real(double value): _value(value) {}
Real::~Real() = default;

double Real::value() const {
    return _value;
}

/*
 * Abstraction accessors
 */
double Real::radius() const {
    return 0;
}
double Real::min() const {
    return _value;
}
double Real::max() const {
    return _value;
}

/*
 * Operations
 */
Real Real::operator+(const Real &right) const {
    return { _value + right._value };
}
Real Real::operator-(const Real &right) const {
    return { _value - right._value };
}
Real Real::operator*(const Real &right) const {
    return { _value * right._value };
}
Real Real::operator/(const Real &right) const {
    return { _value / right._value };
}
bool Real::operator==(const Real &right) const {
    return _value == right._value;
}
bool Real::operator<=(const Real &right) const {
    return _value <= right._value;
}
bool Real::operator<(const Real &right) const {
    return _value < right._value;
}
bool Real::operator>=(const Real &right) const {
    return _value >= right._value;
}
bool Real::operator>(const Real &right) const {
    return _value > right._value;
}
/*
 * Compositional operations
 */
Real Real::union_with(const Real &right) const {
    if (_value != right._value) {
        throw std::runtime_error("Cannot perform union on distinct discrete numbers.");
    }
    return _value;
}
std::vector<Real> Real::split(uint32_t parts) const {
    if (parts == 0) {
        throw std::invalid_argument("Number of parts to split must be positive!");
    }
    return std::vector(parts, *this);
}


std::ostream& operator<<(std::ostream &os, Real rhs) {
    os << std::to_string(rhs.value());
    return os;
}

Real Real::exp() const {
    return { std::exp(_value) };
}
Real Real::tanh() const {
    return { std::tanh(_value) };
}
Real Real::sigmoid() const {
    return { 1 / (1 + std::exp(-_value)) };
}
Real Real::pow(uint32_t power) const {
    return { std::pow(_value, power) };
}
Real Real::abs() const {
    return { std::abs(_value) };
}