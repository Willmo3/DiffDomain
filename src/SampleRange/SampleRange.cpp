//
// Created by will on 3/26/26.
//

#include "SampleRange.hpp"

#include <cmath>
#include <utility>

/*
 * Constructors
 */

SampleRange::SampleRange() : _samples(Eigen::VectorXd::Zero(SAMPLED_VALUES)) {}

SampleRange::SampleRange(double value) : _samples(Eigen::VectorXd::Constant(SAMPLED_VALUES, value)) {}

SampleRange::SampleRange(double min, double max) {
    // Generate random values in [min, max] using Eigen's Random() method
    Eigen::ArrayXd t = (Eigen::ArrayXd::Random(SAMPLED_VALUES) + 1.0) / 2.0;
    // Samples 0 and 1 are pinned to min and max.
    t(0) = 0.0;
    t(1) = 1.0;
    // overflow protection
    const auto width = max - min;
    if (std::isfinite(width)) {
        _samples = min + width * t;
    } else {
        _samples = min * (1.0 - t) + max * t;
    }
}

SampleRange::SampleRange(Eigen::VectorXd samples) : _samples(std::move(samples)) {
    if (_samples.size() != static_cast<int>(SAMPLED_VALUES)) {
        throw std::runtime_error("SampleRange samples must have the same number of values");
    }
}

/*
 * Accessors
 */
double SampleRange::min() const {
    return _samples.minCoeff();
}
double SampleRange::max() const {
    return _samples.maxCoeff();
}
double SampleRange::radius() const {
    // Halve before subtracting: max - min overflows for large bounds.
    return max() / 2.0 - min() / 2.0;
}

/*
 * SampleRange-SampleRange arithmetic
 */

SampleRange SampleRange::operator+(const SampleRange &rhs) const {
    return { _samples + rhs._samples };
}
SampleRange SampleRange::operator-(const SampleRange &rhs) const {
    return { _samples - rhs._samples };
}
SampleRange SampleRange::operator*(const SampleRange &rhs) const {
    return { _samples.cwiseProduct(rhs._samples) };
}
SampleRange SampleRange::operator/(const SampleRange &rhs) const {
    return { _samples.cwiseQuotient(rhs._samples) };
}

/*
 * SampleRange-scalar arithmetic
 */

SampleRange SampleRange::operator+(double scalar) const {
    return { _samples.array() + scalar };
}
SampleRange SampleRange::operator-(double scalar) const {
    return { _samples.array() - scalar };
}
SampleRange SampleRange::operator*(double scalar) const {
    return { _samples * scalar };
}
SampleRange SampleRange::operator/(double scalar) const {
    return { _samples / scalar };
}

/*
 * Unary operations
 */

SampleRange SampleRange::pow(uint32_t n) const {
    if (n == 0) {
        return SampleRange(1.0);
    }
    auto result = _samples.array().pow(static_cast<double>(n));
    return { result };
}
SampleRange SampleRange::sqrt() const {
    return { _samples.array().sqrt() };
}
SampleRange SampleRange::abs() const {
    return { _samples.array().abs() };
}
SampleRange SampleRange::exp() const {
    return { _samples.array().exp() };
}
SampleRange SampleRange::tanh() const {
    return { _samples.array().tanh() };
}
SampleRange SampleRange::sigmoid() const {
    return { 1.0 / (1.0 + (-_samples.array()).exp()) };
}
SampleRange SampleRange::relu() const {
    return { _samples.array().max(0.0) };
}

/*
 * Comparison operators
 */

bool SampleRange::operator==(const SampleRange &rhs) const {
    return _samples.isApprox(rhs._samples, 1e-10);
}
bool SampleRange::operator!=(const SampleRange &rhs) const {
    return !operator==(rhs);
}
bool SampleRange::operator<(const SampleRange &rhs) const {
    return max() < rhs.min();
}
bool SampleRange::operator<=(const SampleRange &rhs) const {
    return max() <= rhs.min();
}
bool SampleRange::operator>(const SampleRange &rhs) const {
    return min() > rhs.max();
}
bool SampleRange::operator>=(const SampleRange &rhs) const {
    return min() >= rhs.max();
}

/*
 * SampleRange-scalar comparison operators
 */

bool SampleRange::operator<(double scalar) const {
    return max() < scalar;
}
bool SampleRange::operator<=(double scalar) const {
    return max() <= scalar;
}
bool SampleRange::operator>(double scalar) const {
    return min() > scalar;
}
bool SampleRange::operator>=(double scalar) const {
    return min() >= scalar;
}

/*
 * Compositional operations
 */
// TODO: reason about how to use union and split over sample ranges
SampleRange SampleRange::union_with(const SampleRange &rhs) const {
    auto union_min = std::min(min(), rhs.min());
    auto union_max = std::max(max(), rhs.max());
    return { union_min, union_max };
}
std::vector<SampleRange> SampleRange::split(uint32_t n_splits) const {
    if (n_splits == 0) {
        throw std::invalid_argument("Attempted to split SampleRange into 0 subranges!");
    }
    std::vector<SampleRange> result;
    result.reserve(n_splits);

    auto min_sample = min();
    auto max_sample = max();
    // Divide before subtracting: max - min overflows for large bounds.
    auto step = max_sample / n_splits - min_sample / n_splits;

    for (auto i = 0; i < n_splits; ++i) {
        auto split_min = min_sample + i * step;
        auto split_max = (i == n_splits - 1) ? max_sample : min_sample + (i + 1) * step;
        result.emplace_back(split_min, split_max);
    }

    return result;
}

SampleRange SampleRange::pointwise_max(const SampleRange &rhs) const {
    return { _samples.cwiseMax(rhs._samples) };
}

/*
 * Utility functions
 */

std::ostream &operator<<(std::ostream &os, const SampleRange &sr) {
    os << "s" << "[" << sr.min() << ", " << sr.max() << "]";
    return os;
}

