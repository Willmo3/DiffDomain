//
// Created by will on 3/26/26.
//

#include "SampleRange.hpp"

#include <utility>

/*
 * Constructors
 */

SampleRange::SampleRange() : _samples(Eigen::VectorXd::Zero(SAMPLED_VALUES)) {}

SampleRange::SampleRange(double value) : _samples(Eigen::VectorXd::Constant(SAMPLED_VALUES, value)) {}

SampleRange::SampleRange(double min, double max) {
    // Generate random values in [min, max] using Eigen's Random() method
    // Random() generates values in [-1, 1], so scale to [0, 1] then to [min, max]
    _samples = min + (max - min) * (Eigen::VectorXd::Random(SAMPLED_VALUES).array() + 1.0) / 2.0;
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
    return (max() - min()) / 2.0;
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

SampleRange SampleRange::abs() const {
    return { _samples.array().abs() };
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
    std::vector<SampleRange> result;
    result.reserve(n_splits);

    auto min_sample = min();
    auto max_sample = max();
    auto step = (max_sample - min_sample) / n_splits;

    for (auto i = 0; i < n_splits; ++i) {
        auto split_min = min_sample + i * step;
        auto split_max = min_sample + (i + 1) * step;
        result.emplace_back(split_min, split_max);
    }

    return result;
}

/*
 * Utility functions
 */

std::ostream &operator<<(std::ostream &os, const SampleRange &sr) {
    os << "s" << "[" << sr.min() << ", " << sr.max() << "]";
    return os;
}

