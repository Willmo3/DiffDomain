//
// Created by will on 3/26/26.
//

#ifndef DIFFDOMAIN_SAMPLERANGE_HPP
#define DIFFDOMAIN_SAMPLERANGE_HPP

#include <cstdint>
#include <ostream>
#include <vector>

#include "Eigen/Dense"

/**
 * Numeric domain representing a range of real values sampled from a distribution.
 * @author Will Morris
 */
class SampleRange final {
public:
    // Change this value to sample more or less data
    static constexpr uint32_t SAMPLED_VALUES = 2048;

    /*
     * Constructors
     */

    SampleRange();
    /**
     * @param value The scalar value to initialize all samples with.
     */
    explicit SampleRange(double value);
    /**
     * Construct a SampleRange uniformly sampled between min and max.
     * @param min Minimum value of the range.
     * @param max Maximum value of the range.
     */
    SampleRange(double min, double max);
    /**
     * @param samples Eigen vector of samples. Size must equal SAMPLED_VALUES
     */
    SampleRange(Eigen::VectorXd samples);

    /*
     * Accessors
     */

    /**
     * @return The minimum value across all samples.
     */
    [[nodiscard]] double min() const;
    /**
     * @return The maximum value across all samples.
     */
    [[nodiscard]] double max() const;
    /**
     * @return The radius (half-width) of the range, i.e., (max - min) / 2.
     */
    [[nodiscard]] double radius() const;

    /*
     * SampleRange-SampleRange arithmetic
     */
    SampleRange operator+(const SampleRange &rhs) const;
    SampleRange operator-(const SampleRange &rhs) const;
    SampleRange operator*(const SampleRange &rhs) const;
    SampleRange operator/(const SampleRange &rhs) const;

    /*
     * SampleRange-scalar arithmetic
     */
    SampleRange operator+(double scalar) const;
    SampleRange operator-(double scalar) const;
    SampleRange operator*(double scalar) const;
    SampleRange operator/(double scalar) const;

    /*
     * Unary operations
     */

    /**
     * @param n Non-negative integer exponent.
     * @return A new SampleRange with each element raised to power n.
     */
    [[nodiscard]] SampleRange pow(uint32_t n) const;
    /**
     * @return A new SampleRange with absolute value of each sample.
     */
    [[nodiscard]] SampleRange abs() const;
    /**
     * @return A new SampleRange with the exponential (e^x) of each sample.
     */
    [[nodiscard]] SampleRange exp() const;

    /*
     * Comparison operators
     */
    bool operator==(const SampleRange &rhs) const;
    bool operator!=(const SampleRange &rhs) const;

    /**
     * Strictly less than => max value of this is less than min of rhs
     */
    bool operator<(const SampleRange &rhs) const;
    /**
     * Leq => max value of this is leq min of rhs
     */
    bool operator<=(const SampleRange &rhs) const;
    /**
     * Greater than => min value of this is greater than max of rhs
     */
    bool operator>(const SampleRange &rhs) const;
    /**
     * Greater-than-or-equal => min value of this is geq max of rhs
     */
    bool operator>=(const SampleRange &rhs) const;

    /*
     * SampleRange-scalar comparison operators
     */
    bool operator<(double scalar) const;
    bool operator<=(double scalar) const;
    bool operator>(double scalar) const;
    bool operator>=(double scalar) const;

    /*
     * Compositional operations
     */

    /**
     * Union with another SampleRange.
     * Returns a new SampleRange with the outer bounds (min of mins, max of maxes)
     * sampled uniformly.
     *
     * @param rhs Other SampleRange to union with.
     * @return A new SampleRange spanning the union bounds.
     */
    [[nodiscard]] SampleRange union_with(const SampleRange &rhs) const;

    /**
     * Split this SampleRange into n_splits roughly equal subranges.
     * Each split is sampled uniformly between its respective min and max bounds.
     *
     * @param n_splits Number of splits to produce.
     * @return Vector of n_splits SampleRanges.
     */
    [[nodiscard]] std::vector<SampleRange> split(uint32_t n_splits) const;

    /*
     * Utility functions
     */
    friend std::ostream &operator<<(std::ostream &os, const SampleRange &sr);

private:
    Eigen::VectorXd _samples;
};

#endif // DIFFDOMAIN_SAMPLERANGE_HPP