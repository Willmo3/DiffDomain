//
// Created by will on 9/15/25.
//

#ifndef WINTERVAL_WINTERVAL_H
#define WINTERVAL_WINTERVAL_H
#include <string>
// Note: cereal root must be in the build path
#include "cereal/archives/json.hpp"

/**
 * @brief Intervals defined with support for operator overloading.
 * Intended for use in simple numeric abstract interpretation.
 *
 * @author Will Morris
 */
class Winterval {
public:
    /**
     * Default constructor for winterval with values [0,0].
     */
    Winterval();

    /**
     * @param min minimum value in interval, inclusive.
     * @param max Maximum value in interval, inclusive.
     */
    Winterval(double min, double max);
    /**
     * @param value Value to form an exact bound around.
     */
    Winterval(double value);
    ~Winterval();

    /*
     * Accessors
     */
    double min() const;
    double max() const;
    /**
     * @return the midpoint of this interval.
     */
    double mid() const;
    /**
     * @return the distance of the min and max to the midpoint.
     */
    double radius() const;

    /*
     * Binary Winterval operations.
     * Note: in order to enable nested operator application, we don't use references.
     * Overhead should be small, but if we need to squeeze out some cycles later, look here!
     */
    Winterval operator+(const Winterval &rhs) const;
    Winterval operator-(const Winterval &rhs) const;
    Winterval operator*(const Winterval &rhs) const;
    Winterval operator/(const Winterval &rhs) const;

    /*
     * Compositional operations
     */

    /**
     * @param rhs Another Winterval to union with this one.
     * @return A new Winterval representing the union of this and rhs.
     */
    Winterval union_with(const Winterval &rhs) const;

    /**
     * @param n_splits Number of splits to perform over this interval.
     * @return An array of subintervals which, when composed, equal this Winterval.
     */
    std::vector<Winterval> split(uint32_t n_splits) const;

    /*
     * Binary relational operations.
     */
    bool operator==(const Winterval &rhs) const;
    /**
     * @return Whether every element in this interval is <= every element in rhs.
     */
    bool operator<=(const Winterval &rhs) const;
    /**
     * @return Whether every element in this interval is < every element in rhs.
     */
    bool operator<(const Winterval &rhs) const;
    /**
     * @return Whether every element in this interval is >= every element in rhs.
     */
    bool operator>=(const Winterval &rhs) const;
    /**
     * @return Whether every element in this interval is > every element in rhs.
     */
    bool operator>(const Winterval &rhs) const;

    /*
     * Scalar comparison operators
     */
    // Evaluated against max bound.
    bool operator<=(double rhs) const;
    bool operator<(double rhs) const;
    // Evaluated against min bound.
    bool operator>=(double rhs) const;
    bool operator>(double rhs) const;

    /*
     * Scalar Winterval operations.
     */
    Winterval operator+(double rhs) const;
    Winterval operator-(double rhs) const;
    Winterval operator*(double rhs) const;
    Winterval operator/(double rhs) const;

    /*
     * Other numeric functions.
     */
    // TODO: as more flux functions are needed, add them here.
    /**
     * Compute the hyperbolic tangent of an interval.
     * Since tanh is monotonic, it is sufficient to apply it to the lower, upper bounds (i.e. optimal conditions on boundary)
     * @return A new interval whose value is the hyperbolic tangent of this interval.
     */
    Winterval tanh() const;

    /**
     * Compute an interval to some constant power.
     * @param power Power to raise the interval to.
     * Note that for intervals, integer powers should always be nonnegative -- we represent as an int to comply with idiomatic Cpp.
     *
     * @return a new interval equal to this one raised to power
     */
    Winterval pow(int power) const;

    /**
     * @return The absolute value of this interval.
     */
    Winterval abs() const;

    /**
     * @brief Inclusive bounds check.
     * @param value Value to check if in interval
     * @return Whether value is in the interval.
     */
    bool contains(double value) const;

    /*
     * Serialization support through cereal.
     */
    template<class Archive>
    void serialize(Archive &archive) {
        archive( cereal::make_nvp("min", _min), cereal::make_nvp("max", _max) );
    }

private:
    double _min;
    double _max;
};

std::ostream& operator<<(std::ostream &os, Winterval rhs);

#endif //WINTERVAL_WINTERVAL_H