//
// Created by will on 10/8/25.
//

#ifndef DUALDOMAIN_MIXEDFORM_H
#define DUALDOMAIN_MIXEDFORM_H

#include <memory>

#include "Caffeine/AffineForm.hpp"
#include "Winterval/Winterval.hpp"

/**
 * Mixed affine form. Tracks both interval and affine representation of data.
 * Credit to Phocus library.
 *
 * @author Will Morris
 */
class MixedForm {
public:
    /*
     * Constructors
     */

    /**
     * Create a new MixedForm from an affine form and an interval.
     * @param affine_rep Affine form to construct from.
     * @param interval_rep Interval to construct from.
     */
    MixedForm(const AffineForm &affine_rep, const Winterval &interval_rep);

    /**
     * Create a new MixedForm from an affine form. Interval is constructed from affine form.
     * @param affine_rep Affine form to construct from.
     */
    explicit MixedForm(const AffineForm &affine_rep);

    /**
     * Create a new MixedForm from an interval. Affine form is constructed from interval.
     * @param interval_rep Interval to construct from.
     */
    explicit MixedForm(const Winterval &interval_rep);

    /**
     * @param value Exact value to center affine and interval representations around
     */
    explicit MixedForm(double value);

    /**
     * Default constructor, useful for serialization.
     * Create an empty MixedForm by calling default constructor on subcomponents.
     */
    MixedForm();

    /**
     * Destructor
     */
    ~MixedForm();

    /*
     * Accessors
     */
    const AffineForm &affine_rep() const;
    const Winterval &interval_bounds() const;
    double min() const;
    double max() const;
    /**
     * @return the radius of the intersection of interval and affine bounds.
     */
    double radius() const;

    /*
     * Wixed-Wixed operations
     */
    MixedForm operator+(const MixedForm &w) const;
    MixedForm operator-(const MixedForm &w) const;
    MixedForm operator*(const MixedForm &w) const;
    MixedForm operator/(const MixedForm &w) const;

    /*
     * Compositional operations
     */
    /**
     * @param w Mixed form to union with this one.
     * @return A new MixedForm representing the union of this and w.
     * Since affine noise symbols will be deleted by unions, we only preserve the interval representation
     * and construct a new affine form from the interval union.
     */
    MixedForm union_with(const MixedForm &w) const;

    /**
     *
     * @param num_splits Number of splits to perform
     * @return A vector of new mixed forms over the intersected interval bounds of this form.
     */
    std::vector<MixedForm> split(uint32_t num_splits) const;

    /*
     * Unary Wixed operations
     */
    MixedForm pow(uint32_t power) const;
    /**
     * @return Absolute value of interval bounds of this form.
     * Affine form is not considered -- abs is poorly defined for affine forms.
     */
    MixedForm abs() const;

    /*
     * Wixed-scalar operations
     */
    MixedForm operator+(double scalar) const;
    MixedForm operator-(double scalar) const;
    MixedForm operator*(double scalar) const;
    MixedForm operator/(double scalar) const;

    /*
     * Binary comparison operations
     */
    /**
     * @return whether rhs has the same affine and interval representation as this form.
     */
    bool operator==(const MixedForm &rhs) const;
    /**
     * @return whether rhs has a different affine or interval representation from this form.
     */
    bool operator!=(const MixedForm &rhs) const;
    /**
     * @return Whether the intersected bounds of this form are strictly less than those of rhs.
     */
    bool operator<(const MixedForm &rhs) const;
    /**
     * @return Whether the intersected bounds of this form are less than or equal to those of rhs.
     */
    bool operator<=(const MixedForm &rhs) const;
    /**
     * @return Whether the intersected bounds of this form are strictly greater than those of rhs.
     */
    bool operator>(const MixedForm &rhs) const;
    /**
     * @return Whether the intersected bounds of this form are greater than or equal to those of rhs.
     */
    bool operator>=(const MixedForm &rhs) const;

    /*
     * Scalar relational operations
     */
    bool operator<(double scalar) const;
    bool operator<=(double scalar) const;
    bool operator>(double scalar) const;
    bool operator>=(double scalar) const;

private:
    /*
     * Internal helpers.
     */
    /**
     * @param a Affine form to intersect
     * @param b Interval to intersect
     * @return An interval enclosing the intersection of these two forms.
     */
    static Winterval interval_intersection(const AffineForm &a, const Winterval &b);

    /*
     * Private fields
     * Note: we use values here because large portions of the object are in memory.
     */
    AffineForm _affine_rep;
    Winterval _intersected_bounds;
};

// Utility pipe to stdout
std::ostream& operator<<(std::ostream &os, const MixedForm &rhs);


#endif //DUALDOMAIN_MIXEDFORM_H