//
// Created by will on 10/8/25.
//

#include "MixedForm.hpp"

#include <complex>
#include <memory>

/*
 * Constructors
 */
MixedForm::MixedForm(const AffineForm &affine_rep, const Winterval &interval_rep) :
    _affine_rep(affine_rep),
    _intersected_bounds(interval_intersection(affine_rep, interval_rep)) {}

MixedForm::MixedForm(const AffineForm &affine_rep) :
    _affine_rep(affine_rep),
    _intersected_bounds(affine_rep.to_interval()) {}

MixedForm::MixedForm(const Winterval &interval_rep) :
    _affine_rep(interval_rep),
    _intersected_bounds(interval_rep) {}

MixedForm::MixedForm(double value):
    _affine_rep(value),
    _intersected_bounds(value) {}

MixedForm::MixedForm() :
    _affine_rep(AffineForm()),
    _intersected_bounds(Winterval()) {}

MixedForm::~MixedForm() = default;

MixedForm MixedForm::union_with(const MixedForm &w) const {
    // Loses noise symbol dependence, but this would happen anyways.
    return MixedForm(interval_bounds().union_with(w.interval_bounds()));
}

std::vector<MixedForm> MixedForm::split(uint32_t num_splits) const {
    std::vector results = std::vector<MixedForm>();
    results.reserve(num_splits);
    // Since the affine form cannot be meaningfully split while preserving noise symbols,
    // we only split over the intersected bounds.
    auto interval_splits = _intersected_bounds.split(num_splits);
    for (const Winterval& split : interval_splits) {
        results.emplace_back(split);
    }
    return results;
}


/*
 * Accessors
 */
const AffineForm &MixedForm::affine_rep() const {
    return _affine_rep;
}
const Winterval &MixedForm::interval_bounds() const {
    return _intersected_bounds;
}
double MixedForm::radius() const {
    return _intersected_bounds.radius();
}
double MixedForm::min() const {
    return _intersected_bounds.min();
}
double MixedForm::max() const {
    return _intersected_bounds.max();
}

/*
 * Scalar operators
 */
MixedForm MixedForm::operator+(const double scalar) const {
    // Must propagate affine form, even if interval tighter, to preserve relationship between vars.
    return { _affine_rep + scalar, _intersected_bounds + scalar };
}
MixedForm MixedForm::operator-(const double scalar) const {
    // Must propagate affine form, even if interval tighter, to preserve relationship between vars.
    return { _affine_rep - scalar, _intersected_bounds - scalar };
}
MixedForm MixedForm::operator*(const double scalar) const {
    // Must propagate affine form, even if interval tighter, to preserve relationship between vars.
    return {
        _affine_rep * scalar,
        interval_intersection(_affine_rep * scalar, _intersected_bounds * scalar)
    };
}
MixedForm MixedForm::operator/(const double scalar) const {
    // Divison by 0 ill defined for affine forms.
    if (scalar == 0) {
        return MixedForm(_intersected_bounds / scalar);
    }

    // Otherwise, propagate affine form, even if interval tighter, to preserve relationship between vars.
    return {
        _affine_rep / scalar,
        interval_intersection(_affine_rep / scalar, _intersected_bounds / scalar)
    };
}

/*
 * Binary mixed relational operations.
 */
bool MixedForm::operator==(const MixedForm &rhs) const {
    return _affine_rep == rhs._affine_rep && _intersected_bounds == rhs._intersected_bounds;
}
bool MixedForm::operator!=(const MixedForm &rhs) const {
    return !operator==(rhs);
}
bool MixedForm::operator>(const MixedForm &right) const {
    return _intersected_bounds > right._intersected_bounds;
}
bool MixedForm::operator>=(const MixedForm &right) const {
    return _intersected_bounds >= right._intersected_bounds;
}
bool MixedForm::operator<(const MixedForm &right) const {
    return _intersected_bounds < right._intersected_bounds;
}
bool MixedForm::operator<=(const MixedForm &right) const {
    return _intersected_bounds <= right._intersected_bounds;
}

/*
 * Scalar relational operations.
 */
bool MixedForm::operator>(double scalar) const {
    return _intersected_bounds > scalar;
}
bool MixedForm::operator>=(double scalar) const {
    return _intersected_bounds >= scalar;
}
bool MixedForm::operator<(double scalar) const {
    return _intersected_bounds < scalar;
}
bool MixedForm::operator<=(double scalar) const {
    return _intersected_bounds <= scalar;
}

/*
 * Wixed-Wixed operations
 */
MixedForm MixedForm::operator-(const MixedForm &right) const {
    return { _affine_rep - right._affine_rep, _intersected_bounds - right._intersected_bounds };
}
MixedForm MixedForm::operator+(const MixedForm &w) const {
    return { _affine_rep + w._affine_rep, _intersected_bounds + w._intersected_bounds };
}
MixedForm MixedForm::operator*(const MixedForm &w) const {
    return { _affine_rep * w._affine_rep, _intersected_bounds * w._intersected_bounds };
}
MixedForm MixedForm::operator/(const MixedForm &w) const {
    return { _affine_rep / w._affine_rep, _intersected_bounds / w._intersected_bounds };
}

/*
 * Unary wixed operations
 */
MixedForm MixedForm::abs() const {
    return { _affine_rep.abs(), _intersected_bounds.abs() };
}
MixedForm MixedForm::pow(uint32_t pow) const {
    // TODO: fix powers to be only unsigned.
    return { _affine_rep.pow(pow), _intersected_bounds.pow(pow) };
}
MixedForm MixedForm::exp() const {
    return { _affine_rep.exp(), _intersected_bounds.exp() };
}

/*
 * Internal helpers
 */
Winterval MixedForm::interval_intersection(const AffineForm &a, const Winterval &b) {
    auto min_intersect = std::max(a.to_interval().min(), b.min());
    auto max_intersect = std::min(a.to_interval().max(), b.max());

    if (min_intersect > max_intersect) {
        // This can happen e.g. when dividing by 0.
        max_intersect = min_intersect;
    }

    return {min_intersect, max_intersect};
}

std::ostream& operator<<(std::ostream &os, const MixedForm &rhs) {
    os << "m" << rhs.interval_bounds();
    return os;
}