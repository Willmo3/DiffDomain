//
// Created by will on 9/17/25.
//

#ifndef PDENCLOSE_REAL_H
#define PDENCLOSE_REAL_H
#include <cstdint>
#include <iosfwd>

#include "cereal/archives/binary.hpp"


/**
 * Wrapper class for real numbers, compliant with template requirements of abstract domain.
 */
class Real {
public:
    Real() = default;
    Real(double value);
    ~Real();

    /*
     * Accessors
     */
    double value() const;

    /*
     * Operations
     */
    Real operator+(const Real &right) const;
    Real operator-(const Real &right) const;
    Real operator*(const Real &right) const;
    Real operator/(const Real &right) const;
    bool operator==(const Real &right) const;
    bool operator<=(const Real &right) const;
    bool operator<(const Real &right) const;
    bool operator>=(const Real &right) const;
    bool operator>(const Real &right) const;

    Real exp() const;
    Real tanh() const;
    Real sigmoid() const;
    Real relu() const;
    Real pow(uint32_t power) const;
    Real abs() const;

    /*
     * Abstraction accessors
     * Tombstones for numeric abstract domain
     */
    double radius() const;
    double min() const;
    double max() const;

    /*
     * Compositional operations
     * Tombstones for numeric abstract domain
     */
    Real union_with(const Real &right) const;
    std::vector<Real> split(uint32_t parts) const;

    /*
     * Serialization support through cereal.
     */
    template<class Archive>
    void serialize(Archive &archive) {
        archive(_value);
    }
private:
    double _value;
};

// Note: cannot use reference for rhs because we want to be able to print shortlived values
// i.e. std::cout << Real(a) + Real(b) << std::endl;
std::ostream& operator<<(std::ostream &os, Real rhs);


#endif //PDENCLOSE_REAL_H