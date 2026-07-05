//
// Created by will on 6/9/26.
//

#ifndef PDENCLOSE_ROOTFINDING_HPP
#define PDENCLOSE_ROOTFINDING_HPP
#include <vector>

/**
 *
 * @param coeffs Coefficients of the polynomial, in order of decreasing degree. For example, for 2x^3 + 3x^2 + x + 5, the coeffs would be {2, 3, 1, 5}.
 * @return A vector of complex-valued roots.
 */
std::vector<std::complex<double>> poly_roots(const std::vector<double> &coeffs);

#endif //PDENCLOSE_ROOTFINDING_HPP