//
// Created by will on 6/9/26.
//

#include "rootfinding.hpp"
#include <Eigen/Eigenvalues>

std::vector<std::complex<double>> poly_roots(const std::vector<double>& coeffs) {
    auto n = coeffs.size() - 1; // degree

    // Build companion matrix
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n, n);
    auto leading = coeffs[0];
    for (int i = 0; i < n; i++) {
        C(i, n - 1) = -coeffs[n - i] / leading;
    }
    for (int i = 1; i < n; i++) {
        C(i, i - 1) = 1.0;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(C);
    auto eigs = solver.eigenvalues();

    std::vector<std::complex<double>> roots(n);
    for (int i = 0; i < n; ++i) {
        roots[i] = eigs[i];
    }
    return roots;
}