//
// Created by will on 2/16/26.
//

#ifndef DIFFDOMAIN_REGRESSION_HPP
#define DIFFDOMAIN_REGRESSION_HPP

#include "Eigen/Dense"
#include "Winterval/Winterval.hpp"

/**
 *
 * @param xvec n length vector of values
 * @param yvec n length vector of values
 * @return n^2 x 2 matrix, the cartesian product of these vectors. Each row represents an (x, y) pair.
 */
Eigen::MatrixXd cartesian_product(const Eigen::VectorXd &xvec, const Eigen::VectorXd &yvec);

/**
 *
 * @param xvec n length vector of values
 * @param yvec n length vector of values
 * @param zvec n length vector of values
 * @return n^3 x 3 matrix, the cartesian product of these vectors. Each row represents an (x, y, z) triple.
 */
Eigen::MatrixXd cartesian_product(const Eigen::VectorXd &xvec, const Eigen::VectorXd &yvec, const Eigen::VectorXd &zvec);

/**
 *
 * @param xvec n length vector of values
 * @param yvec n length vector of values
 * @param zvec n length vector of values
 * @param wvec n length vector of values
 * @return n^4 x 4 matrix, the cartesian product of these vectors. Each row represents an (x, y, z, w) tuple.
 */
Eigen::MatrixXd cartesian_product(const Eigen::VectorXd &xvec, const Eigen::VectorXd &yvec, const Eigen::VectorXd &zvec, const Eigen::VectorXd &wvec);
/**
 *
 * @param xmat Matrix of values to regress over. Each row is a different input dimension.
 * @param yvec Vector of output values to regress against. Each row is a different input dimension.
 * @return The coefficients for the regression, in the form of a vector. The first value is the intercept, and the rest are the coefficients for each input dimension, in order.
 */
Eigen::VectorXd regress(Eigen::MatrixXd &xmat, const Eigen::VectorXd &yvec);

/**
 * Take the cartesian product of the corners of the four intervals.
 * @return Size 16 matrix, representing every permutation of corner values.
 */
Eigen::MatrixXd get_corners(const Winterval &bounds1, const Winterval &bounds2, const Winterval &bounds3, const Winterval &bounds4);

#endif //DIFFDOMAIN_REGRESSION_HPP