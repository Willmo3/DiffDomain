//
// Created by will on 7/3/26.
//
#include "regression.hpp"

Eigen::VectorXd regress(Eigen::MatrixXd& xmat, const Eigen::VectorXd& yvec) {
    Eigen::MatrixXd A(xmat.rows(), xmat.cols() + 1);
    A.col(0) = Eigen::VectorXd::Ones(xmat.rows());
    A.block(0, 1, xmat.rows(), xmat.cols()) = xmat;
    return A.colPivHouseholderQr().solve(yvec);
}

Eigen::MatrixXd cartesian_product(const Eigen::VectorXd& xvec, const Eigen::VectorXd& yvec) {
    Eigen::MatrixXd result(xvec.size() * yvec.size(), 2);
    auto index = 0;

    for (auto i : xvec) {
        for (auto j : yvec) {
            result(index, 0) = i;
            result(index, 1) = j;
            index++;
        }
    }
    return result;
}

Eigen::MatrixXd cartesian_product(const Eigen::VectorXd &xvec, const Eigen::VectorXd &yvec, const Eigen::VectorXd &zvec) {
    Eigen::MatrixXd result(xvec.size() * yvec.size() * zvec.size(), 3);
    auto index = 0;

    for (auto i : xvec) {
        for (auto j : yvec) {
            for (auto k : zvec) {
                result(index, 0) = i;
                result(index, 1) = j;
                result(index, 2) = k;
                index++;
            }
        }
    }
    return result;
}

Eigen::MatrixXd cartesian_product(const Eigen::VectorXd &xvec, const Eigen::VectorXd &yvec, const Eigen::VectorXd &zvec, const Eigen::VectorXd &wvec) {
    Eigen::MatrixXd result(xvec.size() * yvec.size() * zvec.size() * wvec.size(), 4);
    auto index = 0;

    for (auto i : xvec) {
        for (auto j : yvec) {
            for (auto k : zvec) {
                for (auto l : wvec) {
                    result(index, 0) = i;
                    result(index, 1) = j;
                    result(index, 2) = k;
                    result(index, 3) = l;
                    index++;
                }
            }
        }
    }
    return result;
}

Eigen::MatrixXd get_corners(const Winterval &bounds1, const Winterval &bounds2, const Winterval &bounds3, const Winterval &bounds4) {
    auto x1 = Eigen::VectorXd(2);
    auto x2 = Eigen::VectorXd(2);
    auto x3 = Eigen::VectorXd(2);
    auto x4 = Eigen::VectorXd(2);
    x1 << bounds1.min(), bounds1.max();
    x2 << bounds2.min(), bounds2.max();
    x3 << bounds3.min(), bounds3.max();
    x4 << bounds4.min(), bounds4.max();
    return cartesian_product(x1, x2, x3, x4);
}
