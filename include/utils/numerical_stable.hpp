#ifndef UTILS_NUMERICAL_STABLE_HPP
#define UTILS_NUMERICAL_STABLE_HPP
#pragma once

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include "utils/svd_stack.h"

namespace PQMC { class PqmcEngine; }
namespace DQMC { class DqmcEngine; }

namespace Utils {

    class NumericalStable {
        public:
            using Matrix   = Eigen::MatrixXd;
            using SvdStack = SvdStackReal;

            template<class T>
            static void compute_equaltime_greens( SvdStack& left, SvdStack& right, Matrix &gtt ){};

    };

    template<>
    void NumericalStable::compute_equaltime_greens<PQMC::PqmcEngine>( SvdStack& left, SvdStack& right, Matrix &gtt ){
        // define L = P^\dag B(2\theta,t) = VSU^\dag ;
        //        R = B(t,0) P = USV^\dag .
        // L, R are in general rectangular matrices
        gtt = Matrix::Identity(gtt.rows(), gtt.cols()) - ( right.MatrixU() * ( left.MatrixU().transpose()*right.MatrixU() ).inverse() ) * left.MatrixU().transpose();
    }

    template<>
    void NumericalStable::compute_equaltime_greens<DQMC::DqmcEngine>( SvdStack& left, SvdStack& right, Matrix &gtt ){
        // todo
    }

} // namespace Utils

#endif // UTILS_NUMERICAL_STABLE_HPP