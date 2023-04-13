/*
 *   numerical_stable.hpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

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

            template<typename ScalarType>
            using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

            template <class QmcEngine, typename ScalarType>
            static void compute_equaltime_greens( SvdStack<ScalarType>& left, SvdStack<ScalarType>& right, Matrix<ScalarType>& gtt )
            {
                if constexpr ( std::is_same_v< QmcEngine, PQMC::PqmcEngine > ) {
                    // define L = P^\dag B(2\theta,t) = VSU^\dag ;
                    //        R = B(t,0) P = USV^\dag .
                    // L, R are in general rectangular matrices

                    if constexpr ( std::is_same_v< ScalarType, double > ) {
                        gtt = Matrix<double>::Identity(gtt.rows(), gtt.cols()) - ( right.MatrixU() * ( left.MatrixU().transpose()*right.MatrixU() ).inverse() ) * left.MatrixU().transpose();
                    }
                    else if constexpr ( std::is_same_v< ScalarType, std::complex<double> > ) {
                        gtt = Matrix<std::complex<double>>::Identity(gtt.rows(), gtt.cols()) - ( right.MatrixU() * ( left.MatrixU().adjoint()*right.MatrixU() ).inverse() ) * left.MatrixU().adjoint();
                    }
                    else {
                        std::cerr << "Utils::NumericalStable::compute_equaltime_greens<QmcEngine>(): "
                                  << "undefined scalar type in Utils::SvdStack<ScalarType>."
                                  << std::endl;
                        exit(1);
                    }
                }
                else if constexpr ( std::is_same_v< QmcEngine, DQMC::DqmcEngine > ) {
                    // todo
                }
                else {
                    std::cerr << "Utils::NumericalStable::compute_equaltime_greens<QmcEngine>(): "
                              << "undefined QmcEngine type."
                              << std::endl;
                    exit(1);
                }
            }

    };


} // namespace Utils

#endif // UTILS_NUMERICAL_STABLE_HPP