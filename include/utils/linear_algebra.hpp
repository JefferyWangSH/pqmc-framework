#ifndef UTILS_LINEAR_ALGEBRA_HPP
#define UTILS_LINEAR_ALGEBRA_HPP
#pragma once

/*
 *   This source file includes some diagonalizing tools with C++/Eigen interface
 *   for diagonalizing (long) real and (long) complex matrices using mkl-implemented lapack.
 *   including:
 *     1. generalized SVD decomposition for arbitrary M * N matrices ( dgesvd and zgesvd )
 *     2. optimized diagonalizing mechanism for N * N real symmetric matrix ( dsyev and zsyev )
 *   the numerical accuracy and efficiency of the calculations are guaranteed.
 */


// todo: the simple redefinition seems not to work
// #define MKL_Complex16 std::complex<double>

#include <iostream>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <mkl_lapacke.h>


namespace Utils {

    // -------------------------------------  Utils::LinearAlgebra class  ----------------------------------------
    class LinearAlgebra {

        public:


            /**
             *  SVD decomposition of an arbitrary M * N real matrix, using LAPACKE_dgesvd:
             *       A  ->  U * S * V^T
             *  note that V is returned in this subroutine, not V transposed.
             *
             *  @param row -> number of rows.
             *  @param col -> number of cols.
             *  @param mat -> an arbitrary `row` * `col` real matrix to be solved.
             *  @param u -> `row` * `row` orthogonal u matrix of type Eigen::MatrixXd.
             *  @param s -> eigenvalues s of type Eigen::VectorXd, descending sorted.
             *  @param v -> `col` * `col` orthogonal v matrix of type Eigen::MatrixXd.
             */
            static void mkl_lapack_dgesvd( const int& row, 
                                           const int& col, 
                                           const Eigen::MatrixXd& mat, 
                                           Eigen::MatrixXd& u, 
                                           Eigen::VectorXd& s, 
                                           Eigen::MatrixXd& v  ) 
            {
                assert( row == mat.rows() );
                assert( col == mat.cols() );
                // todo: currently, the subroutine fails if the input matrix has different rows and columns
                assert( row == col );

                // info of the input matrix
                MKL_INT matrix_layout = LAPACK_ROW_MAJOR;
                MKL_INT info, lda = row, ldu = row, ldvt = col;

                // local arrays
                double tmp_s[ldu * ldu], tmp_u[ldu * row], tmp_vt[ldvt * col];
                double tmp_a[lda * col];
                double tmp_superb[ldu * lda];

                // convert the eigen matrix to c-style array
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(&tmp_a[0], lda, col) = mat;

                // compute SVD
                info = LAPACKE_dgesvd( matrix_layout, 'A', 'A', row, col, tmp_a, lda, tmp_s, tmp_u, ldu, tmp_vt, ldvt, tmp_superb );

                // check for convergence
                if( info > 0 ) {
                    std::cerr << "Utils::LinearAlgebra::mkl_lapack_dgesvd(): "
                              << "the algorithm computing SVD failed to converge." << std::endl;
                    exit(1);
                }

                // convert the results into the Eigen style
                u = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tmp_u, col, col);
                s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(tmp_s, 1, col);
                v = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(tmp_vt, row, row);
            }


            /**
             *  SVD decomposition of an arbitrary M * N complex matrix, using LAPACKE_zgesvd:
             *       A  ->  U * S * V^H
             *  note that V is returned in this subroutine, not V conjugate transposed.
             *
             *  @param row -> number of rows.
             *  @param col -> number of cols.
             *  @param mat -> an arbitrary `row` * `col` complex matrix to be solved.
             *  @param u -> `row` * `row` unitary u matrix of type Eigen::MatrixXcd, 
             *  @param s -> eigenvalues s of type Eigen::VectorXd, which are real, non-negative and descending sorted.
             *  @param v -> `col` * `col` unitary v matrix of type Eigen::MatrixXcd
             */
            static void mkl_lapack_zgesvd( const int& row, 
                                           const int& col, 
                                           const Eigen::MatrixXcd& mat, 
                                           Eigen::MatrixXcd& u, 
                                           Eigen::VectorXd& s, 
                                           Eigen::MatrixXcd& v  )
            {
                assert( row == mat.rows() );
                assert( col == mat.cols() );
                assert( row == col );

                // info of the input matrix
                MKL_INT matrix_layout = LAPACK_ROW_MAJOR;
                MKL_INT info, lda = col, ldu = row, ldvh = col;

                // local arrays
                double tmp_s[row];
                double tmp_superb[std::min(row,col)-1];
                MKL_Complex16 tmp_u[ldu * row], tmp_vh[ldvh * col];
                MKL_Complex16 tmp_a[lda * row];

                // convert the eigen matrix to c-style array
                for ( auto r = 0; r < row; ++r ){
                    for ( auto c = 0; c < lda; ++c ) {
                        tmp_a[row*r+c].real = mat.real()(r, c);
                        tmp_a[row*r+c].imag = mat.imag()(r, c);
                    }
                }

                // compute SVD
                info = LAPACKE_zgesvd( matrix_layout, 'A', 'A', row, col, tmp_a, lda, tmp_s, tmp_u, ldu, tmp_vh, ldvh, tmp_superb );

                // check for convergence
                if( info > 0 ) {
                    std::cerr << "Utils::LinearAlgebra::mkl_lapack_zgesvd(): "
                              << "the algorithm computing SVD failed to converge." << std::endl;
                    exit(1);
                }

                // convert the results into the Eigen style
                u.resize(row, col);
                v.resize(row, col);
                for ( auto r = 0; r < row; ++r ) {
                    for ( auto c = 0; c < col; ++c ) {
                        u.real()(r, c) =   tmp_u[row*r+c].real;
                        u.imag()(r, c) =   tmp_u[row*r+c].imag;
                        v.real()(r, c) =   tmp_vh[r+ldvh*c].real;
                        v.imag()(r, c) = - tmp_vh[r+ldvh*c].imag;
                    }
                }
                s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(tmp_s, 1, col);
            }



            /**
             *  Given an arbitrary N * N real symmetric matrix, calculate the eigenvalues and eigenstates using LAPACKE_dsyev
             *        A  ->  T^dagger * S * T
             *  where T is the rotation matrix, which is orthogonal;
             *  and S is a diagonal matrix with the eigenvalues being diagonal elements.
             *
             *  @param size -> number of rows/cols.
             *  @param mat -> an arbitrary `size` * `size` real symmetric matrix to be solved.
             *  @param s -> diagonal eigen matrix.
             *  @param t -> rotation matrix, whose columns are the corresponding eigenstates.
             */
            static void mkl_lapack_dsyev( const int& size, 
                                          const Eigen::MatrixXd& mat, 
                                          Eigen::VectorXd& s, 
                                          Eigen::MatrixXd& t  ) 
            {
                assert( mat.rows() == size );
                assert( mat.cols() == size );
                // make sure the input matrix is symmetric
                assert( mat.isApprox(mat.transpose(), 1e-12) );

                // locals params
                MKL_INT matrix_layout = LAPACK_ROW_MAJOR;
                MKL_INT n = size, lda = size, info;
                double tmp_s[n];
                double tmp_a[lda * n];

                // convert the eigen matrix to c-style array
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(&tmp_a[0], lda, n) = mat;

                // solve eigen problem
                info = LAPACKE_dsyev( matrix_layout, 'V', 'U', n, tmp_a, lda, tmp_s );

                // check for convergence
                if( info > 0 ) {
                    std::cerr << "Utils::LinearAlgebra::mkl_lapack_dsyev(): " 
                              << "the algorithm failed to compute eigenvalues." << std::endl;
                    exit(1);
                }

                // convert the results into the Eigen style
                s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(tmp_s, 1, n);
                t = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(tmp_a, n, n);
            }



            /** 
             *   
             * 
             *  @param size 
             *  @param mat
             *  @param s
             *  @param t
             */
            static void mkl_lapack_zsyev( const int& size, 
                                          const Eigen::MatrixXd& mat, 
                                          Eigen::VectorXd& s, 
                                          Eigen::MatrixXd& t  )
            {

                // todo

            }


    };


} // namespace Utils


#endif // UTILS_LINEAR_ALGEBRA_HPP
