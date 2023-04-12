/*
 *   svd_stack.h
 * 
 *     Created on: Apr 12, 2023
 *         Author: Jeffery Wang
 * 
 *   This head file includes SvdData and SvdStack class for 
 *   the stable multiplication of long chains of dense matrices.
 *   LAPACK libraries (MKL implemented) is needed for the SVD decomposition.
 */

#ifndef UTILS_SVD_STACK_H
#define UTILS_SVD_STACK_H

#include <complex>
#include <iostream>
#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

#include "utils/linear_algebra.hpp"


namespace Utils {
    

    // ----------------------------------------------  Utils::SvdData< ScalarType >  ------------------------------------------------
    template< typename ScalarType > class SvdData {
        
        private:

            using uMatrix = Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic >;
            using vMatrix = Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic >;
            using sVector = Eigen::VectorXd;

            uMatrix m_u{};
            vMatrix m_v{};
            sVector m_s{};
        
        public:

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            SvdData() = default;

            uMatrix& MatrixU()        { return this->m_u; }
            vMatrix& MatrixV()        { return this->m_v; }
            sVector& SingularValues() { return this->m_s; }

    };



    // ---------------------------------------------  Utils::SvdStack< ScalarType >  ------------------------------------------------
    
    // udv stack for a chain of matrix products: U * D * V^H = An * ... * A2 * A1 * A0

    template< typename ScalarType > class SvdStack {
        
        private:
            
            using SvdDataVector = std::vector< SvdData<ScalarType> >;
            using Matrix = Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic >;
            using Vector = Eigen::VectorXd;
            using MatShape = std::array<int,2>;

            SvdDataVector m_stack{};             // matrix stacks
            MatShape      m_matrix_shape{};      // shape of the matrices in the stack
            int           m_stack_length{0};     // curreng length of the stack
            Matrix        m_temp_matrix{};       // intermediate variable for computing the svd
            Matrix        m_initial_matrix{};    // initial matrix A0, with shape of `m_matrix_shape`


        public:

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            // ---------------------------------------------  Constructions  ----------------------------------------------------

            SvdStack() = default;

            explicit SvdStack( const MatShape& shape, int max_stack_length, const Matrix& initial_matrix )
                : m_matrix_shape(shape), m_initial_matrix(initial_matrix)
            {
                assert( initial_matrix.rows() == shape[0] && initial_matrix.cols() == shape[1] );
                this->m_stack.reserve( max_stack_length );
                for ( auto i = 0; i < max_stack_length; ++i ) {
                    this->m_stack.emplace_back();
                }
            }


            // -----------------------------------------------  Interfaces  -----------------------------------------------------
            
            bool empty() const { return this->m_stack_length == 0; }
            int  CurrentStackLength() const { return this->m_stack_length; }
            const MatShape MatrixShape() const { return this->m_matrix_shape; }

            // return the SVD decomposition matrices for the current stack
            const Vector SingularValues()
            {
                assert( this->m_stack_length > 0 );
                return this->m_stack[this->m_stack_length-1].SingularValues();
            }

            const Matrix MatrixU()
            {
                assert( this->m_stack_length > 0 );
                return this->m_stack[this->m_stack_length-1].MatrixU();
            }

            const Matrix MatrixV()
            {
                assert( this->m_stack_length > 0 );
                Matrix r = this->m_stack[0].MatrixV();
                for ( auto i = 1; i < this->m_stack_length; ++i ) {
                    r = r * this->m_stack[i].MatrixV();
                }
                return r;
            }
            

            // ----------------------------------  Stack operations push(), pop() and clear()  ----------------------------------

            // multply a matrix to the stack from the left
            void push( const Matrix& matrix )
            {
                assert( matrix.rows() == this->m_matrix_shape[0] && matrix.cols() == this->m_matrix_shape[0] );
                assert( this->m_stack_length < (int)this->m_stack.size() );

                if ( this->m_stack_length == 0 ) {

                    if constexpr ( std::is_same_v< ScalarType, double > ) {
                        // SVD decomposition of real matrix
                        Utils::LinearAlgebra::mkl_lapack_dgesvd
                        (
                            this->m_matrix_shape[0],
                            this->m_matrix_shape[1],
                            matrix * this->m_initial_matrix,
                            this->m_stack[this->m_stack_length].MatrixU(),
                            this->m_stack[this->m_stack_length].SingularValues(),
                            this->m_stack[this->m_stack_length].MatrixV()
                        );
                    }

                    else if constexpr ( std::is_same_v< ScalarType, std::complex<double> > ) {
                        // SVD decomposition of complex matrix
                        Utils::LinearAlgebra::mkl_lapack_zgesvd
                        (
                            this->m_matrix_shape[0],
                            this->m_matrix_shape[1],
                            matrix * this->m_initial_matrix,
                            this->m_stack[this->m_stack_length].MatrixU(),
                            this->m_stack[this->m_stack_length].SingularValues(),
                            this->m_stack[this->m_stack_length].MatrixV()
                        );
                    }

                    else {
                        std::cerr << "Utils::SvdStack<ScalarType>::push(): "
                                  << "undefined scalar type."
                                  << std::endl;
                        exit(1);   
                    }

                }
                else {
                    // important: mind the order of multiplication.
                    // avoid mixing of different numerical scales here
                    this->m_temp_matrix = ( matrix * this->MatrixU() ) * this->SingularValues().asDiagonal();

                    if constexpr ( std::is_same_v< ScalarType, double > ) {
                        // SVD decomposition of real matrix
                        Utils::LinearAlgebra::mkl_lapack_dgesvd
                        (
                            this->m_temp_matrix.rows(),
                            this->m_temp_matrix.cols(),
                            this->m_temp_matrix,
                            this->m_stack[this->m_stack_length].MatrixU(),
                            this->m_stack[this->m_stack_length].SingularValues(),
                            this->m_stack[this->m_stack_length].MatrixV()
                        );
                    }

                    else if constexpr ( std::is_same_v< ScalarType, std::complex<double> > ) {
                        // SVD decomposition of complex matrix
                        Utils::LinearAlgebra::mkl_lapack_zgesvd
                        (
                            this->m_temp_matrix.rows(),
                            this->m_temp_matrix.cols(),
                            this->m_temp_matrix,
                            this->m_stack[this->m_stack_length].MatrixU(),
                            this->m_stack[this->m_stack_length].SingularValues(),
                            this->m_stack[this->m_stack_length].MatrixV()
                        );
                    }

                    else {
                        std::cerr << "Utils::SvdStack<ScalarType>::push(): "
                                  << "undefined scalar type."
                                  << std::endl;
                        exit(1);   
                    }

                }
                this->m_stack_length += 1;
            }

            // remove the most left matrix from the stack
            void pop() 
            {
                // notice that the memory is not actually released
                assert( this->m_stack_length > 0 );
                this->m_stack_length -= 1;
            }

            // clear the stack
            // simply set stack_length = 0, note that the memory is not really deallocated.
            void clear() { this->m_stack_length = 0; }

    };

    // some aliases
    using SvdStackReal = SvdStack<double>;
    using SvdStackCpx  = SvdStack<std::complex<double>>;


} // namespace Utils


#endif // UTILS_SVD_STACK_H