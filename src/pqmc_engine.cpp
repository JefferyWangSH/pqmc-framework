/*
 *   pqmc_engine.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "pqmc_engine.h"
#include "model.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Eigenvalues>

#include <iostream>

namespace PQMC {

    void PqmcEngine::initial( const Model::Hubbard& model )
    {
        Eigen::SelfAdjointEigenSolver<ProjectionMat> eigensolver( model.m_K );
        if ( eigensolver.info() != Eigen::Success ) {
            std::cerr << "PQMC::PqmcEngine::initial(): failure in diagonalizing the hopping matrix, check the input." << std::endl;
            exit(1);
        }
        // for su(2) fermion, we seperate the spin-up and down sector.
        // the trial state is choson as the many-body ground state in the free theory.        
        assert( this->m_particle_num % 2 == 0 );
        assert( this->m_ns >= this->m_particle_num/2 );
        this->m_projection_mat.resize( this->m_ns, this->m_particle_num/2 );
        this->m_projection_mat.setZero();
        for ( auto i = 0; i < this->m_particle_num/2; ++i ){
            this->m_projection_mat.col(i) = eigensolver.eigenvectors().col(i);
        }

        // this->m_svd_stack_left_up  = new Utils::SvdStackReal( this->m_ns, this->m_nt );
        // this->m_svd_stack_left_dn  = new Utils::SvdStackReal( this->m_ns, this->m_nt );
        // this->m_svd_stack_right_up = new Utils::SvdStackReal( this->m_ns, this->m_nt );
        // this->m_svd_stack_right_dn = new Utils::SvdStackReal( this->m_ns, this->m_nt );
    


    }


} // namespace PQMC