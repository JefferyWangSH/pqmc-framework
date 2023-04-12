/*
 *   pqmc_engine.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "pqmc_engine.h"
#include "model.h"
#include "pqmc_params.hpp"
#include "utils/numerical_stable.hpp"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Eigenvalues>

#include <iostream>

namespace PQMC {

    using MatrixStack = Eigen::MatrixXd;

    void PqmcEngine::initial( const PqmcParams& params, const Model::Hubbard& model )
    {
        this->m_theta = params.theta;
        this->m_dt = params.dt;
        this->m_nt = params.nt;
        this->m_nl = params.nl;
        this->m_ns = params.nl * params.nl;
        this->m_np = params.np;
        this->m_stabilization_pace = params.stabilization_pace;
        
        Eigen::SelfAdjointEigenSolver<ProjectionMat> eigensolver( model.m_K );
        if ( eigensolver.info() != Eigen::Success ) {
            std::cerr << "PQMC::PqmcEngine::initial(): diagonalization of the hopping matrix failed, check the input." << std::endl;
            exit(1);
        }
        
        // for SU(2) fermion, we seperate the spin-up and down sector.
        // the trial state is choson as the many-body ground state in the free theory.        
        assert( this->m_np % 2 == 0 );
        // if Ns <= Np, the green's function completely vanishes due to the pauli principle
        assert( 2*this->m_ns > this->m_np );
        this->m_projection_mat.resize( this->m_ns, this->m_np/2 );
        this->m_projection_mat.setZero();
        for ( auto i = 0; i < this->m_np/2; ++i ){
            this->m_projection_mat.col(i) = eigensolver.eigenvectors().col(i);
        }

        // initialize SVD stacks and green's function
        this->m_svd_stack_left_up  = new Utils::SvdStackReal( {this->m_ns, this->m_np/2}, this->m_nt, this->m_projection_mat );
        this->m_svd_stack_left_dn  = new Utils::SvdStackReal( {this->m_ns, this->m_np/2}, this->m_nt, this->m_projection_mat );
        this->m_svd_stack_right_up = new Utils::SvdStackReal( {this->m_ns, this->m_np/2}, this->m_nt, this->m_projection_mat );
        this->m_svd_stack_right_dn = new Utils::SvdStackReal( {this->m_ns, this->m_np/2}, this->m_nt, this->m_projection_mat );
        this->m_green_tt_up = new GreensFunction( this->m_ns, this->m_ns );
        this->m_green_tt_dn = new GreensFunction( this->m_ns, this->m_ns );

        MatrixStack temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
        MatrixStack temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
        
        // initial SVD stacks for sweeping usages
        for ( auto t = 0; t < this->m_nt; ++t ) {
            model.multiply_B_from_left( temp_stack_up, t, +1 );
            model.multiply_B_from_left( temp_stack_dn, t, -1 );

            // stabilize every `stabilization_pace` steps using SVD decomposition
            if ( (t + 1) % this->m_stabilization_pace == 0 ) {
                this->m_svd_stack_right_up->push( temp_stack_up );
                this->m_svd_stack_right_dn->push( temp_stack_dn );
                temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
                temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
            }
        }

        // initialize green's function
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up );
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn );
    }


} // namespace PQMC