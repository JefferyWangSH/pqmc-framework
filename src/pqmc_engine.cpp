/*
 *   pqmc_engine.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "pqmc_engine.h"
#include "model.h"
#include "measure_handler.h"
#include "random.h"
#include "pqmc_params.hpp"
#include "utils/numerical_stable.hpp"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Eigenvalues>

#include <iostream>

namespace PQMC {

    using MatrixStack = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

    void PqmcEngine::set_thermalization( bool is_thermalization ) { this->m_is_thermalization = is_thermalization; }
    
    bool PqmcEngine::isInMeasurementWindow( timeIndex t ) const
    {
        assert( t >= 0 && t < this->m_nt );
        return ( ( t > (this->m_nt-this->m_ntm)/2-1 ) && ( t <= (this->m_nt+this->m_ntm)/2-1 ) );
    }

    timeIndex PqmcEngine::Convert2WindowIndex( timeIndex t ) const
    {
        assert( t >= 0 && t < this->m_nt );
        return t - (int)(this->m_nt-this->m_ntm)/2;
    }


    void PqmcEngine::allocate()
    {
        // release the pointers if assigned before
        if ( this->m_svd_stack_left_up ) { this->m_svd_stack_left_up.reset(); }
        if ( this->m_svd_stack_left_dn ) { this->m_svd_stack_left_dn.reset(); }
        if ( this->m_svd_stack_right_up ) { this->m_svd_stack_right_up.reset(); }
        if ( this->m_svd_stack_right_dn ) { this->m_svd_stack_right_dn.reset(); }
        if ( this->m_green_tt_up ) { this->m_green_tt_up.reset(); }
        if ( this->m_green_tt_dn ) { this->m_green_tt_dn.reset(); }
        if ( this->m_green_t0_up ) { this->m_green_t0_up.reset(); }
        if ( this->m_green_t0_dn ) { this->m_green_t0_dn.reset(); }
        if ( this->m_green_0t_up ) { this->m_green_0t_up.reset(); }
        if ( this->m_green_0t_dn ) { this->m_green_0t_dn.reset(); }
        if ( this->m_vec_green_tt_up ) { this->m_vec_green_tt_up.reset(); }
        if ( this->m_vec_green_tt_dn ) { this->m_vec_green_tt_dn.reset(); }
        if ( this->m_vec_green_t0_up ) { this->m_vec_green_t0_up.reset(); }
        if ( this->m_vec_green_t0_dn ) { this->m_vec_green_t0_dn.reset(); }
        if ( this->m_vec_green_0t_up ) { this->m_vec_green_0t_up.reset(); }
        if ( this->m_vec_green_0t_dn ) { this->m_vec_green_0t_dn.reset(); }
        
        // allocate memory for SvdStack classes
        this->m_svd_stack_left_up  = std::make_unique<SvdStack>( std::array<int,2>({this->m_ns, this->m_np/2}), this->m_nt, this->m_projection_mat );
        this->m_svd_stack_left_dn  = std::make_unique<SvdStack>( std::array<int,2>({this->m_ns, this->m_np/2}), this->m_nt, this->m_projection_mat );
        this->m_svd_stack_right_up = std::make_unique<SvdStack>( std::array<int,2>({this->m_ns, this->m_np/2}), this->m_nt, this->m_projection_mat );
        this->m_svd_stack_right_dn = std::make_unique<SvdStack>( std::array<int,2>({this->m_ns, this->m_np/2}), this->m_nt, this->m_projection_mat );

        // allocate memories for green's functions
        this->m_green_tt_up = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
        this->m_green_tt_dn = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
        if ( this->m_is_equaltime || this->m_is_dynamic ) {
            this->m_vec_green_tt_up = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
            this->m_vec_green_tt_dn = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
        }

        if ( this->m_is_dynamic ) {
            this->m_green_t0_up = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
            this->m_green_t0_dn = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
            this->m_green_0t_up = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
            this->m_green_0t_dn = std::make_unique<GreensFunction>( this->m_ns, this->m_ns );
            this->m_vec_green_t0_up = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
            this->m_vec_green_t0_dn = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
            this->m_vec_green_0t_up = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
            this->m_vec_green_0t_dn = std::make_unique<VecGreensFunction>( this->m_ntm, GreensFunction(this->m_ns, this->m_ns) );
        }
    }


    void PqmcEngine::initial( const PqmcParams& params, const Model::Hubbard& model, const Measure::MeasureHandler& handler )
    {
        this->m_theta = params.theta;
        this->m_beta = params.beta;
        this->m_dt = params.dt;
        this->m_nt = params.nt;
        this->m_ntm = params.ntm;
        this->m_nl = params.nl;
        this->m_ns = params.nl * params.nl;
        this->m_np = params.np;
        this->m_stabilization_pace = params.stabilization_pace;

        this->m_is_thermalization = ( params.sweeps_warmup != 0 );
        this->m_is_equaltime = handler.isEqualTime();
        this->m_is_dynamic = handler.isDynamic();
        
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
        this->allocate();
        MatrixStack temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
        MatrixStack temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
        
        // initial SVD stacks for sweeping usages
        // the sweep will start from 0 up to 2theta, hence we initialize L stack here  
        for ( auto t = this->m_nt-1; t >= 0; --t ) {
            model.multiply_transB_from_left( temp_stack_up, t, +1 );
            model.multiply_transB_from_left( temp_stack_dn, t, -1 );

            // stabilize every `stabilization_pace` steps using SVD decomposition
            if ( t % this->m_stabilization_pace == 0 ) {
                this->m_svd_stack_left_up->push( temp_stack_up );
                this->m_svd_stack_left_dn->push( temp_stack_dn );
                temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
                temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
            }
        }

        // initialize green's function
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up );
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn );
    }


    void PqmcEngine::metropolis_update( Model::Hubbard& model, timeIndex t )
    {
        assert( t >= 0 && t < this->m_nt );
        assert( this->m_current_time_slice == t );
        
        for ( auto i = 0; i < this->m_ns; ++i ) {
            // obtain the ratio of flipping the ising field at (t,i)
            const auto accepting_ratio = model.get_accepting_ratio( *this, t, i );

            if ( std::bernoulli_distribution( std::min(1.0, std::abs(accepting_ratio)) )( Utils::Random::Engine ) )
            {   
                // if accepted, update the green's function
                model.update_greens_function( *this, t, i );

                // update the ising field at (t,i) after updating the green's function
                model.update_ising_field( t, i );
            }
        }
    }


    void PqmcEngine::wrap_from_0_to_2theta( Model::Hubbard& model, timeIndex t )
    {
        assert( t >= 0 && t < this->m_nt );
        model.multiply_B_from_left     ( *this->m_green_tt_up, t, +1 );
        model.multiply_invB_from_right ( *this->m_green_tt_up, t, +1 );
        model.multiply_B_from_left     ( *this->m_green_tt_dn, t, -1 );
        model.multiply_invB_from_right ( *this->m_green_tt_dn, t, -1 );
    }


    void PqmcEngine::wrap_from_2theta_to_0( Model::Hubbard& model, timeIndex t )
    {
        assert( t >= 0 && t < this->m_nt );
        model.multiply_B_from_right   ( *this->m_green_tt_up, t, +1 );
        model.multiply_invB_from_left ( *this->m_green_tt_up, t, +1 );
        model.multiply_B_from_right   ( *this->m_green_tt_dn, t, -1 );
        model.multiply_invB_from_left ( *this->m_green_tt_dn, t, -1 );
    }


    void PqmcEngine::sweep_from_0_to_2theta( Model::Hubbard& model )
    {
        assert( this->m_current_time_slice == 0 );
        assert( this->m_svd_stack_right_up->empty() && this->m_svd_stack_right_dn->empty() );

        MatrixStack temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
        MatrixStack temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );

        for ( auto t = 0; t < this->m_nt; ++t ) {
            
            this->wrap_from_0_to_2theta( model, t );

            this->metropolis_update( model, t );

            if ( !this->m_is_thermalization && this->m_is_equaltime && this->isInMeasurementWindow(t) ) {
                (*this->m_vec_green_tt_up)[this->Convert2WindowIndex(t)] = *this->m_green_tt_up;
                (*this->m_vec_green_tt_dn)[this->Convert2WindowIndex(t)] = *this->m_green_tt_dn;
            }

            model.multiply_B_from_left( temp_stack_up, t, +1 );
            model.multiply_B_from_left( temp_stack_dn, t, -1 );

            // perform stabilization
            if ( (t+1) % this->m_stabilization_pace == 0 || t == this->m_nt-1 ) {
                
                // update SVD stacks
                this->m_svd_stack_left_up->pop();
                this->m_svd_stack_left_dn->pop();
                this->m_svd_stack_right_up->push( temp_stack_up );
                this->m_svd_stack_right_dn->push( temp_stack_dn );

                // collect the wrapping errors
                GreensFunction temp_green_tt_up( this->m_ns, this->m_ns );
                GreensFunction temp_green_tt_dn( this->m_ns, this->m_ns );
                double temp_wrap_error_tt_up = 0.0;
                double temp_wrap_error_tt_dn = 0.0;

                // recompute green's function
                Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, temp_green_tt_up );
                Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, temp_green_tt_dn );

                // compute wrapping errors
                Utils::NumericalStable::matrix_compare_error<ScalarType>(temp_green_tt_up, *this->m_green_tt_up, temp_wrap_error_tt_up);
                Utils::NumericalStable::matrix_compare_error<ScalarType>(temp_green_tt_dn, *this->m_green_tt_dn, temp_wrap_error_tt_dn);
                this->m_wrap_error = std::max( this->m_wrap_error, std::max(temp_wrap_error_tt_up, temp_wrap_error_tt_dn) );

                *this->m_green_tt_up = temp_green_tt_up;
                *this->m_green_tt_dn = temp_green_tt_dn;

                if ( !this->m_is_thermalization && this->m_is_equaltime && this->isInMeasurementWindow(t) ) {
                    (*this->m_vec_green_tt_up)[this->Convert2WindowIndex(t)] = *this->m_green_tt_up;
                    (*this->m_vec_green_tt_dn)[this->Convert2WindowIndex(t)] = *this->m_green_tt_dn;
                }

                temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
                temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
            }

            // finally stop at t = nt
            this->m_current_time_slice++;
        }

        // end with the fresh green's function
        // the measurement window basically does not reach the maximum projection length.
        if ( !this->m_is_thermalization && this->m_is_equaltime && this->isInMeasurementWindow(this->m_nt-1) ) {
            (*this->m_vec_green_tt_up)[this->Convert2WindowIndex(this->m_nt-1)] = *this->m_green_tt_up;
            (*this->m_vec_green_tt_dn)[this->Convert2WindowIndex(this->m_nt-1)] = *this->m_green_tt_dn;
        }
    }


    void PqmcEngine::sweep_from_2theta_to_0( Model::Hubbard& model )
    {
        assert( this->m_current_time_slice == this->m_nt );
        assert( this->m_svd_stack_left_up->empty() && this->m_svd_stack_left_dn->empty() );

        MatrixStack temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
        MatrixStack temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
        
        for ( auto t = this->m_nt-1; t >= 0; --t ) {

            // finally stop at t = 0
            this->m_current_time_slice--;

            // perform stabilization
            if ( (t+1) % this->m_stabilization_pace == 0 && t != this->m_nt-1 ) {

                // update SVD stacks
                this->m_svd_stack_right_up->pop();
                this->m_svd_stack_right_dn->pop();
                this->m_svd_stack_left_up->push( temp_stack_up );
                this->m_svd_stack_left_dn->push( temp_stack_dn );

                // collect the wrapping errors
                GreensFunction temp_green_tt_up( this->m_ns, this->m_ns );
                GreensFunction temp_green_tt_dn( this->m_ns, this->m_ns );
                double temp_wrap_error_tt_up = 0.0;
                double temp_wrap_error_tt_dn = 0.0;

                // recompute green's function
                Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, temp_green_tt_up );
                Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, temp_green_tt_dn );

                // compute wrapping errors
                Utils::NumericalStable::matrix_compare_error<ScalarType>(temp_green_tt_up, *this->m_green_tt_up, temp_wrap_error_tt_up);
                Utils::NumericalStable::matrix_compare_error<ScalarType>(temp_green_tt_dn, *this->m_green_tt_dn, temp_wrap_error_tt_dn);
                this->m_wrap_error = std::max( this->m_wrap_error, std::max(temp_wrap_error_tt_up, temp_wrap_error_tt_dn) );

                *this->m_green_tt_up = temp_green_tt_up;
                *this->m_green_tt_dn = temp_green_tt_dn;

                temp_stack_up = MatrixStack::Identity( this->m_ns, this->m_ns );
                temp_stack_dn = MatrixStack::Identity( this->m_ns, this->m_ns );
            }
            
            this->metropolis_update( model, t );

            if ( !this->m_is_thermalization && this->m_is_equaltime && this->isInMeasurementWindow(t) ) {
                (*this->m_vec_green_tt_up)[this->Convert2WindowIndex(t)] = *this->m_green_tt_up;
                (*this->m_vec_green_tt_dn)[this->Convert2WindowIndex(t)] = *this->m_green_tt_dn;
            }

            model.multiply_transB_from_left( temp_stack_up, t, +1 );
            model.multiply_transB_from_left( temp_stack_dn, t, -1 );

            this->wrap_from_2theta_to_0( model, t );
        }

        // at time slice t = 0
        this->m_svd_stack_right_up->pop();
        this->m_svd_stack_right_dn->pop();
        this->m_svd_stack_left_up->push( temp_stack_up );
        this->m_svd_stack_left_dn->push( temp_stack_dn );
        
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up );
        Utils::NumericalStable::compute_equaltime_greens<PqmcEngine>( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn );
        
        // end with the fresh green's function
        // the measurement window basically does not reach t=0 point ( no projection ).
        if ( !this->m_is_thermalization && this->m_is_equaltime && this->isInMeasurementWindow(0) ) {
            (*this->m_vec_green_tt_up)[this->Convert2WindowIndex(0)] = *this->m_green_tt_up;
            (*this->m_vec_green_tt_dn)[this->Convert2WindowIndex(0)] = *this->m_green_tt_dn;
        }
    }


    void PqmcEngine::sweep_for_dynamic_greens_function( const Model::Hubbard& model )
    {
        // todo
    }


} // namespace PQMC