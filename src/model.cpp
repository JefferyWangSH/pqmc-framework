/*
 *   model.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "model.h"
#include "pqmc_params.hpp"
#include "pqmc_engine.h"
#include "random.h"
#include <unsupported/Eigen/MatrixFunctions>

namespace Model {

    using RowVector = Eigen::Matrix<ScalarType, 1, Eigen::Dynamic>;

    void Hubbard::initial( const PQMC::PqmcParams& params )
    {
        this->m_nt = params.nt;
        this->m_nl = params.nl;
        this->m_ns = params.nl * params.nl;
        this->m_dt = params.dt;

        this->m_t = params.t;
        this->m_u = params.u;
        this->m_alpha = std::acosh( std::exp(0.5 * this->m_dt * this->m_u) );
        // open boundary in imaginary-time direction
        //   t = dt, ..., dt*nt(=2theta).
        // nt time slices in total.
        this->m_ising_fields.resize( this->m_nt, this->m_ns );

        this->m_K.resize( this->m_ns, this->m_ns );
        this->m_K.setZero();
        for ( int x = 0; x < this->m_nl; ++x ) {
            for ( int y = 0; y < this->m_nl; ++y ) {
                const int i = x + y * this->m_nl;
                const int ixplus1 = (x+1)%this->m_nl + y * this->m_nl;
                const int iyplus1 = x + (y+1)%this->m_nl * this->m_nl;
                this->m_K( i, ixplus1 ) += - this->m_t;
                this->m_K( ixplus1, i ) += - this->m_t;
                this->m_K( i, iyplus1 ) += - this->m_t;
                this->m_K( iyplus1, i ) += - this->m_t;
            }
        }
        this->m_expK     = ( -this->m_dt * this->m_K ).exp();
        this->m_inv_expK = ( +this->m_dt * this->m_K ).exp();
    }

    void Hubbard::randomly_initial_ising_fields()
    {
        std::bernoulli_distribution bernoulli_dist(0.5);
        for ( auto t = 0; t < this->m_nt; ++t ) {
            for ( auto i = 0; i < this->m_ns; ++i ) {
                this->m_ising_fields(t, i) = bernoulli_dist( Utils::Random::Engine )? +1.0 : -1.0;
            }
        }
    }

    void Hubbard::update_ising_field( timeIndex t, spaceIndex i )
    {
        assert( t >= 0 && t < this->m_nt );
        assert( i >= 0 && i < this->m_ns );
        this->m_ising_fields(t, i) = - this->m_ising_fields(t, i);
    }

    void Hubbard::update_greens_function( PQMC::PqmcEngine& engine, timeIndex t, spaceIndex i )
    {
        assert( t >= 0 && t < this->m_nt );
        assert( i >= 0 && i < this->m_ns );
        GreensFunction& green_tt_up = *engine.m_green_tt_up;
        GreensFunction& green_tt_dn = *engine.m_green_tt_dn;

        this->m_delta_up = std::exp( -2*this->m_alpha*this->m_ising_fields(t,i) ) - 1;
        this->m_delta_dn = std::exp( +2*this->m_alpha*this->m_ising_fields(t,i) ) - 1;
        const double factor_up = this->m_delta_up / ( 1 + ( 1-green_tt_up(i,i) ) * this->m_delta_up );
        const double factor_dn = this->m_delta_dn / ( 1 + ( 1-green_tt_dn(i,i) ) * this->m_delta_dn );

        green_tt_up -= factor_up * green_tt_up.col(i) * ( RowVector::Unit(this->m_ns, i) - green_tt_up.row(i) );
        green_tt_dn -= factor_dn * green_tt_dn.col(i) * ( RowVector::Unit(this->m_ns, i) - green_tt_dn.row(i) );
    }

    const double Hubbard::get_accepting_ratio( const PQMC::PqmcEngine& engine, timeIndex t, spaceIndex i ) const
    {
        assert( t >= 0 && t < this->m_nt );
        assert( i >= 0 && i < this->m_ns );
        const GreensFunction& green_tt_up = *engine.m_green_tt_up;
        const GreensFunction& green_tt_dn = *engine.m_green_tt_dn;
        return  ( 1 + ( 1-green_tt_up(i,i) ) * ( std::exp( -2*this->m_alpha*this->m_ising_fields(t,i) ) - 1 ) )
              * ( 1 + ( 1-green_tt_dn(i,i) ) * ( std::exp( +2*this->m_alpha*this->m_ising_fields(t,i) ) - 1 ) );
    }

    void Hubbard::multiply_B_from_left( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        assert( t >= 0 && t < this->m_nt );
        green = this->m_expK * green;
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.row(i) *= std::exp( +s * this->m_alpha * this->m_ising_fields(t, i) );
        }
    }

    void Hubbard::multiply_B_from_right( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        assert( t >= 0 && t < this->m_nt );
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.col(i) *= std::exp( +s * this->m_alpha * this->m_ising_fields(t, i) );
        }
        green = green * this->m_expK;
    }

    void Hubbard::multiply_invB_from_left( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        assert( t >= 0 && t < this->m_nt );
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.row(i) *= std::exp( -s * this->m_alpha * this->m_ising_fields(t, i) );
        }
        green = this->m_inv_expK * green;
    }

    void Hubbard::multiply_invB_from_right( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        assert( t >= 0 && t < this->m_nt );
        green = green * this->m_inv_expK;
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.col(i) *= std::exp( -s * this->m_alpha * this->m_ising_fields(t, i) );
        }
    }

    void Hubbard::multiply_transB_from_left( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        assert( t >= 0 && t < this->m_nt );
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.row(i) *= std::exp( +s * this->m_alpha * this->m_ising_fields(t, i) );
        }
        // since K is symmetric, the exponential of K is also symmetric.
        green = this->m_expK * green;
    }

} // namespace Model