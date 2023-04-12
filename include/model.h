/*
 *   model.h
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#ifndef MODEL_H
#define MODEL_H
#pragma once

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace PQMC { class PqmcEngine; struct PqmcParams; }

namespace Model{

    using HoppingMat       = Eigen::MatrixXd;
    using ProjectionMat    = Eigen::MatrixXd;
    using InteractionMat   = Eigen::MatrixXd;
    using GreensFunction   = Eigen::MatrixXd;
    using SpaceTimeLattice = Eigen::MatrixXd;

    using timeIndex  = int;
    using spaceIndex = int;
    using spinIndex  = int;
    
    class Hubbard{
        public:
            int m_nt{};
            int m_nl{};
            int m_ns{};
            double m_dt{};
            
            double m_t{};
            double m_u{};
            double m_alpha{};
            
            HoppingMat m_K{};
            HoppingMat m_expK{};
            HoppingMat m_inv_expK{};

            SpaceTimeLattice m_ising_fields{};
        
            void initial( const PQMC::PqmcParams& params );
            void randomly_initial_ising_fields();

            void update_ising_field         ( timeIndex t, spaceIndex i );
            void update_greens_function     ( PQMC::PqmcEngine& engine, timeIndex t, spaceIndex i );
            const double get_updating_ratio ( const PQMC::PqmcEngine& engine, timeIndex t, spaceIndex i, spinIndex s ) const ;

            void multiply_B_from_left       ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_B_from_right      ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_invB_from_left    ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_invB_from_right   ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
    };

} // namespace Model

#endif // MODEL_H