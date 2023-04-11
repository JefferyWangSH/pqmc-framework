#ifndef MODEL_H
#define MODEL_H
#pragma once

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

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
            int m_nt{80};
            int m_nl{4};
            int m_ns{};
            double m_dt{0.05};
            
            double m_t{1.0};
            double m_u{4.0};
            double m_alpha{};
            
            HoppingMat m_K{};
            HoppingMat m_expK{};
            HoppingMat m_inv_expK{};

            SpaceTimeLattice m_ising_fields{};
        
            void initial();
            void set_ising_fields_to_random();

            void update_ising_field       ( timeIndex t, spaceIndex i );
            void update_greens_function   ( timeIndex t, spaceIndex i, spinIndex s );
            const double get_update_ratio ( timeIndex t, spaceIndex i, spinIndex s ) const ;

            void multiply_B_from_left       ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_B_from_right      ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_invB_from_left    ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
            void multiply_invB_from_right   ( GreensFunction& green, timeIndex t, spinIndex s ) const ;
    };

} // namespace Model

#endif // MODEL_H