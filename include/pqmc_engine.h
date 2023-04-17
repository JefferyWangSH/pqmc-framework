/*
 *   pqmc_engine.h
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#ifndef PQMC_ENGINE_H
#define PQMC_ENGINE_H
#pragma once

#include <memory>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "utils/svd_stack.h"

namespace Model { class Hubbard; }
namespace Measure { class MeasureHandler; }

namespace PQMC {

    struct PqmcParams;

    using ScalarType = double;
    using timeIndex  = int;
    using GreensFunction       = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using VecGreensFunction    = std::vector<GreensFunction>;
    using ptrGreensFunction    = std::unique_ptr<GreensFunction>;
    using ptrVecGreensFunction = std::unique_ptr<VecGreensFunction>;
    using ProjectionMat        = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using SvdStack             = Utils::SvdStack<ScalarType>;
    using ptrSvdStack          = std::unique_ptr<SvdStack>;


    // ------------------------------------------------------------  PQMC::PqmcEngine  ---------------------------------------------------------------------------
    class PqmcEngine {

        public:

            double m_theta{};                    // half of the projection length, 2theta = dt * nt 
            double m_beta{};                     // projection window within which measurments are performed ( theta-beta >> 1 for efficient projection )
            double m_dt{};                       // imaginary-time spacing
            int    m_nt{};                       // imaginary-time slices
            int    m_ntm{};                      // imaginary-time slices in the measurement window
            int    m_nl{};                       // linear size of lattice
            int    m_ns{};                       // number of lattice sites
            int    m_np{};                       // particle number
            int    m_stabilization_pace{};       // pace of numerical stabilization
            int    m_current_time_slice{0};      // ranging from 0 to nt

            // equal-time green's function
            ptrGreensFunction m_green_tt_up{};
            ptrGreensFunction m_green_tt_dn{};
            ptrVecGreensFunction m_vec_green_tt_up{};
            ptrVecGreensFunction m_vec_green_tt_dn{};

            // time-displaced green's functions G(t,0) and G(0,t)
            ptrGreensFunction m_green_t0_up{};
            ptrGreensFunction m_green_t0_dn{};
            ptrGreensFunction m_green_0t_up{};
            ptrGreensFunction m_green_0t_dn{};
            ptrVecGreensFunction m_vec_green_t0_up{};
            ptrVecGreensFunction m_vec_green_t0_dn{};
            ptrVecGreensFunction m_vec_green_0t_up{};
            ptrVecGreensFunction m_vec_green_0t_dn{};

            ProjectionMat m_projection_mat{};

            ptrSvdStack m_svd_stack_left_up{};
            ptrSvdStack m_svd_stack_left_dn{};
            ptrSvdStack m_svd_stack_right_up{};
            ptrSvdStack m_svd_stack_right_dn{};

            ptrSvdStack m_svd_stack_dynamic_up{};
            ptrSvdStack m_svd_stack_dynamic_dn{};

            // bool parameters
            bool m_is_thermalization{};         // whether the MC simulation is in the thermalization phase
            bool m_is_equaltime{};              // whether to perform the equal-time measurements
            bool m_is_dynamic{};                // whether to perform the dynamic measurements

            double m_equaltime_wrap_error{0.0};
            double m_dynamic_wrap_error{0.0};

        public:
        
            void initial( const PqmcParams& params, const Model::Hubbard& model, const Measure::MeasureHandler& handler );

            void set_thermalization( bool is_thermalization );

            void sweep_from_0_to_2theta( Model::Hubbard& model );
            void sweep_from_2theta_to_0( Model::Hubbard& model );
            void sweep_for_dynamic_greens_function( Model::Hubbard& model );

        private:

            void allocate();
            
            void metropolis_update( Model::Hubbard& model, timeIndex t );

            void wrap_from_0_to_2theta( Model::Hubbard& model, timeIndex t );
            void wrap_from_2theta_to_0( Model::Hubbard& model, timeIndex t );

            bool isInMeasurementWindow( timeIndex t ) const;
            timeIndex Convert2WindowIndex( timeIndex t ) const;

    };

} // namespace PQMC

#endif // PQMC_ENGINE_H