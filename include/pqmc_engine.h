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

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "utils/svd_stack.h"

namespace Model { class Hubbard; }

namespace PQMC {

    struct PqmcParams;

    using ScalarType     = double;
    using GreensFunction = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using ProjectionMat  = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
    using SvdStack       = Utils::SvdStack<ScalarType>;
    using timeIndex      = int;

    class PqmcEngine {
        public:
            double m_theta{};
            double m_dt{};
            int    m_nt{};
            int    m_nl{};
            int    m_ns{};
            int    m_np{};
            int    m_stabilization_pace{};
            int    m_current_time_slice{};

            GreensFunction* m_green_tt_up{};
            GreensFunction* m_green_tt_dn{};
            GreensFunction* m_green_t0_up{};
            GreensFunction* m_green_t0_dn{};
            GreensFunction* m_green_0t_up{};
            GreensFunction* m_green_0t_dn{};

            ProjectionMat m_projection_mat{};

            SvdStack* m_svd_stack_left_up{};
            SvdStack* m_svd_stack_left_dn{};
            SvdStack* m_svd_stack_right_up{};
            SvdStack* m_svd_stack_right_dn{};

            double m_wrap_error{};

            void initial( const PqmcParams& params, const Model::Hubbard& model );

            void metropolis_update( Model::Hubbard& model, timeIndex t );

            void wrap_from_0_to_2theta( Model::Hubbard& model, timeIndex t );
            void wrap_from_2theta_to_0( Model::Hubbard& model, timeIndex t );

            void sweep_from_0_to_2theta( Model::Hubbard& model );
            void sweep_from_2theta_to_0( Model::Hubbard& model );

    };

} // namespace PQMC

#endif // PQMC_ENGINE_H