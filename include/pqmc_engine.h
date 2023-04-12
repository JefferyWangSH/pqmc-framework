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

    using GreensFunction = Eigen::MatrixXd;
    using ProjectionMat  = Eigen::MatrixXd;

    class PqmcEngine {
        public:
            double m_theta{};
            double m_dt{};
            int    m_nt{};
            int    m_nl{};
            int    m_ns{};
            int    m_np{};
            int    m_stabilization_pace{};
            int    m_current_time_slice{0};

            GreensFunction* m_green_tt_up{};
            GreensFunction* m_green_tt_dn{};
            GreensFunction* m_green_t0_up{};
            GreensFunction* m_green_t0_dn{};
            GreensFunction* m_green_0t_up{};
            GreensFunction* m_green_0t_dn{};

            ProjectionMat m_projection_mat{};

            Utils::SvdStackReal* m_svd_stack_left_up{};
            Utils::SvdStackReal* m_svd_stack_left_dn{};
            Utils::SvdStackReal* m_svd_stack_right_up{};
            Utils::SvdStackReal* m_svd_stack_right_dn{};

            double m_wrap_error{};

            void initial( const PqmcParams& params, const Model::Hubbard& model );

    };

} // namespace PQMC

#endif // PQMC_ENGINE_H