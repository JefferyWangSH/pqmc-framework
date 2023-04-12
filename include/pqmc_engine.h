#ifndef PQMC_ENGINE_H
#define PQMC_ENGINE_H
#pragma once

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "utils/svd_stack.h"

namespace Model { class Hubbard; }

namespace PQMC {

    using GreensFunction = Eigen::MatrixXd;
    using ProjectionMat  = Eigen::MatrixXd;

    class PqmcEngine {
        public:
            double m_theta{2};
            double m_dt{0.05};
            int    m_nt{80};
            int    m_nl{4};
            int    m_ns{4*4};
            int    m_particle_num{16};
            int    m_current_time_slice{0};
            int    m_stablization_pace{10};

            GreensFunction m_green_tt_up{};
            GreensFunction m_green_tt_dn{};
            GreensFunction m_green_t0_up{};
            GreensFunction m_green_t0_dn{};
            GreensFunction m_green_0t_up{};
            GreensFunction m_green_0t_dn{};

            ProjectionMat m_projection_mat{};

            Utils::SvdStackReal* m_svd_stack_left_up{};
            Utils::SvdStackReal* m_svd_stack_left_dn{};
            Utils::SvdStackReal* m_svd_stack_right_up{};
            Utils::SvdStackReal* m_svd_stack_right_dn{};

            void initial( const Model::Hubbard& model );
            
    };

} // namespace PQMC

#endif // PQMC_ENGINE_H