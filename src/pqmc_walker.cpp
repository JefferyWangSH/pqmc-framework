#include "pqmc_walker.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Eigenvalues>

#include <iostream>

namespace PQMC {

    void PqmcWalker::initial( const Model::Hubbard& model )
    {
        Eigen::SelfAdjointEigenSolver<ProjectionMat> eigensolver( model.m_K );
        if ( eigensolver.info() != Eigen::Success ) {
            std::cerr << "PQMC::PqmcWalker::initial(): failure in diagonalizing the hopping matrix, check the input." << std::endl;
            exit(1);
        }
        // std::cout << eigensolver.eigenvalues() << std::endl;
        // std::cout << eigensolver.eigenvectors() << std::endl;
        this->m_projection_mat = eigensolver.eigenvectors();
    }

    const GreensFunction PqmcWalker::compute_fresh_green_tt( const Model::Hubbard& model, int t )
    {

    }

} // namespace PQMC