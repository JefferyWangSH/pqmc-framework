/*
 *   pqmc_main.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include <iostream>
#include <Eigen/Dense>
#include "model.h"
#include "pqmc_engine.h"
#include "pqmc_params.hpp"
#include "random.h"
#include "utils/numerical_stable.hpp"
// #include "utils/linear_algebra.hpp"
// #include "utils/svd_stack.h"

int main() {

    Model::Hubbard* hubbard = new Model::Hubbard();
    PQMC::PqmcEngine* engine = new PQMC::PqmcEngine();

    PQMC::PqmcParams* params = new PQMC::PqmcParams();
    params->nl = 2;
    params->np = params->nl*params->nl;    // half-filling
    params->nt = 10;
    params->dt = 0.1;
    params->theta = 0.5;
    params->t = 1.0;
    params->u = 4.0;
    params->stabilization_pace = 1;

    Utils::Random::set_seed( time(nullptr) );
    hubbard->initial( *params );
    // std::cout << hubbard->m_K << std::endl;
    // std::cout << hubbard->m_expK << std::endl;
    
    hubbard->randomly_initial_ising_fields();
    // std::cout << "\n" << hubbard->m_ising_fields << std::endl;


    engine->initial( *params, *hubbard );
    // std::cout << engine->m_projection_mat << std::endl;
    // std::cout << *(engine->m_green_tt_up) << std::endl;
    
//     engine->metropolis_update( *hubbard, 0 );
    hubbard->update_ising_field( params->nt, 0 );
    hubbard->update_greens_function( *engine, params->nt, 0 );
    Eigen::MatrixXd mat1 = *(engine->m_green_tt_up);
    
    delete engine;
    engine = new PQMC::PqmcEngine();
    engine->initial( *params, *hubbard );
    Eigen::MatrixXd mat2 = *(engine->m_green_tt_up);
    
    std::cout << mat1 << std::endl;
    std::cout << "\n" << mat2 << std::endl;
    std::cout << (mat1 - mat2).maxCoeff() << std::endl;


//     Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//     for ( auto t = 0; t <= params->nt; ++t ) {
//         hubbard->multiply_B_from_left( mat, t, +1 );
//     }
//     Eigen::MatrixXd gtt = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl) - engine->m_projection_mat * (engine->m_projection_mat.transpose() * mat * engine->m_projection_mat).inverse() * engine->m_projection_mat.transpose() * mat;
//     std::cout << ( gtt - *(engine->m_green_tt_up) ).maxCoeff() << std::endl;









//     Utils::SvdStackReal left({params->nl*params->nl,params->np/2}, params->nt+1, engine->m_projection_mat);
//     Utils::SvdStackReal right({params->nl*params->nl,params->np/2}, params->nt+1, engine->m_projection_mat);
//     for ( auto t = params->nt; t >= 1; --t ) {
//         Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//         hubbard->multiply_transB_from_left( mat, t, +1 );
//         left.push( mat );
//     }
//     Eigen::MatrixXd mat2 = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//     hubbard->multiply_B_from_left( mat2, 0, +1 );
//     right.push(mat2);

//     Eigen::MatrixXd gtt = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//     Utils::NumericalStable::compute_equaltime_greens<PQMC::PqmcEngine>( left, right, gtt );
    
//     engine->wrap_from_0_to_2theta( *hubbard, 0 );
//     std::cout << ( gtt - *(engine->m_green_tt_up) ).maxCoeff() << std::endl;





//     Utils::SvdStackReal stack1({params->nl*params->nl,params->np/2}, params->nt+1, engine->m_projection_mat);
//     Utils::SvdStackReal stack2({params->nl*params->nl,params->np/2}, params->nt+1, engine->m_projection_mat);
//     for ( auto t = 0; t <= params->nt ; ++t ) {
//         Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//         hubbard->multiply_B_from_left( mat, t, +1 );
//         stack1.push( mat );
//     }
//     for ( auto t = params->nt; t >= 0; --t ) {
//         Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl);
//         hubbard->multiply_transB_from_left( mat, t, +1 );
//         stack2.push( mat );
//     }
//     Eigen::MatrixXd gtt = Eigen::MatrixXd::Identity(params->nl*params->nl,params->nl*params->nl) - engine->m_projection_mat * ( engine->m_projection_mat.transpose() * stack1.MatrixU() * stack1.SingularValues().asDiagonal() * stack1.MatrixV().transpose() ).inverse() * stack2.MatrixV() * stack2.SingularValues().asDiagonal() * stack2.MatrixU().transpose();
//     std::cout << gtt << std::endl;
//     std::cout << (gtt - *(engine->m_green_tt_up)).maxCoeff() << std::endl;


    
//     std::cout << engine->m_svd_stack_right_up->CurrentStackLength() << std::endl;
//     std::cout << engine->m_svd_stack_left_up->CurrentStackLength() << std::endl;
//     std::cout << ( engine->m_svd_stack_right_up->MatrixU() * engine->m_svd_stack_right_up->SingularValues().asDiagonal() * engine->m_svd_stack_right_up->MatrixV().transpose() - engine->m_projection_mat ).maxCoeff() << std::endl;




//     Eigen::MatrixXd mat1 = Eigen::MatrixXd::Random(params->nl*params->nl, params->nl*params->nl);
//     Eigen::MatrixXd mat2 = mat1;
//     hubbard->multiply_B_from_left( mat1, 0, +1 );
//     hubbard->multiply_invB_from_left( mat1, 0, +1 );
//     hubbard->multiply_B_from_right( mat1, 0, +1 );
//     hubbard->multiply_invB_from_right( mat1, 0, +1 );
    
//     hubbard->multiply_transB_from_left( mat1, 0, +1 );
//     mat1.transposeInPlace();
//     hubbard->multiply_invB_from_right( mat1, 0, +1 );
//     mat1.transposeInPlace();
//     std::cout << (mat1 - mat2).maxCoeff() << std::endl;






//     int row = 6;
//     int col = 5;
//     Eigen::MatrixXcd mat = Eigen::MatrixXcd::Random( row, col );
//     Eigen::MatrixXcd u, v;
//     Eigen::VectorXd s;

//     Utils::LinearAlgebra::mkl_lapack_zgesvd( row, col, mat, u, s, v );
    
//     std::cout << u.rows() << " " << u.cols() << std::endl;
//     std::cout << v.rows() << " " << v.cols() << std::endl;
//     std::cout << s.rows() << " " << s.cols() << std::endl;

//     std::cout << u << std::endl;
//     std::cout << "\n" << v << std::endl;
//     std::cout << "\n" << s << std::endl;
    
//     std::cout << u * s.asDiagonal() * v.adjoint() << std::endl;
//     std::cout << "\n" << mat << std::endl;
//     std::cout << "\n" << (mat-u * s.asDiagonal() * v.adjoint()).cwiseAbs2().cwiseSqrt().maxCoeff() << std::endl;


//     int n = 4;
//     Eigen::MatrixXcd mat1 = Eigen::MatrixXcd::Random( n, n );
//     Eigen::MatrixXcd mat2 = mat1 + mat1.adjoint();
//     Eigen::MatrixXcd t;
//     Eigen::VectorXd s;

//     Utils::LinearAlgebra::mkl_lapack_zheev( n, mat2, s, t );

//     std::cout << t.adjoint() * s.asDiagonal() * t << std::endl;
//     std::cout << s << std::endl;
//     std::cout << "\n" << t << std::endl;
//     std::cout << "\n" << mat2 << std::endl;
//     std::cout << (mat2 - t.adjoint() * s.asDiagonal() * t).cwiseAbs2().cwiseSqrt().maxCoeff() << std::endl;





//     std::array<int,2> shape = { 4, 5 };
//     int max_stack_length = 2;
//     Eigen::MatrixXcd P = Eigen::MatrixXcd::Random(shape[0], shape[1]);
//     Utils::SvdStackCpx* stack = new Utils::SvdStackCpx( shape, max_stack_length, P );
    
//     Eigen::MatrixXcd temp = P;
//     Eigen::MatrixXcd mat = Eigen::MatrixXcd::Random( shape[0], shape[0] );
//     temp = mat * mat * temp;
//     stack->push( mat );
//     stack->push( mat );

//     std::cout << stack->MatrixU() * stack->SingularValues().asDiagonal() * stack->MatrixV().adjoint() << std::endl;
//     std::cout << "\n" << temp << std::endl;
//     std::cout << "\n" << stack->CurrentStackLength() << std::endl;
//     std::cout << ( temp - stack->MatrixU() * stack->SingularValues().asDiagonal() * stack->MatrixV().adjoint() ).cwiseAbs2().cwiseSqrt().maxCoeff() << std::endl;

//     stack->pop();
//     temp = mat.inverse() * temp;
//     std::cout << stack->MatrixU() * stack->SingularValues().asDiagonal() * stack->MatrixV().adjoint() << std::endl;
//     std::cout << "\n" << temp << std::endl;
//     std::cout << "\n" << stack->CurrentStackLength() << std::endl;
//     std::cout << ( temp - stack->MatrixU() * stack->SingularValues().asDiagonal() * stack->MatrixV().adjoint() ).cwiseAbs2().cwiseSqrt().maxCoeff() << std::endl;
//     std::cout << "\n" << stack->empty() << std::endl;

//     stack->pop();
//     std::cout << "\n" << stack->CurrentStackLength() << std::endl;
//     std::cout << "\n" << stack->empty() << std::endl;





    return 0;
}