/*
 *   pqmc_params.hpp
 * 
 *     Created on: Apr 12, 2023
 *         Author: Jeffery Wang
 * 
 */

#ifndef PQMC_PARAMS_HPP
#define PQMC_PARAMS_HPP
#pragma once

namespace PQMC {

    struct PqmcParams {

        int nl{};                   // linear size of lattice
        int np{};                   // number of particles
        
        int nt{};                   // imaginary-time slices
        double dt{};                // imaginary-time interval
        double theta{};             // half of total projection length in imaginary-time direction ( 2theta = nt*dt )

        double t{};                 // nearest hopping
        double u{};                 // on-site Hubbard interaction

        int stabilization_pace{};   // pace of stabilization

    };

} // namespace PQMC

#endif // PQMC_PARAMS_HPP