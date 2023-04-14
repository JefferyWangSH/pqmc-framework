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

#include <vector>
#include <string>

namespace PQMC {

    struct PqmcParams {

        int nl{};                               // linear size of lattice
        int ns{};                               // number of lattice sites
        int np{};                               // number of particles
        
        int nt{};                               // imaginary-time slices
        double dt{};                            // imaginary-time interval
        double theta{};                         // half of total projection length in imaginary-time direction ( 2theta = nt*dt )

        double t{};                             // nearest hopping
        double u{};                             // on-site Hubbard interaction

        int stabilization_pace{};               // pace of stabilization

        int sweeps_warmup{};                    // MC sweeps for thermalization
        int bin_num{};                          // number of bins
        int bin_capacity{};                     // capacity ( MC samples ) in one bin
        int sweeps_between_bins{};              // MC sweeps between neighbouring bins, avoiding correlation among bins

        double beta{};                          // imaginary-time length within which measurements are performed [ theta-beta, theta+beta ]
        int ntm{};                              // imaginary-time slices for measurements
        
        std::vector<std::string> obs_list{};    // list of observables
        int momentum{};                         // ( index of ) momentum point
        std::vector<int> momentum_list{};       // ( indices of ) momentum list

    };

} // namespace PQMC

#endif // PQMC_PARAMS_HPP