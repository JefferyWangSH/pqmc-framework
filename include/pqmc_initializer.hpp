/*
 *   pqmc_initializer.hpp
 * 
 *     Created on: Apr 15, 2023
 *         Author: Jeffery Wang
 * 
 *   This file defines the pqmc initializer class PQMC::PqmcInitializer,
 *   containing the static member functions to initialize the pqmc modules as a whole. 
 */

#ifndef PQMC_INITIALIZER_HPP
#define PQMC_INITIALIZER_HPP
#pragma once

#include <string_view>
#include "model.h"
#include "pqmc_engine.h"
#include "measure_handler.h"
#include "pqmc_params.hpp"
#include "pqmc_io.hpp"
#include "utils/toml.hpp"

namespace PQMC {

    class PqmcInitializer {

        public:

            // parse simulation parmameters from the toml configuration file
            static void parse_toml_config( std::string_view toml_config, int world_size, PqmcParams& params )
            {
                // parse the configuration file
                auto config = toml::parse_file( toml_config );

                // ---------------------------------------------------------------------------------------------------------------------------
                //                                                  Parse the Model
                // ---------------------------------------------------------------------------------------------------------------------------
                
                params.t = config["HubbardModel"]["t"].value_or(1.0);
                params.u = config["HubbardModel"]["u"].value_or(4.0);
                params.nl = config["SquareLattice"]["nl"].value_or(4);
                params.ns = params.nl * params.nl;
                params.np = config["SquareLattice"]["np"].value_or(16);
                if ( params.np % 2 != 0 ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "particle number should be even. ( su2 symmetric required )"
                              << std::endl;
                    exit(1);
                }

                // todo: read momentum, momentum_list

                // ---------------------------------------------------------------------------------------------------------------------------
                //                                                Parse the PQMC Engine
                // ---------------------------------------------------------------------------------------------------------------------------

                params.theta = config["MonteCarlo"]["theta"].value_or(5.0);
                params.dt = config["MonteCarlo"]["dt"].value_or(0.125);
                params.nt = config["MonteCarlo"]["nt"].value_or(80);
                params.stabilization_pace = config["MonteCarlo"]["stabilization_pace"].value_or(10);
                if ( params.dt * params.nt != 2 * params.theta ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "conflicts in projection length `theta`, imag-time slices `nt` and imag-time spacing `dt`. ( dt * nt = 2 theta )"
                              << std::endl;
                    exit(1);
                }
                if ( params.nt % 2 != 0 ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "imag-time slices `nt` should be even. ( theta/dt = nt/2 belongs to Z )"
                              << std::endl;
                    exit(1);
                }

                params.beta = config["Measure"]["beta"].value_or(1.25);
                params.ntm = config["Measure"]["ntm"].value_or(20);
                if ( params.dt * params.ntm != 2 * params.beta ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "conflicts in measurement window `beta`, imag-time slices `ntm` and imag-time spacing `dt`. ( dt * ntm = 2 beta )"
                              << std::endl;
                    exit(1);
                }
                if ( params.ntm % 2 != 0 ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "imag-time slices `ntm` in measurement window should be even. ( beta/dt = ntm/2 belongs to Z )"
                              << std::endl;
                    exit(1);
                }
                if ( params.beta > params.theta ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "measurement window `beta` exceeds the maximum of projection length `theta`."
                              << std::endl;
                    exit(1);    
                }
                if ( params.ntm > params.nt ) {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): " 
                              << "imag-time slices `ntm` in measurement window exceeds the total slice number `nt`."
                              << std::endl;
                    exit(1);    
                }
                
                // ---------------------------------------------------------------------------------------------------------------------------
                //                                                Parse the Measure Handler
                // ---------------------------------------------------------------------------------------------------------------------------

                params.sweeps_warmup = config["Measure"]["sweeps_warmup"].value_or(100);
                // distribute the measuring tasks to a set of processors
                params.bin_num = std::ceil( config["Measure"]["bin_num"].value_or(50)/world_size );
                params.bin_capacity = config["Measure"]["bin_capacity"].value_or(50);
                params.sweeps_between_bins = config["Measure"]["sweeps_between_bins"].value_or(50);

                // parse obervable lists
                std::vector<std::string> observables;
                toml::array* observable_arr = config["Measure"]["observables"].as_array();
                if ( observable_arr && observable_arr->is_homogeneous<std::string>() ) {
                    observables.reserve( observable_arr->size() );
                    for ( auto&& el : *observable_arr ) {
                        observables.emplace_back( el.value_or("") );
                    }
                }
                else {
                    std::cerr << "PQMC::PqmcInitializer::parse_toml_config(): "
                              << "undefined observables, check the config." << std::endl;
                    exit(1);
                }
                // deal with special keywords ( all/All , none/None )
                if ( observables.size() == 1 ) {
                    if ( observables[0] == "all" || observables[0] == "All" ) {
                        observables = Measure::MeasureHandler::ObservableAll;
                    }
                    else if ( observables[0] == "none" || observables[0] == "None" ) {
                        observables = {};
                    }
                }
                params.observables = observables;
            }


            static void initialize_pqmc( PqmcEngine& engine, Model::Hubbard& model, Measure::MeasureHandler& handler,
                                         const PqmcParams& params, std::string field_config = "" )
            {
                model.initial( params );
                if ( field_config.empty() ) { model.randomly_initial_ising_fields(); }  // randomly initialize
                else { PqmcIO::read_ising_fields_from_file( field_config, model ); }    // read field configs from file

                handler.initial( params );
                engine.initial( params, model, handler );
            }

    };

} // namespace PQMC

#endif // PQMC_INITIALIZER_HPP