/*
 *   measure_methods.hpp
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 *   This file includes definition of measuring methods for sampling the physical observables using QMC.
 *   Both equal-time and dynamic measurements for fermionic ( or bosonic, if needed ) observables are supported.
 */

#ifndef MEASURE_METHODS_HPP
#define MEASURE_METHODS_HPP
#pragma once

#include "observable.h"
#include "pqmc_engine.h"
#include "model.h"
#include "measure_handler.h"

namespace Measure {

    class Methods {

        public:

            // definitions of the measuring methods

            // Equal-time Measurements:
            //    1. Double Occupation D
            //    2. Kinetic Energy Ek

            static void measure_double_occupation( Observable::ScalarObs& double_occupation,
                                                   const Measure::MeasureHandler& handler,
                                                   const PQMC::PqmcEngine& engine,
                                                   const Model::Hubbard& model                  )
            {
                const int ntm = engine.m_ntm;
                const int ns = engine.m_ns;

                for ( auto t = 0; t < ntm; ++t ) {
                    const PQMC::GreensFunction& gu = (*engine.m_vec_green_tt_up)[t];
                    const PQMC::GreensFunction& gd = (*engine.m_vec_green_tt_dn)[t];

                    double temp_double_occupation = 0.0;
                    for ( auto i = 0; i < ns; ++i ) {
                        temp_double_occupation += ( 1.0 - gu(i,i) ) * ( 1.0 - gd(i,i) );
                    }
                    double_occupation.temp_value() += temp_double_occupation / ns;
                    ++double_occupation;
                }
            }

            static void measure_kinetic_energy( Observable::ScalarObs& kinetic_energy,
                                                const Measure::MeasureHandler& handler,
                                                const PQMC::PqmcEngine& engine,
                                                const Model::Hubbard& model                  )
            {
                const int ntm = engine.m_ntm;
                const int nl = engine.m_nl;
                const int ns = engine.m_ns;

                for ( auto t = 0; t < ntm; ++t ) {
                    const PQMC::GreensFunction& gu = (*engine.m_vec_green_tt_up)[t];
                    const PQMC::GreensFunction& gd = (*engine.m_vec_green_tt_dn)[t];

                    double temp_kinetic_energy = 0.0;
                    for ( auto x = 0; x < nl; ++x ) {
                        for ( auto y = 0; y < nl; ++y ) {
                            const auto i = x + y * nl;
                            const auto xplus1 = (x+1) % nl + y * nl;
                            const auto yplus1 = x + ((y+1) % nl) * nl;
                            temp_kinetic_energy += gu(i,xplus1) + gd(i,xplus1) + gu(i,yplus1) + gd(i,yplus1);
                        }
                    }
                    kinetic_energy.temp_value() += 2 * model.m_t * temp_kinetic_energy / ns;
                    ++kinetic_energy;
                }
            }


            // Dynamic Measurements:
            //    1. Momentum-resolved dynamic green's functions: G(k,t) = < c(k,t) * c^+(k,0) >
            //    2. Density of states on imaginary-time grids: D(t) = 1/N \sum i < c(i,t) * c^+(i,0) >

            static void measure_greens_functions( Observable::MatrixObs& greens_functions, 
                                                  const Measure::MeasureHandler& handler,
                                                  const PQMC::PqmcEngine& engine,
                                                  const Model::Hubbard& model                  )
            {

            }

            static void measure_density_of_states( Observable::MatrixObs& density_of_states, 
                                                   const Measure::MeasureHandler& handler,
                                                   const PQMC::PqmcEngine& engine,
                                                   const Model::Hubbard& model                  )
            {
                const int ntm = engine.m_ntm;
                const int ns = engine.m_ns;

                for (auto t = 0; t < ntm; ++t) {
                    const PQMC::GreensFunction& gt0up = (*engine.m_vec_green_t0_up)[t];
                    const PQMC::GreensFunction& gt0dn = (*engine.m_vec_green_t0_dn)[t];
                    density_of_states.temp_value()(0,t) += gt0up.trace() / ns;
                    density_of_states.temp_value()(1,t) += gt0dn.trace() / ns;
                }
                ++density_of_states;
            }

    };

} // namespace Measure

#endif // MEASURE_METHODS_HPP