/*
 *   measure_methods.hpp
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 *   This file includes definition of measuring methods for sampling the physical observables using QMC.
 *   Both equal-time and dynamic measurements for fermionic ( or bosonic, if needed ) observables are supported.
 */

#ifndef MEASURE_METHODS_H
#define MEASURE_METHODS_H
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

            static void measure_double_occupation ( Observable::ScalarObs& double_occupation,
                                                    const Measure::MeasureHandler& meas_handler,
                                                    const PQMC::PqmcEngine& walker,
                                                    const Model::Hubbard& model                  )
            {
                
            }

            static void measure_kinetic_energy    ( Observable::ScalarObs& kinetic_energy,
                                                    const Measure::MeasureHandler& meas_handler,
                                                    const PQMC::PqmcEngine& walker,
                                                    const Model::Hubbard& model                  )
            {
                
            }


            // Dynamic Measurements:
            //    1. Momentum-resolved dynamic green's functions: G(k,t) = < c(k,t) * c^+(k,0) >
            //    2. Density of states on imaginary-time grids: D(t) = 1/N \sum i < c(i,t) * c^+(i,0) >

            static void measure_greens_functions  ( Observable::MatrixObs& greens_functions, 
                                                    const Measure::MeasureHandler& meas_handler,
                                                    const PQMC::PqmcEngine& walker,
                                                    const Model::Hubbard& model                  )
            {

            }

            static void measure_density_of_states ( Observable::VectorObs& density_of_states, 
                                                    const Measure::MeasureHandler& meas_handler,
                                                    const PQMC::PqmcEngine& walker,
                                                    const Model::Hubbard& model                  )
            {

            }

    };

} // namespace Measure

#endif // MEASURE_METHODS_H