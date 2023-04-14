/*
 *   pqmc_main.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "model.h"
#include "pqmc_engine.h"
#include "pqmc_walker.h"
#include "measure_handler.h"
#include "pqmc_params.hpp"
#include "random.h"

#include <iostream>

// the main program
int main( int argc, char* argv[] ) {

    Model::Hubbard* hubbard = new Model::Hubbard();
    PQMC::PqmcEngine* engine = new PQMC::PqmcEngine();
    Measure::MeasureHandler* handler = new Measure::MeasureHandler();
    PQMC::PqmcParams* params = new PQMC::PqmcParams();

    params->nl = 4;
    params->ns = params->nl*params->nl;
    params->np = params->nl*params->nl;    // half-filling
    params->nt = 80;
    params->dt = 0.05;
    params->theta = 2.0;
    params->t = 1.0;
    params->u = 4.0;
    params->stabilization_pace = 10;
    params->sweeps_warmup = 1000;
    params->bin_num = 100;
    params->bin_capacity = 20;
    params->sweeps_between_bins = 20;
    params->obs_list = { "DoubleOccupation", };
    params->momentum = 0;
    params->momentum_list = { 0 };
    params->beta = 0.5;
    params->ntm = 20;

    Utils::Random::set_seed( time(nullptr) );

    hubbard->initial( *params );
    hubbard->randomly_initial_ising_fields();
    handler->initial( *params );
    engine->initial( *params, *hubbard, *handler );    
    
    PQMC::PqmcWalker::thermalize( *engine, *hubbard, *handler );
    PQMC::PqmcWalker::measure( *engine, *hubbard, *handler );
    PQMC::PqmcWalker::analyse( *handler );
    
    std::cout << engine->m_wrap_error << std::endl;
    std::cout << handler->find<Observable::ScalarObs>( "DoubleOccupation" ).description() << std::endl;
    std::cout << handler->find<Observable::ScalarObs>( "DoubleOccupation" ).mean_value() << std::endl;
    std::cout << handler->find<Observable::ScalarObs>( "DoubleOccupation" ).std_error() << std::endl;
    
    return 0;
}