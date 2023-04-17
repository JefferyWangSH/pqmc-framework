/*
 *   pqmc_main.cpp
 * 
 *     Created on: Apr 11, 2023
 *         Author: Jeffery Wang
 * 
 */

#include <string>
#include <iostream>
#include <fstream>

#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/program_options.hpp>

#include "model.h"
#include "pqmc_engine.h"
#include "pqmc_walker.h"
#include "measure_handler.h"
#include "pqmc_params.hpp"
#include "pqmc_initializer.hpp"
#include "pqmc_io.hpp"
#include "random.h"
#include "utils/mpi.hpp"

// the main program
int main( int argc, char* argv[] ) {

    // -------------------------------------------------------------------------------------------------------------
    //                                         Initialize MPI environment
    // -------------------------------------------------------------------------------------------------------------
    
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    const int master = 0;
    const int rank = world.rank();

    
    // -------------------------------------------------------------------------------------------------------------
    //                                               Program options
    // -------------------------------------------------------------------------------------------------------------
    
    std::string config_file;
    std::string output_path;
    std::string ising_fields_file;
    
    // read parameters from the command line 
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
        (   "help,h", "display this information" )
        (   "config,c", 
            boost::program_options::value<std::string>(&config_file)->default_value("../example/config.toml"), 
            "path of the toml configuration file, default: ../example/config.toml" )
        (   "output,o",
            boost::program_options::value<std::string>(&output_path)->default_value("../example"), 
            "folder path which stores the output of measuring results, default: ../example" )
        (   "ising-fields,f",
            boost::program_options::value<std::string>(&ising_fields_file), 
            "path of the configurations of ising fields, if not assigned the fields are to be set randomly." );
    
    // parse the command line options
    try {
        boost::program_options::store( parse_command_line( argc, argv, opts ), vm );
    }
    catch (...) {
        std::cerr << "main(): undefined options got from command line." << std::endl;
        exit(1);
    }
    boost::program_options::notify(vm);

    // show the helping messages
    if (vm.count("help")) {
        std::cerr << argv[0] << "\n" << opts << std::endl;
        return 0;
    }

    // initialize the output folder, create if not exist
    if ( rank == master ) {
        if ( access(output_path.c_str(), 0) != 0 ) {
            const std::string command = "mkdir -p " + output_path;
            if ( system(command.c_str()) != 0 ) {
                std::cerr << boost::format("main(): fail to creat folder at %s .") % output_path 
                          << std::endl;
                exit(1);
            }
        }
    }

    
    // -------------------------------------------------------------------------------------------------------------
    //                                        Output current date and time
    // -------------------------------------------------------------------------------------------------------------

    if ( rank == master ) {
        const auto current_time = boost::posix_time::second_clock::local_time();
        std::cout << boost::format(">> Current time: %s\n") % current_time << std::endl;
    }

    
    // -------------------------------------------------------------------------------------------------------------
    //                                        Output MPI and hardware info
    // -------------------------------------------------------------------------------------------------------------

    if ( rank == master ) {
        // print the MPI and hardware information
        boost::format fmt_mpi(">> Distribute tasks to %s processors, with the master processor being %s.\n");
        std::cout << fmt_mpi % world.size() % env.processor_name() << std::endl;
    }

    
    // -------------------------------------------------------------------------------------------------------------
    //                                              PQMC simulation
    // -------------------------------------------------------------------------------------------------------------

    // set up random seeds for different processes
    Utils::Random::set_seed( std::time(nullptr) + rank );

    // --------------------------------------------  Initializations  ----------------------------------------------

    // create PQMC module objects
    Model::Hubbard* model = new Model::Hubbard();
    PQMC::PqmcEngine* engine = new PQMC::PqmcEngine();
    Measure::MeasureHandler* handler = new Measure::MeasureHandler();
    PQMC::PqmcParams* params = new PQMC::PqmcParams();

    // parse parmas from the configuation file
    PQMC::PqmcInitializer::parse_toml_config( config_file, world.size(), *params );
    
    // initialize PQMC, preparing for the simulation
    PQMC::PqmcInitializer::initialize_pqmc( *engine, *model, *handler, *params, ising_fields_file );
    
    // output the initialization info
    if ( rank == master ) {
        if ( !ising_fields_file.empty() ) { std::cout << ">> Configurations of ising fields read from file.\n" << std::endl; }
        else { std::cout << ">> Configurations of ising fields initialized randomly.\n" << std::endl; }
        std::cout << ">> Initialization finished. \n\n" 
                  << ">> The simulation is going to get started with parameters shown below:\n"
                  << std::endl;
        PQMC::PqmcIO::output_initialization_info( std::cout, world.size(), *params, *handler );
    }

    // set up progress bar
    PQMC::PqmcWalker::show_progress_bar( (rank == master) );
    PQMC::PqmcWalker::progress_bar_format( 60, '=', ' ' );
    PQMC::PqmcWalker::set_refresh_rate( 10 );

    // ----------------------------------------  Crucial simulation steps  -----------------------------------------
    
    // the PQMC simulation start
    PQMC::PqmcWalker::timer_begin();
    PQMC::PqmcWalker::thermalize( *engine, *model, *handler );
    PQMC::PqmcWalker::measure( *engine, *model, *handler );
    
    // gather observable objects from other processes
    Utils::MPI::mpi_gather( world, *handler );

    // perform the analysis
    PQMC::PqmcWalker::analyse( *handler );

    // end the timer
    PQMC::PqmcWalker::timer_end();

    // output the ending info
    if ( rank == master ) {
        PQMC::PqmcIO::output_ending_info( std::cout, *engine );
    }
    
    // ----------------------------------------  Output measuring results  -----------------------------------------

    if ( rank == master ) {

        std::ofstream outfile;
        const auto ising_fields_out = ( ising_fields_file.empty() )? output_path + "/ising.fields" : ising_fields_file;
        outfile.open( ising_fields_out, std::ios::trunc );
        PQMC::PqmcIO::output_ising_fields( outfile, *model );
        outfile.close();

        outfile.open( output_path + "/tgrids.out", std::ios::trunc );
        PQMC::PqmcIO::output_imaginary_time_grids( outfile, *engine );
        outfile.close();

        if ( handler->find( "DoubleOccupation" ) ) {
            PQMC::PqmcIO::output_observable( std::cout, handler->find<Observable::ScalarObs>("DoubleOccupation") );
            outfile.open( output_path + "/double_occu.bins.out", std::ios::trunc );
            PQMC::PqmcIO::output_observable_in_bins( outfile, handler->find<Observable::ScalarObs>("DoubleOccupation"), true );
            outfile.close();
        }        
        
        if ( handler->find( "KineticEnergy" ) ) {
            PQMC::PqmcIO::output_observable( std::cout, handler->find<Observable::ScalarObs>("KineticEnergy") );
            outfile.open( output_path + "/kinetic_energy.bins.out", std::ios::trunc );
            PQMC::PqmcIO::output_observable_in_bins( outfile, handler->find<Observable::ScalarObs>("KineticEnergy"), true );
            outfile.close();
        }

        if ( handler->find( "DensityOfStates" ) ) {
            outfile.open( output_path + "/dos.bins.out", std::ios::trunc );
            PQMC::PqmcIO::output_observable_in_bins( outfile, handler->find<Observable::MatrixObs>("DensityOfStates"), true );
            outfile.close();
        }

        if ( handler->find( "ProjectionBenchmark" ) ) {
            outfile.open( output_path + "/benchmark.out", std::ios::trunc );
            PQMC::PqmcIO::output_observable( outfile, handler->find<Observable::VectorObs>("ProjectionBenchmark") );
            outfile.close();
        }
        
    }
    
    return 0;
}