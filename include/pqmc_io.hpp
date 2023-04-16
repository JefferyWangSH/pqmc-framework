/*
 *   pqmc_io.hpp
 * 
 *     Created on: Apr 15, 2023
 *         Author: Jeffery Wang
 * 
 *   This file defines PQMC::PqmcIO class 
 *   containing the basic IO interfaces for the input/output of QMC data
 */

#ifndef PQMC_IO_HPP
#define PQMC_IO_HPP
#pragma once

#include <fstream>
#include <string>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "model.h"
#include "pqmc_engine.h"
#include "pqmc_walker.h"
#include "observable.h"
#include "measure_handler.h"
#include "pqmc_params.hpp"

namespace PQMC {

    class PqmcIO {
        public:

            template< typename StreamType >
            static void output_initialization_info( StreamType& ostream, int world_size, const PqmcParams& params, const Measure::MeasureHandler& handler )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_init_info(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // output formats
                    boost::format fmt_param_str("%| 30s|%| 7s|%| 24s|\n");
                    boost::format fmt_param_int("%| 30s|%| 7s|%| 24d|\n");
                    boost::format fmt_param_double("%| 30s|%| 7s|%| 24.3f|\n");
                    const std::string_view joiner = "->";
                    auto bool2str = [](bool b) {if (b) return "True"; else return "False";};

                    // -----------------------------------------------------------------------------------------------------------------
                    //                                          Output model information
                    // -----------------------------------------------------------------------------------------------------------------

                    ostream << "   Model: Repulsive Hubbard\n\n"
                            << fmt_param_double % "Nearest hopping \'t\'" % joiner % params.t
                            << fmt_param_double % "Hubbard interaction \'U\'" % joiner % params.u
                            << std::endl;
                    
                    // -----------------------------------------------------------------------------------------------------------------
                    //                                         Output lattice information
                    // -----------------------------------------------------------------------------------------------------------------

                    boost::format fmt_cell("%d * %d");
                    boost::format fmt_momentum("(%.2f, %.2f) pi");

                    ostream << "   Lattice: Square lattice\n\n"
                            << fmt_param_str % "Size of cell" % joiner % ( fmt_cell % params.nl % params.nl )
                            << fmt_param_int % "Particle number" % joiner % params.np
                            << std::endl;

                    // -----------------------------------------------------------------------------------------------------------------
                    //                                          Output MonteCarlo Params
                    // -----------------------------------------------------------------------------------------------------------------
                    ostream << "   MonteCarlo Params:\n\n"
                            << fmt_param_double % "Projection length 'theta'" % joiner % params.theta
                            << fmt_param_int % "Imag-time slices" % joiner % params.nt
                            << fmt_param_double % "Imag-time spacing" % joiner % params.dt
                            << fmt_param_int % "Stabilization pace" % joiner % params.stabilization_pace
                            << std::endl;

                    // -----------------------------------------------------------------------------------------------------------------
                    //                                           Output Measuring Params
                    // -----------------------------------------------------------------------------------------------------------------
                    ostream << "   Measuring Params:\n\n"
                            << fmt_param_double % "Measurement window 'beta'" % joiner % params.beta
                            << fmt_param_int % "Imag-time slices for Meas." % joiner % params.ntm
                            << fmt_param_str % "Warm up" % joiner % bool2str( handler.isWarmUp() )
                            << fmt_param_str % "Equal-time measure" % joiner % bool2str( handler.isEqualTime() )
                            << fmt_param_str % "Dynamical measure" % joiner % bool2str( handler.isDynamic() )
                            << std::endl;
                    
                    ostream << fmt_param_int % "Sweeps for warmup" % joiner % params.sweeps_warmup
                            << fmt_param_int % "Number of bins" % joiner % ( params.bin_num * world_size )
                            << fmt_param_int % "MC Samples per bin" % joiner % params.bin_capacity
                            << fmt_param_int % "Sweeps between bins" % joiner % params.sweeps_between_bins
                            << std::endl;
                }
            }


            template< typename StreamType >
            static void output_ending_info( StreamType& ostream, const PqmcEngine& engine )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_ending_info(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // parse the duration
                    const int day = std::floor((double)PqmcWalker::timer() / 86400000);
                    const int hour = std::floor(((double)PqmcWalker::timer()/1000 - day * 86400) / 3600);
                    const int minute = std::floor(((double)PqmcWalker::timer()/1000 - day * 86400 - hour * 3600) / 60);
                    const double sec = (double)PqmcWalker::timer()/1000 - 86400 * day - 3600 * hour - 60 * minute;

                    // output the time cost of simulation
                    if ( day ) { ostream << boost::format("\n>> The simulation finished in %d d %d h %d m %.2f s.\n") % day % hour % minute % sec << std::endl; }
                    else if ( hour ) { ostream << boost::format("\n>> The simulation finished in %d h %d m %.2f s.\n") % hour % minute % sec << std::endl; }
                    else if ( minute ) { ostream << boost::format("\n>> The simulation finished in %d m %.2f s.\n") % minute % sec << std::endl; }
                    else { ostream << boost::format("\n>> The simulation finished in %.2f s.\n") % sec << std::endl; }

                    // output wrapping errors of the evaluations of Green's functions
                    ostream << boost::format(">> Maximum of the wrapping error: %.5e\n") % engine.m_wrap_error << std::endl;
                }
            }


            template< typename StreamType, typename ObsType >
            static void output_observable( StreamType& ostream, const Observable::Observable<ObsType>& obs )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_observable(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // standard screen output
                    if constexpr ( std::is_same_v<StreamType, std::ostream> ) {
                        // for scalar observables
                        if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                            boost::format fmt_scalar_obs("%| 28s|%| 7s|%| 20.12f|  pm  %.12f");
                            const std::string joiner = "->";
                            ostream << fmt_scalar_obs % obs.description() % joiner % obs.mean_value() % obs.std_error() << std::endl;
                        }
                        // for vector observables
                        else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                            // todo: currently not used
                        }
                        // for matrix observables
                        else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                            // todo
                        }
                        // other observable type, raising errors
                        else {
                            std::cerr << "PQMC::PqmcIO::output_observable(): "
                                      << "undefined observable type." << std::endl;
                            exit(1);
                        }
                    }

                    // standard file output
                    else if constexpr ( std::is_same_v<StreamType, std::ofstream> ) {
                        
                        // for scalar observables
                        if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                            // for specfic scalar observable, output the mean value, std_error and relative error in order. 
                            boost::format fmt_scalar_obs("%| 20.10f|%| 20.10f|%| 20.10f|");
                            ostream << fmt_scalar_obs % obs.mean_value() % obs.std_error() % (obs.std_error()/obs.mean_value()) << std::endl;
                        }

                        // for vector observables
                        else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                            // output vector observable
                            boost::format fmt_size_info("%| 20d|");
                            boost::format fmt_vector_obs("%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

                            const int size = obs.mean_value().size();
                            const auto relative_error = (obs.std_error().array()/obs.mean_value().array()).matrix();
                            ostream << fmt_size_info % size << std::endl;
                            for ( auto i = 0; i < size; ++i ) {
                                // output the mean value, error bar and relative error in order. 
                                ostream << fmt_vector_obs % i % obs.mean_value()(i) % obs.std_error()(i) % relative_error(i) << std::endl;
                            }
                        }

                        // for matrix observables
                        else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                            // output matrix observable 
                            boost::format fmt_size_info("%| 20d|%| 20d|");
                            boost::format fmt_matrix_obs("%| 20d|%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

                            const int row = obs.mean_value().rows();
                            const int col = obs.mean_value().cols();
                            const auto relative_error = (obs.std_error().array()/obs.mean_value().array()).matrix();
                            ostream << fmt_size_info % row % col << std::endl;
                            for ( auto i = 0; i < row; ++i ) {
                                for ( auto j = 0; j < col; ++j ) {
                                    // output the mean value, error bar and relative error in order. 
                                    ostream << fmt_matrix_obs % i % j % obs.mean_value()(i,j) % obs.std_error()(i,j) % relative_error(i,j) << std::endl;
                                }
                            }
                        }

                        // other observable types, raising errors
                        else {
                            std::cerr << "PQMC::PqmcIO::output_observable(): "
                                      << "undefined observable type." << std::endl;
                            exit(1);
                        }
                    }

                    // others stream types, raising errors
                    else {
                        std::cerr << "PQMC::PqmcIO::output_observable(): "
                                  << "unsupported type of output stream." << std::endl;
                        exit(1);
                    }
                }
            }


            template< typename StreamType, typename ObsType >
            static void output_observable_in_bins( StreamType& ostream, const Observable::Observable<ObsType>& obs, bool show_header = true )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_observable_in_bins(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // for scalar observables
                    if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                        // output bin data of scalar observable
                        boost::format fmt_size_info("%| 20d|");
                        boost::format fmt_scalar_obs("%| 20d|%| 20.10f|");

                        const int number_of_bins = obs.bin_num();
                        if ( show_header ) { ostream << fmt_size_info % number_of_bins << std::endl; }
                        for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                            ostream << fmt_scalar_obs % bin % obs.bin_data(bin) << std::endl;
                        }
                    }

                    // for vector observables
                    else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                        // output bin data of vector observable
                        boost::format fmt_size_info("%| 20d|%| 20d|");
                        boost::format fmt_vector_obs("%| 20d|%| 20d|%| 20.10f|");

                        const int number_of_bins = obs.bin_num();
                        const int size = obs.mean_value().size();
                        if ( show_header ) { ostream << fmt_size_info % number_of_bins % size << std::endl; }
                        for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                            for ( auto i = 0; i < size; ++i ) {
                                ostream << fmt_vector_obs % bin % i % obs.bin_data(bin)(i) << std::endl;
                            }
                        }
                    }

                    // for matrix observables
                    else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                        // output bin data of matrix observable
                        boost::format fmt_size_info("%| 20d|%| 20d|%| 20d|");
                        boost::format fmt_matrix_obs("%| 20d|%| 20d|%| 20d|%| 20.10f|");

                        const int number_of_bins = obs.bin_num();
                        const int row = obs.mean_value().rows();
                        const int col = obs.mean_value().cols();
                        if ( show_header ) { ostream << fmt_size_info % number_of_bins % row % col << std::endl; }
                        for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                            for ( auto i = 0; i < row; ++i ) {
                                for ( auto j = 0; j < col; ++j ) {
                                    ostream << fmt_matrix_obs % bin % i % j % obs.bin_data(bin)(i,j) << std::endl;
                                }
                            }
                        }
                    }

                    // other observable types, raising errors
                    else {
                        std::cerr << "PQMC::PqmcIO::output_observable_in_bins(): "
                                  << "undefined observable type." << std::endl;
                        exit(1);
                    }
                }
            }


            template< typename StreamType >
            static void output_k_stars( StreamType& ostream )
            {
                // todo
            }


            template< typename StreamType >
            static void output_imaginary_time_grids( StreamType& ostream, const PqmcEngine& engine )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_imaginary_time_grids(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // output the imaginary-time grids within measurement window
                    boost::format fmt_tgrids_info("%| 20d|%| 20.5f|%| 20.5f|");
                    boost::format fmt_tgrids("%| 20d|%| 20.10f|");
                    ostream << fmt_tgrids_info % engine.m_ntm % (2*engine.m_beta) % engine.m_dt << std::endl;
                    for ( auto t = 0; t < engine.m_ntm; ++t ) {
                        ostream << fmt_tgrids % t % ( t * engine.m_dt ) << std::endl;
                    }
                }
            }


            template< typename StreamType >
            static void output_ising_fields( StreamType& ostream, const Model::Hubbard& model )
            {
                if ( !ostream ) {
                    std::cerr << "PQMC::PqmcIO::output_ising_fields(): "
                              << "the ostream failed to work, check the input." << std::endl;
                    exit(1);
                }
                else {
                    // output current configuration of auxiliary ising fields
                    boost::format fmt_fields_info("%| 20d|%| 20d|");
                    boost::format fmt_fields("%| 20d|%| 20d|%| 20d|");
                    const int nt = model.m_nt;
                    const int ns = model.m_ns;

                    ostream << fmt_fields_info % nt % ns << std::endl;
                    for ( auto t = 0; t < nt; ++t ) {
                        for ( auto i = 0; i < ns; ++i ) {
                            ostream << fmt_fields % t % i % model.m_ising_fields(t,i) << std::endl;
                        }
                    }
                }
            }


            static void read_ising_fields_from_file( std::string filename, Model::Hubbard& model )
            {
                std::ifstream infile( filename, std::ios::in );

                // check whether the ifstream works well
                if ( !infile.is_open() ) {
                    std::cerr << "PQMC::PqmcIO::read_ising_fields_from_file(): "
                              << "fail to open file \'" << filename << "\'." << std::endl;
                    exit(1);
                }

                // temporary parameters
                std::string line;
                std::vector<std::string> data;

                // consistency check of the model parameters
                // read the first line which contains the model information
                getline( infile, line );
                boost::split( data, line, boost::is_any_of(" "), boost::token_compress_on );
                data.erase( std::remove( std::begin(data), std::end(data), "" ), std::end(data) );
                const int nt = boost::lexical_cast<int>(data[0]);
                const int ns = boost::lexical_cast<int>(data[1]);
                if ( ( nt != model.m_nt ) || ( ns != model.m_ns ) ) {
                    std::cerr << "PQMC::PqmcIO::read_ising_fields_from_file(): "
                              << "inconsistency between model settings and input configs (time or space size). " 
                              << std::endl;
                    exit(1);
                }

                // read in the configurations of ising fields
                int t, i;
                while( getline( infile, line ) ) {
                    boost::split( data, line, boost::is_any_of(" "), boost::token_compress_on );
                    data.erase( std::remove( std::begin(data), std::end(data), "" ), std::end(data) );
                    t = boost::lexical_cast<int>( data[0] );
                    i = boost::lexical_cast<int>( data[1] );
                    model.m_ising_fields(t, i) = boost::lexical_cast<double>( data[2] );
                }

                // close the file stream
                infile.close();
            }


    };

} // namespace PQMC

#endif // PQMC_IO_HPP