/*
 *   pqmc_walker.cpp
 * 
 *     Created on: Apr 13, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "pqmc_walker.h"
#include "pqmc_engine.h"
#include "model.h"
#include "measure_handler.h"
#include "utils/progress_bar.hpp"

namespace PQMC {

    // definitions of the static members
    bool PqmcWalker::m_show_progress_bar{true};
    unsigned int PqmcWalker::m_progress_bar_width{70};
    unsigned int PqmcWalker::m_refresh_rate{10};
    char PqmcWalker::m_progress_bar_complete_char{'='}, PqmcWalker::m_progress_bar_incomplete_char{' '};
    std::chrono::steady_clock::time_point PqmcWalker::m_begin_time{}, PqmcWalker::m_end_time{};

    // set up whether to show the process bar or not
    void PqmcWalker::show_progress_bar( bool show_progress_bar ) { PqmcWalker::m_show_progress_bar = show_progress_bar; }

    // set up the format of the progress bar
    void PqmcWalker::progress_bar_format( unsigned int width, char complete, char incomplete )
    {
        PqmcWalker::m_progress_bar_width = width;
        PqmcWalker::m_progress_bar_complete_char = complete;
        PqmcWalker::m_progress_bar_incomplete_char = incomplete;
    }

    // set up the rate of refreshing the progress bar
    void PqmcWalker::set_refresh_rate( unsigned int refresh_rate )
    {
        assert( refresh_rate != 0 );
        PqmcWalker::m_refresh_rate = refresh_rate;
    }

    // timer functions
    void PqmcWalker::timer_begin() { PqmcWalker::m_begin_time = std::chrono::steady_clock::now(); }
    void PqmcWalker::timer_end()   { PqmcWalker::m_end_time = std::chrono::steady_clock::now(); }
    const double PqmcWalker::timer() { 
        return std::chrono::duration_cast<std::chrono::milliseconds>(PqmcWalker::m_end_time - PqmcWalker::m_begin_time).count(); 
    }


    // ----------------------------------------------------------------------------------------------------------------------
    //
    //                       Crucial static member functions for organizing the pqmc simulations
    //
    // ----------------------------------------------------------------------------------------------------------------------

    void PqmcWalker::sweep_forth_and_back( PqmcEngine& engine, Model::Hubbard& model, Measure::MeasureHandler& meas_handler )
    {
        // sweep forth from 0 to 2theta
        engine.sweep_from_0_to_2theta( model );
        if ( meas_handler.isEqualTime() ) {
            meas_handler.equaltime_measure( engine, model );
        }

        // sweep back from 2theta to 0
        engine.sweep_from_2theta_to_0( model );
        if ( meas_handler.isEqualTime() ) {
            meas_handler.equaltime_measure( engine, model );
        }
    }
    

    void PqmcWalker::thermalize( PqmcEngine& engine, Model::Hubbard& model, const Measure::MeasureHandler& meas_handler )
    {
        if ( meas_handler.isWarmUp() ) {

            // create progress bar
            progresscpp::ProgressBar progressbar( std::ceil(meas_handler.WarmUpSweeps()/2),       // total loops 
                                                  PqmcWalker::m_progress_bar_width,               // bar width
                                                  PqmcWalker::m_progress_bar_complete_char,       // complete character
                                                  PqmcWalker::m_progress_bar_incomplete_char      // incomplete character
                                                );

            // display the progress bar
            if ( PqmcWalker::m_show_progress_bar ) {
                std::cout << " Warming up "; progressbar.display();
            }

            for ( auto sweep = 1; sweep <= std::ceil(meas_handler.WarmUpSweeps()/2); ++sweep ) {
                // sweep forth and back without measurments
                engine.sweep_from_0_to_2theta( model );
                engine.sweep_from_2theta_to_0( model );

                // record the tick
                ++progressbar;
                if ( PqmcWalker::m_show_progress_bar && ( sweep % PqmcWalker::m_refresh_rate == 0 ) ) {
                    std::cout << " Warming up "; progressbar.display();
                }
            }
                
            // progress bar finish
            if ( PqmcWalker::m_show_progress_bar ) {
                std::cout << " Warming up "; progressbar.done();
            }
        }
    }


    void PqmcWalker::measure( PqmcEngine& engine, Model::Hubbard& model, Measure::MeasureHandler& meas_handler )
    {
        if ( meas_handler.isEqualTime() || meas_handler.isDynamic() ) {
            
            // create progress bar
            const int total_ticks = meas_handler.BinNum() * ( std::ceil(meas_handler.BinCapacity()/2) + std::ceil(meas_handler.SweepsBetweenBins()/2) );
            progresscpp::ProgressBar progressbar( total_ticks,
                                                  PqmcWalker::m_progress_bar_width,
                                                  PqmcWalker::m_progress_bar_complete_char,
                                                  PqmcWalker::m_progress_bar_incomplete_char
                                                );
            
            // display the progress bar
            if ( PqmcWalker::m_show_progress_bar ) {
                std::cout << " Measuring  "; progressbar.display();
            }

            for ( auto bin = 0; bin < meas_handler.BinNum(); ++bin ) {
                // avoid correlations between adjoining bins
                for ( auto sweep = 1; sweep <= std::ceil(meas_handler.SweepsBetweenBins()/2); ++sweep ) {
                    engine.sweep_from_0_to_2theta( model );
                    engine.sweep_from_2theta_to_0( model );
                    
                    // record the tick
                    ++progressbar;
                    const int current_tick = sweep + bin * ( std::ceil(meas_handler.BinCapacity()/2) + std::ceil(meas_handler.SweepsBetweenBins()/2) );
                    if ( PqmcWalker::m_show_progress_bar && ( current_tick % PqmcWalker::m_refresh_rate == 0 ) ) {
                        std::cout << " Measuring  "; progressbar.display();
                    }
                }

                for ( auto sweep = 1; sweep <= std::ceil(meas_handler.BinCapacity()/2); ++sweep ) {
                    // update and measure
                    PqmcWalker::sweep_forth_and_back( engine, model, meas_handler );

                    // record the tick
                    ++progressbar;
                    const int current_tick = sweep + bin * std::ceil(meas_handler.BinCapacity()/2) + (bin+1) * std::ceil(meas_handler.SweepsBetweenBins()/2);
                    if ( PqmcWalker::m_show_progress_bar && ( current_tick % PqmcWalker::m_refresh_rate == 0 ) ) {
                        std::cout << " Measuring  "; progressbar.display();
                    }
                }

                // store the collected data in the MeasureHandler
                meas_handler.normalize_stats();
                meas_handler.write_stats_to_bins( bin );
                meas_handler.clear_temporary();
            }

            // progress bar finish
            if ( PqmcWalker::m_show_progress_bar ) {
                std::cout << " Measuring  "; progressbar.done();
            }
        }
    }


    void PqmcWalker::analyse( Measure::MeasureHandler& meas_handler )
    {
        // analyse the collected data after the measuring process
        meas_handler.analyse_stats();
    }


} // namespace PQMC