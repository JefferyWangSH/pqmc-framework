/*
 *   pqmc_walker.h
 * 
 *     Created on: Apr 13, 2023
 *         Author: Jeffery Wang
 * 
 */

#ifndef PQMC_WALKER_H
#define PQMC_WALKER_H
#pragma once

#include <chrono>

namespace Model { class Hubbard; }
namespace Measure { class MeasureHandler; }

namespace PQMC {

    class PqmcEngine;

    class PqmcWalker {
        
        public:

            static void show_progress_bar( bool show_progress_bar );
            static void progress_bar_format( unsigned int width, char complete, char incomplete );
            static void set_refresh_rate( unsigned int refresh_rate );

            static const double timer();
            static void timer_begin();
            static void timer_end();

            static void thermalize( PqmcEngine& engine, Model::Hubbard& model, const Measure::MeasureHandler& meas_handler );
            static void measure( PqmcEngine& engine, Model::Hubbard& model, Measure::MeasureHandler& meas_handler );
            static void analyse( Measure::MeasureHandler& meas_handler );
        
        private:

            static bool m_show_progress_bar;
            static unsigned int m_progress_bar_width;
            static unsigned int m_refresh_rate;
            static char m_progress_bar_complete_char, m_progress_bar_incomplete_char;

            static std::chrono::steady_clock::time_point m_begin_time, m_end_time;

            static void sweep_forth_and_back( PqmcEngine& engine, Model::Hubbard& model, Measure::MeasureHandler& meas_handler );

    };

} // namespace Pqmc

#endif // PQMC_WALKER_H