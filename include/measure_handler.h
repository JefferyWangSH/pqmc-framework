/*
 *   measure_handler.h
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 *   This header file defines the measuring handler class Measure::MeasureHandler 
 *   to handle with the measuring process in QMC simulation,
 *   The handler class is derived from the Observable::ObservableHandler class.
 */

#ifndef MEASURE_HANDLER_H
#define MEASURE_HANDLER_H
#pragma once

#include "observable_handler.h"

namespace PQMC { class PqmcEngine; class PqmcParams; }
namespace Model { class Hubbard; }
namespace Utils { class MPI; }

namespace Measure {

    using ObsList = std::vector<std::string>;
    using MomentumIndex = int;
    using MomentumIndexList = std::vector<int>;

    // ----------------------------------------------  Measure::MeasureHandler  --------------------------------------------------
    
    class MeasureHandler : public Observable::ObservableHandler {

        private:

            bool m_is_warmup{};             // whether to warm up the system
            bool m_is_equaltime{};          // whether to perform equal-time measurements
            bool m_is_dynamic{};            // whether to perform dynamic measurements

            int m_sweeps_warmup{};          // number of the MC sweeps for the warm-up process
            int m_bin_num{};                // number of measuring bins 
            int m_bin_capacity{};           // number of samples in one measuring bin
            int m_sweeps_between_bins{};    // number of the MC sweeps between two adjoining bins

            ObsList m_observables{};        // list of observables to be measured
            
            // lattice momentum for the momentum-dependent measurements
            MomentumIndex m_momentum{};
            MomentumIndexList m_momentum_list{};

        public:

            MeasureHandler() = default;

            // ------------------------------------------  Initializations  -------------------------------------------------

            void initial( const PQMC::PqmcParams& params );
            
            
            // --------------------------------------------  Interfaces  ----------------------------------------------------

            const bool isWarmUp() const;
            const bool isEqualTime() const ;
            const bool isDynamic() const ;

            const int WarmUpSweeps() const ;
            const int SweepsBetweenBins() const;
            const int BinNum() const;
            const int BinCapacity() const;

            const MomentumIndex& Momentum() const;
            const MomentumIndex& MomentumList( int i ) const;
            const MomentumIndexList& MomentumList() const;

            
            // the following interfaces have been implemented 
            // in the base class Observable::ObservableHandler.
            // and can be directly called in the MeasureHandler class.

            // check if certain observable exists
            // bool find(const ObsName& obs_name);

            // return certain type of the observable class
            // template<typename ObsType> const ObsType find(const ObsName& obs_name);
            
            
            // -------------------------------  Subroutines for measuring the observables  ------------------------------------

            // perform one step of MC sampling for the measurements
            void equaltime_measure( const PQMC::PqmcEngine& engine, const Model::Hubbard& model );
            void dynamic_measure  ( const PQMC::PqmcEngine& engine, const Model::Hubbard& model );

            // normalize the observable samples
            void normalize_stats();
            
            // push collections of the observable samples in one bin
            void write_stats_to_bins( int bin );

            // analyse the statistics by calculating means and errors
            void analyse_stats();

            // clear the temporary data
            void clear_temporary();


            // -----------------------------------------  Friend class Utils::MPI  --------------------------------------------
            
            // for collecting the measuring data among set of MPI processes
            friend Utils::MPI;
    
    };

} // namespace Measure

#endif // MEASURE_HANDLER_H