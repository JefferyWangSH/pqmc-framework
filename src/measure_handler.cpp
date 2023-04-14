/*
 *   measure_handler.cpp
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "measure_handler.h"
#include "pqmc_engine.h"
#include "model.h"
#include "pqmc_params.hpp"

namespace Measure {

    const bool MeasureHandler::isWarmUp() const { return this->m_is_warmup; }
    const bool MeasureHandler::isEqualTime() const { return this->m_is_equaltime; }
    const bool MeasureHandler::isDynamic() const { return this->m_is_dynamic; }

    const int MeasureHandler::WarmUpSweeps() const { return this->m_sweeps_warmup; }
    const int MeasureHandler::SweepsBetweenBins() const { return this->m_sweeps_between_bins; }
    const int MeasureHandler::BinNum() const { return this->m_bin_num; }
    const int MeasureHandler::BinCapacity() const { return this->m_bin_capacity; }

    const MomentumIndex& MeasureHandler::Momentum() const { return this->m_momentum; }
    const MomentumIndexList& MeasureHandler::MomentumList() const { return this->m_momentum_list; }
    const MomentumIndex& MeasureHandler::MomentumList( int i ) const 
    { 
        assert( i >= 0 && i < (int)this->m_momentum_list.size() );
        return this->m_momentum_list[i]; 
    }


    void MeasureHandler::initial( const PQMC::PqmcParams& params )
    {
        assert( params.sweeps_warmup >= 0 );
        assert( params.bin_num >= 0 );
        assert( params.bin_capacity >= 0 );
        assert( params.sweeps_between_bins >= 0 );

        this->m_sweeps_warmup = params.sweeps_warmup;
        this->m_bin_num = params.bin_num;
        this->m_bin_capacity = params.bin_capacity;
        this->m_sweeps_between_bins = params.sweeps_between_bins;
        
        this->m_obs_list = params.obs_list;
        this->m_momentum = params.momentum;
        this->m_momentum_list = params.momentum_list;

        // initialize ObservableHandler
        Observable::ObservableHandler::initial( params );

        this->m_is_warmup = ( this->m_sweeps_warmup != 0 );
        this->m_is_equaltime = ( !this->m_eqtime_scalar_obs.empty() || !this->m_eqtime_vector_obs.empty() || !this->m_eqtime_matrix_obs.empty() );
        this->m_is_dynamic = ( !this->m_dynamic_scalar_obs.empty() || !this->m_dynamic_vector_obs.empty() || !this->m_dynamic_matrix_obs.empty() );
    }

    
    void MeasureHandler::equaltime_measure( const PQMC::PqmcEngine& engine, const Model::Hubbard& model )
    {
        for ( auto& scalar_obs : this->m_eqtime_scalar_obs ) { scalar_obs->measure( *this, engine, model ); }
        for ( auto& vector_obs : this->m_eqtime_vector_obs ) { vector_obs->measure( *this, engine, model ); }
        for ( auto& matrix_obs : this->m_eqtime_matrix_obs ) { matrix_obs->measure( *this, engine, model ); }
    }

    
    void MeasureHandler::dynamic_measure( const PQMC::PqmcEngine& engine, const Model::Hubbard& model )
    {
        for ( auto& scalar_obs : this->m_dynamic_scalar_obs ) { scalar_obs->measure( *this, engine, model ); }
        for ( auto& vector_obs : this->m_dynamic_vector_obs ) { vector_obs->measure( *this, engine, model ); }
        for ( auto& matrix_obs : this->m_dynamic_matrix_obs ) { matrix_obs->measure( *this, engine, model ); }
    }

    
    void MeasureHandler::normalize_stats()
    {   
        if ( this->m_is_equaltime ) {
            // normalize observables by the countings
            for ( auto& scalar_obs : this->m_eqtime_scalar_obs ) { scalar_obs->temp_value() /= scalar_obs->counts(); }
            for ( auto& vector_obs : this->m_eqtime_vector_obs ) { vector_obs->temp_value() /= vector_obs->counts(); }
            for ( auto& matrix_obs : this->m_eqtime_matrix_obs ) { matrix_obs->temp_value() /= matrix_obs->counts(); }
        }

        if ( this->m_is_dynamic ) {
            for ( auto& scalar_obs : this->m_dynamic_scalar_obs ) { scalar_obs->temp_value() /= scalar_obs->counts(); }
            for ( auto& vector_obs : this->m_dynamic_vector_obs ) { vector_obs->temp_value() /= vector_obs->counts(); }
            for ( auto& matrix_obs : this->m_dynamic_matrix_obs ) { matrix_obs->temp_value() /= matrix_obs->counts(); }
        }
    }


    void MeasureHandler::write_stats_to_bins( int bin )
    {
        if ( this->m_is_equaltime ) {
            for ( auto& scalar_obs : this->m_eqtime_scalar_obs ) { scalar_obs->bin_data(bin) = scalar_obs->temp_value(); }
            for ( auto& vector_obs : this->m_eqtime_vector_obs ) { vector_obs->bin_data(bin) = vector_obs->temp_value(); }
            for ( auto& matrix_obs : this->m_eqtime_matrix_obs ) { matrix_obs->bin_data(bin) = matrix_obs->temp_value(); }
        }

        if ( this->m_is_dynamic ) {
            for ( auto& scalar_obs : this->m_dynamic_scalar_obs ) { scalar_obs->bin_data(bin) = scalar_obs->temp_value(); }
            for ( auto& vector_obs : this->m_dynamic_vector_obs ) { vector_obs->bin_data(bin) = vector_obs->temp_value(); }
            for ( auto& matrix_obs : this->m_dynamic_matrix_obs ) { matrix_obs->bin_data(bin) = matrix_obs->temp_value(); }
        }
    }


    void MeasureHandler::analyse_stats()
    {
        if ( this->m_is_equaltime ) {
            for ( auto& scalar_obs : this->m_eqtime_scalar_obs ) { scalar_obs->analyse(); }
            for ( auto& vector_obs : this->m_eqtime_vector_obs ) { vector_obs->analyse(); }
            for ( auto& matrix_obs : this->m_eqtime_matrix_obs ) { matrix_obs->analyse(); }
        }

        if ( this->m_is_dynamic ) {
            for ( auto& scalar_obs : this->m_dynamic_scalar_obs ) { scalar_obs->analyse(); }
            for ( auto& vector_obs : this->m_dynamic_vector_obs ) { vector_obs->analyse(); }
            for ( auto& matrix_obs : this->m_dynamic_matrix_obs ) { matrix_obs->analyse(); }
        }
    }


    void MeasureHandler::clear_temporary()
    {   
        // clear temporary data for all observables
        if ( this->m_is_equaltime ) {
            for ( auto& scalar_obs : this->m_eqtime_scalar_obs ) { scalar_obs->clear_temporary(); }
            for ( auto& vector_obs : this->m_eqtime_vector_obs ) { vector_obs->clear_temporary(); }
            for ( auto& matrix_obs : this->m_eqtime_matrix_obs ) { matrix_obs->clear_temporary(); }
        }

        if ( this->m_is_dynamic ) {
            for ( auto& scalar_obs : this->m_dynamic_scalar_obs ) { scalar_obs->clear_temporary(); }
            for ( auto& vector_obs : this->m_dynamic_vector_obs ) { vector_obs->clear_temporary(); }
            for ( auto& matrix_obs : this->m_dynamic_matrix_obs ) { matrix_obs->clear_temporary(); }
        }
    }


} // namespace Measure