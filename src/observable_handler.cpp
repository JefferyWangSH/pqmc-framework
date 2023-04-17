/*
 *   observable_handler.cpp
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 */

#include "observable_handler.h"
#include "measure_methods.hpp"
#include "pqmc_params.hpp"
#include <iostream>

namespace Observable {

    // definitions of static members
    // tables of all supported physical observables for measurements
    // equal-time observables
    ObsTable ObservableHandler::m_eqtime_obs_table  = {
                                    "DoubleOccupation",
                                    "KineticEnergy",
                                };

    // dynamic observables
    ObsTable ObservableHandler::m_dynamic_obs_table = {
                                    "GreensFunctions",
                                    "DensityOfStates",
                                    "ProjectionBenchmark",
                                };

    // public member for external calls
    ObsTable ObservableHandler::ObservableAll = {
                                    "DoubleOccupation",
                                    "KineticEnergy",
                                    "GreensFunctions",
                                    "DensityOfStates",
                                    "ProjectionBenchmark",
                                };


    bool ObservableHandler::is_eqtime( const ObsName& obs_name ) const 
    {
        return ( std::find( this->m_eqtime_obs_table.begin(), 
                            this->m_eqtime_obs_table.end(), 
                            obs_name  ) 
                != this->m_eqtime_obs_table.end() );
    }


    bool ObservableHandler::is_dynamic( const ObsName& obs_name ) const 
    {
        return ( std::find( this->m_dynamic_obs_table.begin(), 
                            this->m_dynamic_obs_table.end(), 
                            obs_name  ) 
                != this->m_dynamic_obs_table.end() );
    }


    bool ObservableHandler::find( const ObsName& obs_name ) 
    {
        return ( this->m_obs_map.find(obs_name) != this->m_obs_map.end() );
    }


    bool ObservableHandler::check_validity( const ObsNameList& obs_list ) const 
    {
        // preprocessing
        // remove redundant input
        ObsNameList tmp_list = obs_list;
        std::sort( tmp_list.begin(), tmp_list.end() );
        tmp_list.erase( unique(tmp_list.begin(), tmp_list.end() ), tmp_list.end() );

        // check the validity of the input
        for ( const auto& obs : tmp_list ) {
            if ( !this->is_eqtime(obs) && !this->is_dynamic(obs) ) {
                return false;
            }
        }
        // otherwise return true
        return true;
    }


    void ObservableHandler::deallocate() 
    {
        this->m_obs_map.clear();

        this->m_eqtime_scalar_obs.clear();
        this->m_eqtime_vector_obs.clear();
        this->m_eqtime_matrix_obs.clear();
        this->m_dynamic_scalar_obs.clear();
        this->m_dynamic_vector_obs.clear();
        this->m_dynamic_matrix_obs.clear();

        this->m_eqtime_scalar_obs.shrink_to_fit();
        this->m_eqtime_vector_obs.shrink_to_fit();
        this->m_eqtime_matrix_obs.shrink_to_fit();
        this->m_dynamic_scalar_obs.shrink_to_fit();
        this->m_dynamic_vector_obs.shrink_to_fit();
        this->m_dynamic_matrix_obs.shrink_to_fit();
    }


    void ObservableHandler::initial( const PQMC::PqmcParams& params )
    {
        // release memory if previously initialized
        this->deallocate();

        // check the validity of the input
        if ( !this->check_validity( params.observables ) ) {
            // unsupported observables found, throw errors
            std::cerr << "Observable::ObservableHandler::initial(): "
                      << "unsupported observable type from the input." << std::endl;
            exit(1);
        }

        for ( const auto& obs_name : params.observables ) {
            // allocate memory for observables
            // caution that to this stage only properties like name or method is assigned.
            // other info, such as m_size_of_bin and the dimension of m_zero_elem,
            // is kept unassigned until the MeasureHandler class is specialized.

            // -------------------------------------------------------------------------------------------------------------
            //                                          Equal-time Observables
            // -------------------------------------------------------------------------------------------------------------

            // ----------------------------------------  Double Occupation D -----------------------------------------------
            if ( obs_name == "DoubleOccupation" ) {
                ptrScalarObs double_occupation = std::make_shared<ScalarObs>();
                double_occupation->set_name_and_description( obs_name, "Double Occupation" );
                double_occupation->add_method( Measure::Methods::measure_double_occupation );
                double_occupation->set_zero_element( 0.0 );
                double_occupation->set_number_of_bins( params.bin_num );
                double_occupation->allocate();
                this->m_eqtime_scalar_obs.emplace_back( double_occupation );
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>( double_occupation );
            }
            
            // -----------------------------------------  Kinetic Energy Ek ------------------------------------------------
            if ( obs_name == "KineticEnergy" ) {
                ptrScalarObs kinetic_energy = std::make_shared<ScalarObs>();
                kinetic_energy->set_name_and_description( obs_name, "Kinetic Energy" );
                kinetic_energy->add_method( Measure::Methods::measure_kinetic_energy );
                kinetic_energy->set_zero_element( 0.0 );
                kinetic_energy->set_number_of_bins( params.bin_num );
                kinetic_energy->allocate();
                this->m_eqtime_scalar_obs.emplace_back( kinetic_energy );
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>( kinetic_energy );
            }

            // adding new methods here


            // -------------------------------------------------------------------------------------------------------------
            //                                        Time-displaced Observables
            // -------------------------------------------------------------------------------------------------------------

            // --------------------------------------  Green's functions G(k,t)  -------------------------------------------
            if ( obs_name == "GreensFunctions" ) {
                ptrMatrixObs greens_functions = std::make_shared<MatrixObs>();
                greens_functions->set_name_and_description( obs_name, "Green's functions" );
                greens_functions->add_method( Measure::Methods::measure_greens_functions );
                greens_functions->set_zero_element( MatrixType::Zero(params.momentum_list.size(), params.ntm) );
                greens_functions->set_number_of_bins( params.bin_num );
                greens_functions->allocate();
                this->m_dynamic_matrix_obs.emplace_back( greens_functions );
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>( greens_functions );
            }

            // ----------------------------------------  Density of states D(t)  -------------------------------------------
            if ( obs_name == "DensityOfStates" ) {
                ptrMatrixObs density_of_states = std::make_shared<MatrixObs>();
                density_of_states->set_name_and_description( obs_name, "Density of states" );
                density_of_states->add_method( Measure::Methods::measure_density_of_states );
                density_of_states->set_zero_element( MatrixType::Zero(2, params.ntm) );     // two for su(2) fermion
                density_of_states->set_number_of_bins( params.bin_num );
                density_of_states->allocate();
                this->m_dynamic_matrix_obs.emplace_back( density_of_states );
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>( density_of_states );
            }

            // -------------------------------------  Projection length benchmark  -----------------------------------------
            if ( obs_name == "ProjectionBenchmark" ) {
                ptrVectorObs projection_benchmark = std::make_shared<VectorObs>();
                projection_benchmark->set_name_and_description( obs_name, "Projection benchmark" );
                projection_benchmark->add_method( Measure::Methods::projection_benchamrk_measure );
                projection_benchmark->set_zero_element( VectorType::Zero( params.ntm ) );
                projection_benchmark->set_number_of_bins( params.bin_num );
                projection_benchmark->allocate();
                this->m_dynamic_vector_obs.emplace_back( projection_benchmark );
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>( projection_benchmark );
            }

            // add new methods here
        }
    }

} // namespace Observable