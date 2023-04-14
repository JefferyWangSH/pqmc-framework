/*
 *   observable.h
 * 
 *     Created on: Apr 14, 2023
 *         Author: Jeffery Wang
 * 
 *   This head file includes the abstract base class Observable::ObservableBase 
 *   and its template derived class Observable::Observable<ObsType> designed for measuring physical observables in QMC simualtions.
 *   Supported types of observables include scalar, vector and matrix kinds.
 */

#ifndef OBSERVABLE_H
#define OBSERVABLE_H
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <numeric>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace PQMC { class PqmcEngine; }
namespace Model { class Hubbard; }
namespace Measure{ class MeasureHandler; }

namespace Observable {

    // data types of observables
    using ScalarType = double;
    using VectorType = Eigen::VectorXd;
    using MatrixType = Eigen::MatrixXd;

    
    // -------------------------------------  Abstract base class Observable::ObservableBase  --------------------------------------
    // this class should not be instantiated in any case,
    // and it only serves as a pointer to its derived class
    class ObservableBase {
        protected:
            ObservableBase() = default;
            
        public:
            virtual ~ObservableBase(){};
    };


    
    // ----------------------------------  Derived template class Observable::Observable<ObsType>  ---------------------------------

    template< typename ObsType > class Observable : public ObservableBase {
        
        private:
            
            using MeasureMethod = void( Observable<ObsType>&, 
                                        const Measure::MeasureHandler&, const PQMC::PqmcEngine&, const Model::Hubbard& );

            ObsType m_mean_value{};                     // statistical mean value
            ObsType m_std_error{};                      // estimated standard error
            ObsType m_temp_value{};                     // temporary value used in the collections of samples
            ObsType m_zero_elem{};                      // zero element

            std::string m_name{};                       // observable name
            std::string m_desc{};                       // description of the observable, used as output message
            int m_count{0};                             // countings
            int m_bin_num{0};                           // total number of bins
            std::vector<ObsType> m_bin_data{};          // collected data in bins

            std::function<MeasureMethod> m_method{};    // user-defined measuring method

        
        public:
           
            Observable() = default;
            
            explicit Observable(int bin_num) { this->set_number_of_bins(bin_num); }

            // overload operator ++
            int operator++() { return ++this->m_count; }


            // -----------------------------------------  Interface functions  ----------------------------------------------
            
            int& counts() { return this->m_count; }
            int  counts() const { return this->m_count; }
            int  bin_num() const { return this->m_bin_num; }
            const std::string name() const { return this->m_name; }
            const std::string description() const { return this->m_desc; }
            
            const ObsType& zero_element() const { return this->m_zero_elem; } 
            const ObsType& mean_value() const { return this->m_mean_value; }
            const ObsType& std_error() const { return this->m_std_error; }
            
            const ObsType& temp_value() const { return this->m_temp_value; }
            ObsType& temp_value() { return this->m_temp_value; }

            const ObsType& bin_data( int bin ) const {
                assert( bin >= 0 && bin < this->m_bin_num );
                return this->m_bin_data[bin];
            }
            
            ObsType& bin_data( int bin ) {
                assert( bin >= 0 && bin < this->m_bin_num );
                return this->m_bin_data[bin];
            }

            std::vector<ObsType>& bin_data() { return this->m_bin_data; }
            

            // -------------------------------------  Set up parameters and methods  ----------------------------------------

            void set_number_of_bins( int bin_num ) { this->m_bin_num = bin_num; }
            void set_zero_element( const ObsType& zero_elem ) { this->m_zero_elem = zero_elem; }
            void set_name_and_description( const std::string& name, const std::string& desc ) { this->m_name = name; this->m_desc = desc; }
            void add_method( const std::function<MeasureMethod>& method ) { this->m_method = method; }


            // ----------------------------------------  Other member functions  --------------------------------------------
            
            // perform one step of measurement
            void measure( const Measure::MeasureHandler& meas_handler, 
                          const PQMC::PqmcEngine& walker,
                          const Model::Hubbard& model )
            {
                this->m_method( *this, meas_handler, walker, model ); 
            }

            // allocate memory
            void allocate() {
                this->m_mean_value = this->m_zero_elem;
                this->m_std_error = this->m_zero_elem;
                this->m_temp_value = this->m_zero_elem;

                std::vector<ObsType>().swap( this->m_bin_data );
                this->m_bin_data.reserve( this->m_bin_num );
                for ( auto i = 0; i < this->m_bin_num; ++i ) {
                    this->m_bin_data.emplace_back( this->m_zero_elem );
                }
            }

            // clear statistical data, preparing for a new measurement
            void clear_stats() {
                this->m_mean_value = this->m_zero_elem;
                this->m_std_error = this->m_zero_elem;
            }
            
            // clear temporary data
            void clear_temporary() {
                this->m_temp_value = this->m_zero_elem;
                this->m_count = 0;
            }

            // clear data of bin collections
            void clear_bin_data() {
                for ( auto& bin_data : this->m_bin_data ) {
                    bin_data = this->m_zero_elem;
                }
            }

            // perform data analysis, especially computing the mean value and error bar
            void analyse() {
                this->clear_stats();
                this->calculate_mean_value();
                this->calculate_std_error();
            }


        private:

            // calculating mean value of the measurement
            void calculate_mean_value() {
                this->m_mean_value = std::accumulate( this->m_bin_data.begin(), this->m_bin_data.end(), this->m_zero_elem );
                this->m_mean_value /= this->bin_num();
            }
            
            // estimate the error bar of the measurement
            void calculate_std_error() {
                // for observables with Scalar type
                if constexpr ( std::is_same_v< ObsType, ScalarType > ) {
                    for ( const auto& bin_data : this->m_bin_data ) {
                        this->m_std_error += std::pow( bin_data, 2 );
                    }
                    this->m_std_error /= this->bin_num();
                    this->m_std_error = std::sqrt( this->m_std_error - std::pow(this->m_mean_value, 2) ) / std::sqrt( this->bin_num()-1 );
                }

                // for observables with Vector and Matrix types
                else if constexpr ( std::is_same_v< ObsType, VectorType > || std::is_same_v< ObsType, MatrixType > ) {
                    for ( const auto& bin_data : this->m_bin_data ) {
                        this->m_std_error += bin_data.cwiseAbs2();
                    }
                    this->m_std_error /= this->bin_num();
                    this->m_std_error = ( this->m_std_error - this->m_mean_value.cwiseAbs2() ).cwiseSqrt() / std::sqrt( this->bin_num()-1 );
                }

                // others observable type, raising errors
                else {
                    std::cerr << "Observable::Observable<ObsType>::calculate_std_error(): "
                              << "undefined ObsType in Observable::Observable<ObsType>."
                              << std::endl;
                    exit(1);
                }
            }

    };


    // some aliases
    using ScalarObs = Observable<ScalarType>;
    using VectorObs = Observable<VectorType>;
    using MatrixObs = Observable<MatrixType>;


} // namespace Observable

#endif // OBSERVABLE_H