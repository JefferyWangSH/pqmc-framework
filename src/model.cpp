#include "model.h"
#include "random.h"
#include <unsupported/Eigen/MatrixFunctions>

namespace Model {

    void Hubbard::initial()
    {
        this->m_ns = this->m_nl * this->m_nl;
        this->m_alpha = std::acosh( std::exp(0.5 * this->m_dt * this->m_u) );
        this->m_ising_fields.resize( this->m_nt, this->m_nl*this->m_nl );

        this->m_K.resize( this->m_ns, this->m_ns );
        for ( int x = 0; x < this->m_nl; ++x ) {
            for ( int y = 0; y < this->m_nl; ++y ) {
                const int i = x + y * this->m_nl;
                const int ixplus1 = (x+1)%this->m_nl + y * this->m_nl;
                const int iyplus1 = x + (y+1)%this->m_nl * this->m_nl;
                this->m_K( i, ixplus1 ) += - this->m_t;
                this->m_K( ixplus1, i ) += - this->m_t;
                this->m_K( i, iyplus1 ) += - this->m_t;
                this->m_K( iyplus1, i ) += - this->m_t;
            }
        }
        this->m_expK     = ( -this->m_dt * this->m_K ).exp();
        this->m_inv_expK = ( +this->m_dt * this->m_K ).exp();
    }

    void Hubbard::set_ising_fields_to_random()
    {
        std::bernoulli_distribution bernoulli_dist(0.5);
        for ( auto t = 0; t < this->m_nt; ++t ) {
            for ( auto i = 0; i < this->m_ns; ++i ) {
                this->m_ising_fields(t, i) = bernoulli_dist( Utils::Random::Engine )? +1.0 : -1.0;
            }
        }
    }

    void Hubbard::update_ising_field( timeIndex t, spaceIndex i )
    {
        this->m_ising_fields(t, i) = - this->m_ising_fields(t, i);
    }

    void Hubbard::update_greens_function( timeIndex t, spaceIndex i, spinIndex s )
    {

    }

    const double Hubbard::get_update_ratio ( timeIndex t, spaceIndex i, spinIndex s ) const
    {

    }

    void Hubbard::multiply_B_from_left( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        green = this->m_expK * green;
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.row(i) *= std::exp( +s * this->m_alpha * this->m_ising_fields(t, i) );
        }
    }

    void Hubbard::multiply_B_from_right( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.col(i) *= std::exp( +s * this->m_alpha * this->m_ising_fields(t, i) );
        }
        green = green * this->m_expK;
    }

    void Hubbard::multiply_invB_from_left( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.row(i) *= std::exp( -s * this->m_alpha * this->m_ising_fields(t, i) );
        }
        green = this->m_inv_expK * green;
    }

    void Hubbard::multiply_invB_from_right( GreensFunction& green, timeIndex t, spinIndex s ) const
    {
        green = green * this->m_inv_expK;
        for ( auto i = 0; i < this->m_ns; ++i ) {
            green.col(i) *= std::exp( -s * this->m_alpha * this->m_ising_fields(t, i) );
        }
    }

} // namespace Model