#include <iostream>
#include "model.h"
#include "pqmc_walker.h"
#include "random.h"

int main() {

    Model::Hubbard* hubbard = new Model::Hubbard();
    PQMC::PqmcWalker* walker = new PQMC::PqmcWalker();

    Utils::Random::set_seed( time(nullptr) );
    hubbard->initial();
    std::cout << hubbard->m_K << std::endl;
    std::cout << hubbard->m_expK << std::endl;
    
    hubbard->set_ising_fields_to_random();
    std::cout << "\n" << hubbard->m_ising_fields << std::endl;


    walker->initial( *hubbard );
    
    return 0;
}