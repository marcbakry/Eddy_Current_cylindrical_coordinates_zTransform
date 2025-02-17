#include <iostream>

#include "eddysolver.hpp"

int main()
{
    std::cout << "+--------------------------------------+" << std::endl;
    std::cout << "|   CALCUL D'UN PROBLEME COURANTS DE   |" << std::endl;
    std::cout << "| FOUCAULT EN COORDONNEES CYLINDRIQUES |" << std::endl;
    std::cout << "+--------------------------------------+" << std::endl;
    // auto phy_pars = {PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0)};
    // auto source_pars = {0.0,0.0,1.0,-1.0};
    // try
    // {
    // EddySimulation es(phy_pars,source_pars,false);
    // es.run();
    // }
    // catch(const std::exception& e)
    // {
    //     std::cerr << e.what() << '\n';
    //     return EXIT_FAILURE;
    }
    //
    return EXIT_SUCCESS;
}