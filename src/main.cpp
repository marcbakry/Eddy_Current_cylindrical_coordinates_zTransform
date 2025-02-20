#include <iostream>

#include "physical_parameters.hpp"
#include "source_parameters.hpp"
#include"ec_ztransform.hpp"

int main()
{
    std::cout << "+--------------------------------------+" << std::endl;
    std::cout << "|   CALCUL D'UN PROBLEME COURANTS DE   |" << std::endl;
    std::cout << "| FOUCAULT EN COORDONNEES CYLINDRIQUES |" << std::endl;
    std::cout << "+--------------------------------------+" << std::endl;
    // creating physical parameters
    auto phy_pars = std::vector<PhysicalParameters>({PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0)});
    // source intensities
    auto source_pars = std::vector<SourceParameters>({
        SourceParameters({0.0},{0.0}),
        SourceParameters({0.0},{0.0}),
        SourceParameters({1.0},{0.0}),
        SourceParameters({-1.0},{0.0}),
    });
    // run the computation
    try
    {
        auto nt = 200;
        auto tf = 20.0*1e-3;
        auto nz = 2*nt;
        auto r  = 1.0;
        // observation nodes
        auto obsp = std::vector<dealii::Point<2>>({dealii::Point<2>(0.1,0.0)});
        // initialize the solver
        auto eczt = ECZTransform(nt,tf,nz,r,phy_pars,source_pars,obsp,true);
        // run computations
        eczt.run();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    // the end
    return EXIT_SUCCESS;
}