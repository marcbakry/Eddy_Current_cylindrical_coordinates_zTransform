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
    auto sigma_air    = 0.0;
    auto sigma_piece  = 5000000.0;
    auto sigma_bobine = 50000000.0;

    auto mu0 = 4.0*std::acos(-1.0)*1e-7;
    auto mur = 1.0;

    auto mu_piece  = mu0*mur;
    auto mu_air    = mu0;
    auto mu_bobine = mu0;

    auto phy_pars = std::vector<PhysicalParameters>({
        PhysicalParameters(sigma_piece,1.0/mu_piece),
        PhysicalParameters(sigma_air,1.0/mu_air),
        PhysicalParameters(sigma_air,1.0/mu_air),
        PhysicalParameters(sigma_bobine,1.0/mu_bobine)
    });
    // source intensities
    auto source_pars = std::vector<SourceParameters>({
        SourceParameters({0.0},{0.0}), // plate
        SourceParameters({0.0},{0.0}), // air
        SourceParameters({0.0},{0.0}), // lower coil == air in our configuration
        SourceParameters({1.0},{0.0}), // upper coil
    });
    // run the computation
    try
    {
        auto mm = 1e-3;
        auto ms = 1e-3;
        auto nt = 200;
        auto tf = 0.5*ms;
        auto nz = 2*nt;
        auto r  = 1.0;
        // observation nodes
        auto obsp = std::vector<dealii::Point<2>>({dealii::Point<2>(4.65*mm,-3.0*mm)});
        // initialize the solver
        auto eczt = ECZTransform(nt,tf,nz,r,phy_pars,source_pars,obsp,true,false);
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