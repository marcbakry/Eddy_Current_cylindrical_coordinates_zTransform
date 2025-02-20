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
    auto phy_pars = std::vector<PhysicalParameters>({PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0),PhysicalParameters(1000.0,1.0)});
    // auto source_pars = {0.0,0.0,1.0,-1.0};
    auto source_pars = std::vector<SourceParameters>({
        SourceParameters(std::vector<double>({0.0}),std::vector<double>({0.0})),
        SourceParameters(std::vector<double>({0.0}),std::vector<double>({0.0})),
        SourceParameters(std::vector<double>({1.0}),std::vector<double>({0.0})),
        SourceParameters(std::vector<double>({-1.0}),std::vector<double>({0.0})),
    });
    try
    {
        auto nt = 20;
        auto tf = 20.0*1e-3;
        auto nz = 2*nt;
        auto r  = 1.0;
        // creation des points d'observation
        auto obsp = std::vector<dealii::Point<2>>();
        // creation du solveur temporel
        auto eczt = ECZTransform(nt,tf,nz,r,phy_pars,source_pars,obsp);
        // run computations
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    //
    return EXIT_SUCCESS;
}