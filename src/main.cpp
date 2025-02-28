#include <iostream>
#include <chrono>

#include "physical_parameters.hpp"
#include "source_parameters.hpp"
#include"ec_ztransform.hpp"

void run_computation(ECZTransform &_eczt, int _nt, int _nz, std::string _output_folder);

int main()
{
    std::cout << "+--------------------------------------+" << std::endl;
    std::cout << "|   CALCUL D'UN PROBLEME COURANTS DE   |" << std::endl;
    std::cout << "| FOUCAULT EN COORDONNEES CYLINDRIQUES |" << std::endl;
    std::cout << "+--------------------------------------+" << std::endl;
    // creating physical parameters
    auto sigma_air    = 0.0;
    auto sigma_piece  = 5000000.0;
    auto sigma_bobine = 0.0;

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
        SourceParameters({-1.0},{0.0}), // upper coil
    });
    // run the computation
    try
    {
        auto mm = 1e-3;
        auto ms = 1e-3;
        auto nt = 300;
        auto tf = 0.5*ms;
        auto nz = 4*nt;
        auto r  = 1.0;
        // observation nodes
        auto obsp = std::vector<dealii::Point<2>>({dealii::Point<2>(4.65*mm,-3.0*mm)});
        // initialize the solver
        auto eczt = ECZTransform(nt,tf,nz,r,phy_pars,source_pars,obsp,true,false);
        // run computations
        for(auto coef: {2}) run_computation(eczt,nt,coef*nt,"../output/");
        // for(auto coef: {1,2,4,6}) run_computation(eczt,nt,coef*nt,"../output/");
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    // the end
    return EXIT_SUCCESS;
}

// The purpose of this function is to encapsulate single computations when performing
// a parametric study on the time and "frequency" parameters: number of time steps,
// number of quadrature node.
void run_computation(ECZTransform &_eczt, int _nt, int _nz, std::string _output_folder) {
    // reinitialize solver with inputs
    _eczt.set_n_time_steps(_nt);
    _eczt.set_n_z_nodes(_nz);
    // build output name
    auto oname = _output_folder + "/fields_nt" + std::to_string(_nt) + "_nz" + std::to_string(_nz);
    _eczt.set_output_file_name(oname);
    // run 
    std::cout << "----------------------------------------------------------------" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    _eczt.run();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "COMPUTATION TIME: " << duration.count() << " (s)" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
}