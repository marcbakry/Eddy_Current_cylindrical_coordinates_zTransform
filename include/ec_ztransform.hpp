#ifndef EC_ZTRANSFORM_HPP
#define EC_ZTRANSFORM_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include "defs.hpp"
#include "defs_dealii.hpp"
#include "physical_parameters.hpp"
#include "helmholtz_solver.hpp"

class ECZTransform {
public:
private:
    // getters
    int get_n_time_steps() const;
    double get_final_time() const;
    double get_time_step() const;

    int get_n_z_quadrature_nodes() const;
    double get_z_integration_radius() const;
    std::vector<CDOUBLE> get_z_quadrature_nodes() const;

    // setters
    void set_n_time_steps(int _nt);
    void set_final_time(double _tf);
    
    void set_n_z_nodes(int _nz);
    void set_integration_radius(double _lambda);

    // inner tools
    void compute_time_step();
    void compute_z_quadrature_nodes();
    std::vector<CDOUBLE> compute_source_z_transform_for_z(CDOUBLE _z) const;

private:
    // ------------------
    // private attributes
    // ------------------
    // everything model-related
    HelmholtzSolver m_hs; // the harmonic solver for each z-value
    // std::vector<double> m_source_pars; // for each subdomain, the source term
    std::vector<PhysicalParameters> m_physical_pars;

    // everything time-discretisation-related
    // time-related
    int m_nt;    // number of time steps, EXCLUDING t=0
    double m_tf; // final time in seconds
    double m_dt; // time step in seconds
    // z-transform related
    int m_nz; // number of z-quadrature nodes
    double m_lambda; // the radius of the integration contour
    std::vector<CDOUBLE> m_z; // the z-quadrature nodes
};

#endif