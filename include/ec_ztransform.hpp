#ifndef EC_ZTRANSFORM_HPP
#define EC_ZTRANSFORM_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "defs.hpp"
#include "defs_dealii.hpp"
#include "source_parameters.hpp"
#include "physical_parameters.hpp"
#include "helmholtz_solver.hpp"

class ECZTransform {
public:
    ECZTransform(int _nt, double _tf, int _nz, double _radius, std::vector<PhysicalParameters> &_pp, std::vector<SourceParameters> &_sp, std::vector<dealii::Point<2>> &_obsp = std::vector<dealii::Point<2>>(), bool _verbose=false, bool _debug_info=false);

    // 
    void run();
    void display_solver_info() const;
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
    std::vector<CDOUBLE> compute_source_z_transform_for_z(CDOUBLE _z);
    dealii::Vector<CDOUBLE> solve_for_z(CDOUBLE _z);
    void solve_for_all_z();

    dealii::Vector<double> compute_solution_nth_time_step(int _n) const; // compute the n-th time step
    void compute_solution_time_domain();

    std::vector<std::tuple<double,double,double>> compute_observable_nth_time_step(int _n) const; // compute the observables at the n-th time step
    void compute_observable_time_domain();

    void write_observables(std::string _filename) const; // write fields in *.csv format with double precision

    void reinitialize(); // reset all inner data, for example when setting a new time value

    // debug functions
    void write_quadrature_nodes() const;

private:
    // ------------------
    // private attributes
    // ------------------
    // everything model-related
    HelmholtzSolver m_hs; // the harmonic solver for each z-value
    std::vector<SourceParameters> m_source_pars; // for each subdomain, the source term
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
    std::vector<CDOUBLE> m_z_symmetrical; // symmetrical nodes with respect to the real axis
    std::vector<CDOUBLE> m_z_nonsymmetrical; // non symmetrical nodes

    // solutions in z- and time-domain
    std::vector<dealii::Vector<double>> m_time_solution;
    std::vector<dealii::Vector<CDOUBLE>> m_z_solution_symmetrical; // see above
    std::vector<dealii::Vector<CDOUBLE>> m_z_solution_nonsymmetrical;
    // observables in z- and time-domain
    std::vector<dealii::Point<2>> m_observation_points;
    std::vector<std::vector<std::tuple<double,double,double>>> m_time_observable;
    std::vector<std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>> m_z_observable_symmetrical;
    std::vector<std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>> m_z_observable_nonsymmetrical;


    bool m_verbose;
    bool m_debug_info;
    bool m_is_z_computed;
    bool m_is_time_computed;
    bool m_is_time_observable_computed;
};

#endif