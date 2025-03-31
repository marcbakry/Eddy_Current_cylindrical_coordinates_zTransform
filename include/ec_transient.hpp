#ifndef EC_TRANSIENT_HPP
#define EC_TRANSIENT_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "defs.hpp"
#include "defs_dealii.hpp"
#include "source_parameters.hpp"
#include "physical_parameters.hpp"
#include "helmholtz_solver.hpp"

class ECTransient {
public: // public functions
    ECTransient(int _nt, double _tf, std::vector<PhysicalParameters> &_pp, std::vector<SourceParameters> &_sp, std::vector<dealii::Point<2>> &_obsp = std::vector<dealii::Point<2>>(), bool _verbose=false);

    // 
    void run();
    void display_solver_info() const;

    // getters
    int get_n_time_steps() const;
    double get_final_time() const;
    double get_time_step() const;

    // setters
    void set_n_time_steps(int _nt);
    void set_final_time(double _tf);
    
    void set_output_file_name(std::string _oname);

private: // private functions
    // inner tools
    void initialize(); // initialize everything
    void compute_time_step(); // compute dt
    std::tuple<dealii::Vector<CDOUBLE>,std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>> execute_time_step(int _it); // compute the current solution
    void write_observables() const;
private: // private attributes
    // time discretisation-related
    int m_nt; // number of time steps
    double m_tf; // final time
    double m_dt; // time step

    // model- and solver-related
    HelmholtzSolver m_hs; // the harmonic solver for each z-value
    std::vector<SourceParameters> m_source_pars; // for each subdomain, the source term
    std::vector<PhysicalParameters> m_physical_pars;

    std::string m_oname; // output file

    // results in time domain
    std::vector<dealii::Point<2>> m_observation_points;
    std::vector<dealii::Vector<CDOUBLE>> m_time_solution;
    std::vector<std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>> m_time_observable;

    bool m_verbose;
    bool m_is_initialized;
    bool m_is_computed;
};

#endif