#include "ec_ztransform.hpp"

// -----------
// CONSTRUCTOR
// -----------
ECZTransform::ECZTransform(int _nt, double _tf, int _nz, double _radius, std::vector<PhysicalParameters> &_pp, std::vector<SourceParameters> &_sp, std::vector<dealii::Point<2>> &_obsp, bool _verbose): m_nt(_nt), m_tf(_tf), m_nz(_nz), m_lambda(_radius), m_physical_pars(_pp), m_source_pars(_sp), m_observation_points(_obsp), m_hs(_pp,std::vector<CDOUBLE>(_sp.size(),CDOUBLE(0.0,0.0))), m_verbose(_verbose) {
    // initialize time step and z quadrature rule
    if(m_verbose) std::cout << "ECZT: Initializing solver" << std::endl;
    compute_time_step();
    compute_z_quadrature_nodes();

    m_is_time_computed = false;
    m_is_time_computed = false;
    m_is_time_observable_computed = false;
    display_solver_info();
}

// -------
// METHODS
// -------
void ECZTransform::run() {
    if(m_verbose) std::cout << "ECZT: Running..." << std::endl;
    // loop over all z
    solve_for_all_z();
    // gather results in time domain
    compute_observable_time_domain();
    // write outputs
    write_observables("../data/time_ecdata");
}

void ECZTransform::display_solver_info() const {
    // 
    std::cout << "- Time span         : " << m_tf << " (s)" << std::endl;
    std::cout << "- Time step         : " << m_dt << " (s)" << std::endl;
    std::cout << "- Nb. of time steps : " << m_nt << std::endl;
    std::cout << "- Nb. of quad nodes : " << m_nz << std::endl;
    std::cout << "- Integration radius: " << m_lambda << std::endl;
    std::cout << "- Nb. of obs. points: " << m_observation_points.size() << std::endl;
}

void ECZTransform::solve_for_all_z() {
    if(m_verbose) std::cout << "ECZT: Solving all Helmholtz problems" << std::endl;
    // clear solution vectors
    m_z_solution_symmetrical.clear();
    m_z_solution_nonsymmetrical.clear();
    m_z_solution_symmetrical.reserve(m_z_symmetrical.size());
    m_z_solution_nonsymmetrical.reserve(m_z_nonsymmetrical.size());

    m_z_observable_symmetrical.clear();
    m_z_observable_nonsymmetrical.clear();
    m_z_observable_symmetrical.reserve(m_z_symmetrical.size());
    m_z_observable_nonsymmetrical.reserve(m_z_nonsymmetrical.size());
    // compute for all symmetrical z
    if(m_verbose) std::cout << "ECZT: Solving for symmetrical z" << std::endl;
    int i=0;
    for(auto z: m_z_symmetrical) {
        if(m_verbose) std::cout << "- " << i++ << " / " << m_z_symmetrical.size() << std::endl;
        m_z_solution_symmetrical.push_back(solve_for_z(z));
        m_z_observable_symmetrical.push_back(m_hs.compte_A_and_B_at(m_observation_points));
    }
    // compute for all non symmetrical z
    if(m_verbose) std::cout << "ECZT: Solving for non-symmetrical z" << std::endl;
    i = 0;
    for(auto z: m_z_nonsymmetrical) {
        if(m_verbose) std::cout << "- " << i++ << " / " << m_z_nonsymmetrical.size() << std::endl;
        m_z_solution_nonsymmetrical.push_back(solve_for_z(z));
        m_z_observable_nonsymmetrical.push_back(m_hs.compte_A_and_B_at(m_observation_points));
    }
    // 
    m_is_z_computed = true;
}

void ECZTransform::compute_solution_time_domain() {
    if(m_verbose) std::cout << "ECZT: Compute solution in the time domain" << std::endl;
    if(!m_is_z_computed) {
        std::cout << "ERROR: ECZTransform::compute_solution_time_domain(): solutions in the z-domain must first be computed before going back to time domain. Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // now loop over all time steps
    m_time_solution.clear(); m_time_solution.reserve(m_nt+1);
    for(auto n=0; n<=m_nt; ++n) m_time_solution.push_back(compute_solution_nth_time_step(n));
    // 
    m_is_time_computed = true;
}

dealii::Vector<double> ECZTransform::compute_solution_nth_time_step(int _n) const {
    // case initial step
    if(_n == 0) {
        return dealii::Vector<double>(m_hs.get_n_dofs());
    }
    // other
    auto coef = 1.0/static_cast<double>(m_nz);
    auto res_cmplx = dealii::Vector<CDOUBLE>(m_hs.get_n_dofs());
    // add non-symmetrical
    for(auto iz=0; iz<m_z_nonsymmetrical.size(); ++iz) {
        auto zn = std::pow(m_z_nonsymmetrical.at(iz),_n);
        res_cmplx.add(1.0/zn, m_z_solution_nonsymmetrical.at(iz));
    }
    // add symmetrical
    for(auto iz=0; iz<m_z_symmetrical.size(); ++iz) {
        auto zn = std::pow(m_z_nonsymmetrical.at(iz),_n);
        // compute complex conjugate
        auto zn_bar = std::conj(zn);
        auto sol_bar = dealii::Vector<CDOUBLE>(m_hs.get_n_dofs());
        for(auto i=0; i<m_hs.get_n_dofs(); ++i) sol_bar(i) = std::conj(m_z_solution_symmetrical.at(iz)(i));
        res_cmplx.add(1.0/zn, m_z_solution_symmetrical.at(iz), 1.0/zn_bar, sol_bar);
    }
    // 
    auto res = dealii::Vector<double>(m_hs.get_n_dofs());
    for(auto i=0; i<m_hs.get_n_dofs(); ++i) {
        res(i) = res_cmplx(i).real()*coef;
    }
    // return
    return res;
}

void ECZTransform::compute_observable_time_domain() {
    if(m_verbose) std::cout << "ECZT: Compute fields in the time domain" << std::endl;
    if(!m_is_z_computed) {
        std::cout << "ERROR: ECZTransform::compute_observable_time_domain(): solutions in the z-domain must first be computed before going back to time domain. Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // now loop over all time steps
    m_time_observable.clear(); m_time_observable.reserve(m_nt+1);
    for(auto n=0; n<=m_nt; ++n) m_time_observable.push_back(compute_observable_nth_time_step(n));
    // 
    m_is_time_observable_computed = true;
}

std::vector<std::tuple<double,double,double>> ECZTransform::compute_observable_nth_time_step(int _n) const {
    auto np = m_observation_points.size();
    // case initial step
    if(_n == 0) {
        return std::vector<std::tuple<double,double,double>>(np,std::make_tuple(0.0,0.0,0.0));
    }
    // other
    auto coef = 1.0/static_cast<double>(m_nz);
    auto obs_cmplx = std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>(np);
    // add non-symmetrical
    for(auto iz=0; iz<m_z_nonsymmetrical.size(); ++iz) {
        auto zn = std::pow(m_z_nonsymmetrical.at(iz),_n);
        // loop over all nodes
        for(auto i=0; i<np; ++i) {
            auto &obs = m_z_observable_nonsymmetrical.at(iz).at(i);
            std::get<0>(obs_cmplx[i]) += std::get<0>(obs)/zn;
            std::get<1>(obs_cmplx[i]) += std::get<1>(obs)/zn;
            std::get<2>(obs_cmplx[i]) += std::get<2>(obs)/zn;
        }
    }
    // add symmetrical
    for(auto iz=0; iz<m_z_symmetrical.size(); ++iz) {
        auto zn = std::pow(m_z_symmetrical.at(iz),_n);
        // compute complex conjugate
        auto zn_bar = std::conj(zn);
        // loop over all nodes
        for(auto i=0; i<np; ++i) {
            auto &obs = m_z_observable_symmetrical.at(iz).at(i);
            // compute complex conjugate
            auto obs_bar = std::make_tuple(std::conj(std::get<0>(obs)),std::conj(std::get<1>(obs)),std::conj(std::get<2>(obs)));
            std::get<0>(obs_cmplx[i]) += std::get<0>(obs)/zn + std::get<0>(obs_bar)/zn_bar;
            std::get<1>(obs_cmplx[i]) += std::get<1>(obs)/zn + std::get<1>(obs_bar)/zn_bar;
            std::get<2>(obs_cmplx[i]) += std::get<2>(obs)/zn + std::get<2>(obs_bar)/zn_bar;
        }
    }
    // get real part
    auto obs = std::vector<std::tuple<double,double,double>>();obs.clear(); obs.reserve(np);
    for(auto i=0; i<np; ++i) {
        obs.push_back(std::make_tuple(std::get<0>(obs_cmplx.at(i)).real()*coef,std::get<1>(obs_cmplx.at(i)).real()*coef,std::get<2>(obs_cmplx.at(i)).real()*coef));
    }
    // 
    return obs;
}

dealii::Vector<CDOUBLE> ECZTransform::solve_for_z(CDOUBLE _z) {
    // initialize  z-dependent data
    std::cout << "compute source" << std::endl;
    auto sp   = compute_source_z_transform_for_z(_z);
    std::cout << "compute coef" << std::endl;
    auto coef = (CDOUBLE(1.0,0.0) - _z)/m_dt;
    // update the solver
    std::cout << "set source" << std::endl;
    m_hs.set_source_parameters(sp);
    std::cout << "set coef" << std::endl;
    m_hs.set_coef(coef);
    // compute the solution and extract it
    std::cout << "solve" << std::endl;
    m_hs.run();
    std::cout << "extract solution" << std::endl;
    return m_hs.get_solution();
}

// -----------
// INNER TOOLS
// -----------
void ECZTransform::compute_time_step() {
    if(m_nt == 0) {
        std::cout << "ERROR: ECZTransform::compute_time_step(): invalid number of time steps (nt == 0). Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if(m_tf <= 0.0) {
        std::cout << "ERROR: ECZTransform::compute_time_step(): invalid end-time value (tf <= 0.0). Program will exit" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    m_dt = m_tf/static_cast<double>(m_nt);
}

void ECZTransform::compute_z_quadrature_nodes() {
    if(m_nz == 0) {
        std::cout << "ERROR: ECZTransform::compute_z_quadrature_nodes(): invalid number of quadrature nodes (nz == 0). Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if(m_lambda <= 0.0) {
        std::cout << "ERROR: ECZTransform::compute_z_quadrature_nodes(): invalid integration radius (r <= 0.0). Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    auto twopi = 2.0*std::acos(-1.0);
    m_z = std::vector<CDOUBLE>(m_nz);
    for(auto iz=0; iz<m_nz; ++iz) {
        m_z[iz] = m_lambda*std::exp(CDOUBLE(0.0,twopi*static_cast<double>(iz)/static_cast<double>(m_nz)));
    }
    // now split into symmetrical and non-symmetrical nodes
    m_z_nonsymmetrical.clear(); m_z_symmetrical.clear();
    if(m_nz%2 == 0) { // if even
        auto mid = m_nz/2;
        m_z_nonsymmetrical.push_back(m_z.at(0)); m_z_nonsymmetrical.push_back(m_z.at(mid));
        m_z_symmetrical.reserve(mid-1);
        for(auto iz=1; iz<mid; ++iz) m_z_symmetrical.push_back(m_z.at(iz));
    } else { // if odd
        auto mid = (m_nz+1)/2;
        m_z_nonsymmetrical.push_back(m_z.at(0));
        m_z_symmetrical.reserve(mid-1);
        for(auto iz=1; iz<mid; ++iz) m_z_symmetrical.push_back(m_z.at(iz));
    }
}

std::vector<CDOUBLE> ECZTransform::compute_source_z_transform_for_z(CDOUBLE _z) {
    // 
    auto source_zt = std::vector<CDOUBLE>(m_source_pars.size());
    auto zc = CDOUBLE(1.0,0.0);
    // loop over all time steps
    for(auto it=0; it<=m_nt; ++it) {
        auto t = static_cast<double>(it)*m_dt; // current time
        // loop over all sources
        for(auto i=0; i<m_source_pars.size(); ++i) {
            source_zt[i] += m_source_pars[i].compute(t)*zc;
        }
        zc *= _z;
    }
    // 
    return source_zt;
}

void ECZTransform::reinitialize() {
    // clear vectors
    m_z_solution_nonsymmetrical.clear();
    m_z_solution_symmetrical.clear();
    m_time_solution.clear();

    m_z_observable_symmetrical.clear();
    m_z_observable_nonsymmetrical.clear();
    m_time_observable.clear();

    // reset states
    m_is_time_observable_computed = false;
    m_is_time_computed = false; 
    m_is_z_computed = false; 
}

void ECZTransform::write_observables(std::string _filename) const {
    // The output will be writtent in as many files as there are observation points
    // first we open all the files with a nX_ suffix where X is the number of the node
    if(m_verbose) std::cout << "ECZT: Writing outputs" << std::endl;
    auto nn = m_observation_points.size();
    auto files = std::vector<std::ofstream>(); files.clear();
    for(auto i=0; i<nn; ++i) {
        auto new_filename = _filename + "_n" + std::to_string(i) + ".csv";
        files.push_back(std::ofstream());
        files.back().open(new_filename);
        if(!files.back().is_open()) {
            std::cout << "ERROR: ECZTransform::write_observables(): could not open file '" << new_filename << "'. The program will exit..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        files.back() << "A; Br; Bz" << std::endl;
    }
    // Now loop over all time steps
    for(const auto &e: m_time_observable) { // e is a vector of tuples, one per observation node
        for(auto i=0; i<nn; ++i) { // 0: A, 1: Br, 2: Bz
            files[i] << std::get<0>(e[i]) << ";" << std::get<1>(e[i]) << ";" << std::get<2>(e[i]) << std::endl;
        }
    }
    // Close all files
    for(auto i=0; i<nn; ++i) files[i].close();
    // Done
}

// -------
// SETTERS
// -------
void ECZTransform::set_n_time_steps(int _nt) {
    m_nt = _nt;
    compute_time_step();
}

void ECZTransform::set_final_time(double _tf) {
    m_tf = _tf;
    compute_time_step();
}

void ECZTransform::set_n_z_nodes(int _nz) {
    m_nz = _nz;
    compute_z_quadrature_nodes();
}

void ECZTransform::set_integration_radius(double _lambda) {
    m_lambda = _lambda;
    compute_z_quadrature_nodes();
}

// -------
// GETTERS
// -------
int ECZTransform::get_n_time_steps() const {
    return m_nt;
}

double ECZTransform::get_final_time() const {
    return m_tf;
}

double ECZTransform::get_time_step() const {
    return m_dt;
}

int ECZTransform::get_n_z_quadrature_nodes() const {
    return m_nz;
}

double ECZTransform::get_z_integration_radius() const {
    return m_lambda;
}

std::vector<CDOUBLE> ECZTransform::get_z_quadrature_nodes() const {
    return m_z;
}