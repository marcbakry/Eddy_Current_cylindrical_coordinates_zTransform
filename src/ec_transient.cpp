#include "ec_transient.hpp"

// -----------
// CONSTRUCTOR
// -----------
ECTransient::ECTransient(int _nt, double _tf, std::vector<PhysicalParameters> &_pp, std::vector<SourceParameters> &_sp, std::vector<dealii::Point<2>> &_obsp, bool _verbose): m_nt(_nt), m_tf(_tf), m_physical_pars(_pp), m_source_pars(_sp), m_observation_points(_obsp), m_hs(_pp,std::vector<CDOUBLE>(_sp.size(),CDOUBLE(0.0,0.0))), m_verbose(_verbose), m_oname("../output/ref_time_ecdata") {
    // initialize time step
    if(m_verbose) std::cout << "ECTransient: Initializing solver" << std::endl;
    // 
    initialize();
}

// -------
// METHODS
// -------
void ECTransient::run() {
    if(m_verbose) std::cout << "ECTransient: Running..." << std::endl;
    // initialize
    if(!m_is_initialized) {
        std::cout << "ERROR: ECZTransient::run(): Cannot run a computation if the transient solver has not been initialized. Program will now exit..." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // 
    for(auto it=1; it<=m_nt; ++it){
        std::cout << "- " << it << " / " << m_nt << std::endl;
        auto t = it*m_dt; // current time
        auto [sol,obs] = execute_time_step(it);
        m_time_solution.push_back(sol); m_time_observable.push_back(obs);
    }
    // 
    m_is_computed = true;
    // 
    write_observables();
}

void ECTransient::display_solver_info() const {
    // 
    std::cout << "- Time span         : " << m_tf << " (s)" << std::endl;
    std::cout << "- Time step         : " << m_dt << " (s)" << std::endl;
    std::cout << "- Nb. of time steps : " << m_nt << std::endl;
    std::cout << "- Nb. of obs. points: " << m_observation_points.size() << std::endl;
    std::cout << "Warning: this solver uses internally a complex Helmholtz solver which is less efficient from the RAM perspective." << std::endl;
}

// -------
// SETTERS
// -------
void ECTransient::set_n_time_steps(int _nt) {
    m_nt = _nt;
    compute_time_step();
}

void ECTransient::set_final_time(double _tf) {
    m_tf = _tf;
    compute_time_step();
}

void ECTransient::set_output_file_name(std::string _oname) {
    m_oname = _oname;
}

// -------
// GETTERS
// -------
int ECTransient::get_n_time_steps() const {
    return m_nt;
}

double ECTransient::get_final_time() const {
    return m_tf;
}

double ECTransient::get_time_step() const {
    return m_dt;
}

// -----------
// INNER TOOLS
// -----------
std::tuple<dealii::Vector<CDOUBLE>,std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>> ECTransient::execute_time_step(int _it) {
    // case _it == 0
    if(_it <= 0) {
        return std::make_tuple(
            dealii::Vector<CDOUBLE>(m_hs.get_n_dofs()),
            std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>(m_observation_points.size(),std::make_tuple(C0,C0,C0))
        );
    }
    // first, set the new source value
    auto sp = std::vector<CDOUBLE>(m_source_pars.size());
    for(auto i=0; i<sp.size(); ++i) {
        sp[i] = CDOUBLE(m_source_pars.at(i).compute(_it*m_dt),0.0);
    }
    m_hs.set_source_parameters(sp);
    // second, set solution from previous time step
    m_hs.set_rhs_vector(m_time_solution.at(_it-1));
    // third, run and get solution
    m_hs.run();
    auto sol = m_hs.get_solution();
    // fourth, get observable
    auto obs = m_hs.compute_A_and_B_at(m_observation_points);
    // end of time step
    return std::make_tuple(sol,obs);
}

void ECTransient::initialize() {
    // 
    compute_time_step();
    m_time_solution.clear();
    m_time_observable.clear();
    // 
    m_time_solution.reserve(m_nt+1);
    m_time_observable.reserve(m_nt+1);
    // initial solution is 0
    m_time_solution.push_back(dealii::Vector<CDOUBLE>(m_hs.get_n_dofs()));
    m_time_observable.push_back(std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>(m_observation_points.size(),std::make_tuple(C0,C0,C0)));
    // 
    m_hs.set_coef(CDOUBLE(1.0/m_dt,0.0));
    // 
    m_is_initialized = true;
    m_is_computed = false;
}

void ECTransient::compute_time_step() {
    if(m_nt == 0) {
        std::cout << "ERROR: ECZTransient::compute_time_step(): invalid number of time steps (nt == 0). Program will exit." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if(m_tf <= 0.0) {
        std::cout << "ERROR: ECZTransient::compute_time_step(): invalid end-time value (tf <= 0.0). Program will exit" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    m_dt = m_tf/static_cast<double>(m_nt);
}

void ECTransient::write_observables() const {
    // The output will be written in as many files as there are observation points
    // first we open all the files with a nX_ suffix where X is the number of the node
    if(m_verbose) std::cout << "ECTransient: Writing outputs" << std::endl;
    auto nn = m_observation_points.size();
    auto files = std::vector<std::ofstream>(); files.clear();
    for(auto i=0; i<nn; ++i) {
        auto new_filename = m_oname + "_n" + std::to_string(i) + ".csv";
        files.push_back(std::ofstream());
        files.back().open(new_filename);
        if(!files.back().is_open()) {
            std::cout << "ERROR: ECZTransient::write_observables(): could not open file '" << new_filename << "'. The program will exit..." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        files.back() << "Time; A; Br; Bz" << std::endl;
    }
    // Now loop over all time steps
    auto current_time = 0.0;
    for(const auto &e: m_time_observable) { // e is a vector of tuples, one per observation node
        for(auto i=0; i<nn; ++i) { // 0: A, 1: Br, 2: Bz
            files[i] << current_time << ";" << std::get<0>(e[i]).real() << ";" << std::get<1>(e[i]).real() << ";" << std::get<2>(e[i]).real() << std::endl;
        }
        current_time += m_dt;
    }
    // Close all files
    for(auto i=0; i<nn; ++i) files[i].close();
    // Done
}