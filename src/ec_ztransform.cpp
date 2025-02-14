#include "ec_ztransform.hpp"

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
    for(const auto iz=0; iz<m_nz; ++iz) {
        m_z[iz] = m_lambda*std::exp(CDOUBLE(0.0,twopi*static_cast<double>(iz)/static_cast<double>(m_nz)));
    }
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