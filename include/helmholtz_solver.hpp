#ifndef HELMHOLTZ_SOLVER_HPP
#define HELMHOLTZ_SOLVER_HPP

#include <iostream>
#include <string>
#include <vector>

#include "defs.hpp"
#include "defs_dealii.hpp"
#include "physical_parameters.hpp"

class HelmholtzSolver {
public:
    HelmholtzSolver(const std::vector<PhysicalParameters> &_ppars, const std::vector<CDOUBLE> &_spars, const CDOUBLE _coef=CDOUBLE(0.0,0.0), const bool _print_mesh=false);

    // setters
    void set_print_mesh(bool _pm);
    void set_coef(CDOUBLE _c);
    void set_source_parameters(const std::vector<double> &_sp);

    // getters
    dealii::Vector<CDOUBLE> get_solution() const;

private: // private functions
    void print_mesh_info() const;
    void load_mesh(); // load mesh from file
    void setup_system(); // initialize mass and stiffness matrices
    void assemble_mass_and_stiffness_matrices(); // assemble the mass and stiffness matrices
    void assemble_matrix();
    void assemble_rhs();
    void assemble_system(); // assemble the system to be solved by combining the mass and stiffness matrices, the right-hand-side and applying the boundary conditions
    void solve(); // solve the system sol=lhs\rhs

private: 
    // private dealii attributes
    dealii::Triangulation<2> m_triangulation;
    dealii::FE_SimplexP<2> m_fe;
    dealii::DoFHandler<2> m_dof_handler;
    dealii::SparsityPattern m_sparsity_pattern;
    std::map<dealii::types::global_dof_index,CDOUBLE> m_boundary_values;
    dealii::SparseMatrix<CDOUBLE> m_mass_matrix, m_stiffness_matrix; // for fast rebuilding when updating the "frequency"

    dealii::SparseMatrix<CDOUBLE> m_lhs; // left-hand-side
    dealii::Vector<CDOUBLE> m_rhs; // right-hand-side
    dealii::Vector<CDOUBLE> m_sol; // solution

    // other attributes
    CDOUBLE m_coef; // coefficient
    std::vector<CDOUBLE> m_source_pars; // source parameters
    std::vector<PhysicalParameters> m_physical_pars;

    bool m_print_mesh; // print mesh for debug or not
    bool m_is_solved;
};

#endif