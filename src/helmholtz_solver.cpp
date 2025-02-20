#include "helmholtz_solver.hpp"

// HelmholtzSolver::HelmholtzSolver(const std::vector<PhysicalParameters> &_ppars, const std::vector<CDOUBLE> &_spars, const CDOUBLE _coef, const bool _print_mesh): m_physical_pars(_ppars), m_source_pars(_spars), m_coef(_coef), m_print_mesh(_print_mesh), m_fe(2), m_dof_handler(m_triangulation), m_is_solved(false) {
HelmholtzSolver::HelmholtzSolver(std::vector<PhysicalParameters> &_ppars, std::vector<CDOUBLE> &_spars, CDOUBLE _coef, bool _print_mesh): m_fe(2), m_dof_handler(m_triangulation), m_is_solved(false) {
    m_coef = _coef;
    m_print_mesh = _print_mesh;
    m_physical_pars = _ppars; 
    m_source_pars = _spars;

    if(m_source_pars.size() != m_physical_pars.size())
    {
        throw std::invalid_argument("\nERROR: INVALID PARAMETER DEFINITION:\n\n        Physical parameters vector and source vector should have the same length");
    }

    load_mesh();
    setup_system();
    assemble_mass_and_stiffness_matrices();
}

void HelmholtzSolver::print_mesh_info() const
{
    std::cout << "MESH INFO" << std::endl;
    std::cout << "- no. of cells: " << m_triangulation.n_active_cells() << std::endl;
    std::map<dealii::types::boundary_id,unsigned int> boundary_count;
    for(const auto &face: m_triangulation.active_face_iterators())
    {
        if(face->at_boundary()) boundary_count[face->boundary_id()]++;
    }
    std::cout << "- boundary indicators: ";
    for(const auto &pair: boundary_count)
    {
        std::cout << pair.first << " (" << pair.second << " times)" << std::endl;
    }
    std::cout << std::endl;
}

void HelmholtzSolver::load_mesh()
{
    auto mesh_file = std::string("../data/maillage_v3.vtk");
    std::cout << "- Loading mesh file :'" << mesh_file << "'. ";
    dealii::GridIn gridin(m_triangulation);
    std::ifstream file_stream(mesh_file);
    try
    {
        gridin.read_vtk(file_stream);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    file_stream.close();
    std::cout << "Done. ";

    // set boundary id to 0 for all boundary cells
    for(auto &cell: m_triangulation.active_cell_iterators())
    {
        for(auto &face: cell->face_iterators())
        {
            if(face->at_boundary())
            {
                face->set_all_boundary_ids(0);
            }
        }
    }

    if(m_print_mesh)
    {
        auto check_file = std::string("../data/mesh_check.vtk");
        std::cout << "Writing mesh for check: '" << check_file << "' ";
        std::ofstream output_stream(check_file);
        dealii::GridOut gridout;

        gridout.write_vtk(m_triangulation,output_stream);
        std::cout << "Done.";
    }
    std::cout << std::endl;
    print_mesh_info();
}

void HelmholtzSolver::set_print_mesh(bool _pm) {
    m_print_mesh = _pm;
}

void HelmholtzSolver::set_coef(CDOUBLE _c) {
    m_coef = _c;
    // since the coefficient has been changed, the pb is new
    m_is_solved = false;
}

void HelmholtzSolver::set_source_parameters(const std::vector<CDOUBLE> &_sp) {
    if(_sp.size() != m_physical_pars.size())
    {
        throw std::invalid_argument("\nERROR: INVALID PARAMETER DEFINITION:\n\n        Physical parameters vector and source vector should have the same length");
    }
    m_source_pars = _sp;
    // since the source has been changed, one must solve again
    m_is_solved = false;
}

dealii::Vector<CDOUBLE> HelmholtzSolver::get_solution() const {
    if(!m_is_solved) {
        std::cout << "ERROR: HelmholtzSolver::get_solution(): trying to access a solution which has not yet been computed. Program will stop." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return m_sol;
}

int HelmholtzSolver::get_n_dofs() const {
    return m_dof_handler.n_dofs();
}

void HelmholtzSolver::setup_system() {
    // distribute dofs - compute connectivity
    m_dof_handler.distribute_dofs(m_fe);
    // compute the sparsity pattern
    dealii::DynamicSparsityPattern dsp(m_dof_handler.n_dofs());
    dealii::DoFTools::make_sparsity_pattern(m_dof_handler,dsp);
    m_sparsity_pattern.copy_from(std::move(dsp));
    // initialize the mass matrix and stiffness matrix,
    // without computing the entries
    m_mass_matrix.reinit(m_sparsity_pattern);
    m_stiffness_matrix.reinit(m_sparsity_pattern);
    // interpolate the boundary values
    dealii::VectorTools::interpolate_boundary_values(m_dof_handler,0,dealii::Functions::ZeroFunction<2,CDOUBLE>(),m_boundary_values);
}

void HelmholtzSolver::assemble_mass_and_stiffness_matrices() {
    // initialize quadrature rule and the tool 
    // to compute the fe values (shape functions, jacobian, etc.)
    dealii::QGaussSimplex<2> quad(m_fe.degree+1);
    dealii::FEValues<2> fe_values(m_fe,quad,dealii::update_values | dealii::update_gradients | dealii::update_JxW_values | dealii::update_quadrature_points);
    // vector to store the local->global dof
    const auto dofs_per_cell = m_fe.n_dofs_per_cell();
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    // initialize the local mass matrix
    dealii::FullMatrix<CDOUBLE> cell_mass_matrix(dofs_per_cell,dofs_per_cell);
    dealii::FullMatrix<CDOUBLE> cell_stiffness_matrix(dofs_per_cell,dofs_per_cell);
    // loop over all cells
    for(const auto &cell: m_dof_handler.active_cell_iterators()) {
        // compute the useful values at the quadrature nodes
        fe_values.reinit(cell);
        cell_mass_matrix = CDOUBLE(0.0,0.0);
        cell_stiffness_matrix = CDOUBLE(0.0,0.0);

        auto sigma = m_physical_pars[cell->material_id()-1].sigma;
        auto nu    = m_physical_pars[cell->material_id()-1].nu;
        // loop over quadrature nodes
        for(const auto q: fe_values.quadrature_point_indices()) {
            auto xq = fe_values.quadrature_point(q); // pt of quadrature
            auto JxW = fe_values.JxW(q); // jacobian times quad weight
            // loop over test functions
            for(const auto i: fe_values.dof_indices()) {
                // loop over basis functions
                for(const auto j: fe_values.dof_indices()) {
                    // compute local mass matrix
                    cell_mass_matrix(i,j) += fe_values.shape_value(j,q)*fe_values.shape_value(i,q)*xq[0]*JxW*sigma;
                    // compute the local stiffness matrix
                    cell_stiffness_matrix(i,j) += (fe_values.shape_grad(i,q)[1]*fe_values.shape_grad(j,q)[1] + 1.0/(xq[0]*xq[0])*(fe_values.shape_value(i,q)+xq[0]*fe_values.shape_grad(i,q)[0])*(fe_values.shape_value(j,q)+xq[0]*fe_values.shape_grad(j,q)[0]))*xq[0]*JxW*nu;
                }
            }
        }
        // add into global matrix
        cell->get_dof_indices(local_dof_indices);
        for(const auto i: fe_values.dof_indices()) {
            for(const auto j: fe_values.dof_indices()) {
                m_mass_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_mass_matrix(i,j));
                m_stiffness_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_stiffness_matrix(i,j));
            }
        }
    }
    // 
}

void HelmholtzSolver::assemble_matrix() {
    // regroup mass and stiffness matrix using the coefficient
    m_lhs.reinit(m_sparsity_pattern);
    m_lhs.add(m_coef,m_mass_matrix); // mass matrix
    m_lhs.add(CDOUBLE(1.0,0.0),m_stiffness_matrix); // add stiffness matrix
}

void HelmholtzSolver::assemble_rhs() {
    // initialize the right-hand-side
    m_rhs.reinit(m_dof_handler.n_dofs());
    // quadrature formula
    dealii::QGaussSimplex<2> quadrature_formula(m_fe.degree+1);
    dealii::FEValues<2> fe_values(m_fe,quadrature_formula,dealii::update_values | dealii::update_JxW_values | dealii::update_quadrature_points);
    const auto dofs_per_cell = m_fe.n_dofs_per_cell();
    dealii::Vector<CDOUBLE> cell_rhs(dofs_per_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    // loop over all cells
    for(const auto &cell: m_dof_handler.active_cell_iterators()) {
        // compute values at quadrature nodes for the current cell
        fe_values.reinit(cell);
        cell_rhs = CDOUBLE(0.0,0.0);
        // loop over quadrature nodes
        for(const auto q: fe_values.quadrature_point_indices()) {
            auto xq = fe_values.quadrature_point(q);
            // loop over test functions
            for(const auto i: fe_values.dof_indices()) {
                cell_rhs(i) += m_source_pars[cell->material_id()-1]*fe_values.shape_value(i,q)*fe_values.quadrature_point(q)[0]*xq[0]*fe_values.JxW(q);
            }
        }
        // put value in the global right-hand-side
        cell->get_dof_indices(local_dof_indices);
        for(const auto i: fe_values.dof_indices()) {
            m_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
}

void HelmholtzSolver::assemble_system() {
    assemble_matrix();
    assemble_rhs();
    // apply bc
    m_sol.reinit(m_dof_handler.n_dofs());
    dealii::MatrixTools::apply_boundary_values(m_boundary_values,m_lhs,m_sol,m_rhs);
}

void HelmholtzSolver::solve() {
    // create solver based on UMFPack library
    dealii::SparseDirectUMFPACK solver;
    // factorize the matrix and discard ownership of lhs, m_lhs value is unvalidated
    solver.factorize(std::move(m_lhs));
    // solve
    solver.solve(m_rhs);
    // move result to the solution vector, m_rhs value is invalidated
    m_sol = std::move(m_rhs);
    // pb has been solver
    m_is_solved = true;
}

void HelmholtzSolver::run() {
    assemble_system();
    solve();
}

std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>> HelmholtzSolver::compte_A_and_B_at(std::vector<dealii::Point<2>> &_points) {
    // check that the solution has been computed
    if(!m_is_solved) {
        std::cout << "ERROR: HelmholtzSolver::compute_A_and_B_at(): the solution has not yet been computed. The program will exit..." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // initialize the outputs
    auto output = std::vector<std::tuple<CDOUBLE,CDOUBLE,CDOUBLE>>();
    output.clear();
    output.reserve(_points.size());
    // version 1 using  dealii tools
    for(auto ip=0; ip<_points.size(); ip++) {
        auto &p = _points.at(ip);
        // compute potential
        auto valA = dealii::VectorTools::point_value(m_dof_handler,m_sol,p);
        // compute Br and Bz
        auto sol_grad = dealii::VectorTools::point_gradient(m_dof_handler,m_sol,p);
        auto valBr = -sol_grad[1]; // Br = -\partial_z A
        auto valBz = 1.0/p[0]*(valA + p[0]*sol_grad[0]); // Bz = 1/r(A + r*\partial_r A)
        output.push_back(std::make_tuple(valA,valBr,valBz));
    }
    // maybe another version by projecting the node on the cell then evaluating
    // the local basis functions and derivatives and combining with the local
    // nodal values
    // 
    return output;
}