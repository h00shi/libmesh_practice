// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fstream>


// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"


// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Exact solution
#include "libmesh/exact_solution.h"

// Get pot
#include "libmesh/getpot.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// some global variables
bool b_sym_sol;
bool b_exact_bdry;

// straighten the boundary of a mesh
void mesh_straighten(MeshBase& );
// assemble the system
void assemble_poisson(EquationSystems& es, const std::string& system_name);
void find_error(EquationSystems& es);

Real exact_solution (const Real x, const Real y, const Real z = 0.){
	const double r2 = x*x + y*y;
	const double r = sqrt(r2);
	const double cs = x/r;

	if (b_sym_sol) return -r2/4. + 0.25;
	else return (-r2/4.) * exp(cs);
}
Real right_hand_side(const Real x, const Real y, const Real z = 0.){
	const double r2 = x*x + y*y;
	const double r = sqrt(r2);
	const double cs = x/r;
	const double sn = y/r;

	if (b_sym_sol) return 1;
	else return exp(cs) * ( (sn*sn-cs)/4. + 1);
}
Real boundary_value(const Real x, const Real y, const Real z = 0.){
	const double r2 = x*x + y*y;
	const double r = sqrt(r2);
	const double cs = x/r;

	if (b_exact_bdry) return exact_solution(x,y,z);
	else
	{
		if (b_sym_sol) return 0;
		else return -exp(cs)/4.;
	}
}
Number exact_solution_libmesh(const Point& p,
                      const Parameters&,  // parameters, not needed
                      const std::string&, // sys_name, not needed
                      const std::string&) // unk_name, not needed
{
  return exact_solution (p(0), p(1));
}

Gradient exact_deriv_libmesh(const Point& p,
                      const Parameters&,  // parameters, not needed
                      const std::string&, // sys_name, not needed
                      const std::string&) // unk_name, not needed
{
	Gradient grad;
	const double r = sqrt(p(0)*p(0) + p(1)*p(1));
	const double cs = p(0)/r;
	const double sn = p(1)/r;

	if (b_sym_sol){
		grad(0) = -p(0)/2;
		grad(1) = -p(1)/2;
	} else{
		const double dudr = -r/2. * exp(cs);
		const double dudt = r/4. * sn * exp(cs);
		grad(0) = dudr * cs - dudt * sn;
		grad(1) = dudr * sn + dudt * cs;
	}

	return grad;
}


int main (int argc, char** argv)
{
	// Initialize libraries, like in example 2.
	LibMeshInit init (argc, argv);
	GetPot getpot(argc, argv);
	b_exact_bdry= getpot("exact_bdry", false);
	b_sym_sol = getpot("sym_sol",true);
	bool b_straight = getpot("l",false);
	std::string s_out_name = getpot("out", "loglog.txt");


	std::cout << "Running " << argv[0];
	for (int i=1; i<argc; i++)
		std::cout << " " << argv[i];
	std::cout << std::endl << std::endl;

	// Skip this 2D example if libMesh was compiled as 1D-only.
	if(argc < 2) libmesh_error_msg("Mesh file not given.");

	// create a mesh
	Mesh mesh(init.comm());
	GmshIO gmsh(mesh);
	{
		char gmsh_name[200];
		sprintf(gmsh_name, "mesh/%s.msh", argv[1]);
		gmsh.read(gmsh_name);
	}
	mesh.prepare_for_use();
	if(b_straight){
		std::cout << "--- Straightening mesh ---" << std::endl;
		mesh_straighten(mesh);
	}

	// print the mesh before editing
	mesh.print_info();
	{
		char gmsh_name[200];
		sprintf(gmsh_name, "out/%s_curved.msh", argv[1]);
		gmsh.write(gmsh_name);
	}

	// define the system for u
	EquationSystems equation_systems (mesh);
	equation_systems.add_system<LinearImplicitSystem> ("Poisson");
	equation_systems.get_system("Poisson").add_variable("u", SECOND, LAGRANGE);
	equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
	// find the error
	equation_systems.add_system<ExplicitSystem> ("Error");
	equation_systems.get_system("Error").add_variable("e", CONSTANT, MONOMIAL);
	// solve and find the error
	equation_systems.init();
	equation_systems.print_info();
	equation_systems.get_system("Poisson").solve();
	find_error(equation_systems);


	// find the exact solution
	ExactSolution exact_sol(equation_systems);
	exact_sol.attach_exact_value(exact_solution_libmesh);
	exact_sol.attach_exact_deriv(exact_deriv_libmesh);
	exact_sol.extra_quadrature_order(1);
	exact_sol.compute_error("Poisson", "u");

	std::fstream file_out;
	file_out.open(("out/" + s_out_name).c_str(), std::fstream::app | std::fstream::out);
	file_out <<  exact_sol.l2_error("Poisson", "u") << " "
			<< exact_sol.h1_error("Poisson", "u") << std::endl;
	std::cout <<  exact_sol.l2_error("Poisson", "u") << " "
			<< exact_sol.h1_error("Poisson", "u") << std::endl;
	file_out.close();

	// After solving the system write the solution
	// to a VTK-formatted plot file.
	{
		char vtk_name[200];
		sprintf(vtk_name, "out/%s_%d.exo", argv[1], (int)b_straight+(int)b_exact_bdry);
		ExodusII_IO exo_io(mesh, /*single_precision=*/true);

		// First plot the displacement field using a nodal plot
		std::set<std::string> system_names;
		system_names.insert("Poisson");
		exo_io.write_equation_systems(vtk_name,equation_systems,&system_names);
		// then append element-based discontinuous plots of the stresses
		exo_io.write_element_data(equation_systems);

	}

	// All done.
	return 0;
}

void mesh_straighten(MeshBase& mesh){

	// iterate over all elements
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el ; ++el){
		const Elem* elem = *el;

		for (unsigned int iside=0; iside<elem->n_sides(); iside++)
			if (elem->neighbor(iside) == NULL)
			{
				// get pointer to each node on the side
				UniquePtr<Elem> elem_side = elem->build_side(iside);
				libmesh_assert(elem_side->n_nodes() == 3);
				Node &n1 = *elem_side->get_node(0);
				Node &n2 = *elem_side->get_node(1);
				Node &n3 = *elem_side->get_node(2);

				// change node coordinates
				n3(0) = (n1(0) + n2(0))/2.;
				n3(1) = (n1(1) + n2(1))/2.;
			}
	}

}

void assemble_poisson(EquationSystems& es,
		const std::string& system_name)
{

	libmesh_assert_equal_to (system_name, "Poisson");

	const MeshBase& mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	LinearImplicitSystem& system = es.get_system<LinearImplicitSystem> ("Poisson");

	// dof_map
	const DofMap& dof_map = system.get_dof_map();

	//fe_type
	FEType fe_type = dof_map.variable_type(0);

	// FEBase for elements
	UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, FIFTH);
	fe->attach_quadrature_rule (&qrule);

	// FEBase for edges
	UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qface(dim-1, FIFTH);
	fe_face->attach_quadrature_rule (&qface);

	// The element Jacobian * quadrature weight at each integration point.
	const std::vector<Real>& JxW = fe->get_JxW();
	// The physical XY locations of the quadrature points on the element.
	const std::vector<Point>& q_point = fe->get_xyz();
	// The element shape functions evaluated at the quadrature points.
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	// The element shape function gradients evaluated at the quadrature
	// points.
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	// Local LHS mat and RHS vec
	DenseMatrix<Number> Ke;
	DenseVector<Number> Fe;

	// global dof for an element
	std::vector<dof_id_type> dof_indices;


	/*
	 * Loop over all the elements
	 */
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el ; ++el)
	{
		// fancy alias
		const Elem* elem = *el;

		// get dof indices.
		dof_map.dof_indices (elem, dof_indices);
		// Compute the element-specific data for the current element
		fe->reinit (elem);

		// Zero the element matrix and right-hand side before
		// summing them.
		Ke.resize (dof_indices.size(),
				dof_indices.size());
		Fe.resize (dof_indices.size());

		// Now loop over the quadrature points.
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{

			// double loop for integration - K
			for (unsigned int i=0; i<phi.size(); i++)
				for (unsigned int j=0; j<phi.size(); j++)
				{
					Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
				}

			// single loop for integration - F
			const Real x = q_point[qp](0);
			const Real y = q_point[qp](1);
			const Real fxy = right_hand_side(x,y);
			for (unsigned int i=0; i<phi.size(); i++)
				Fe(i) += JxW[qp]*fxy*phi[i][qp];
		}

		// loop over sides for boundary conditions
		for (unsigned int side=0; side<elem->n_sides(); side++)
			if (elem->neighbor(side) == NULL)
			{
				const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
				const std::vector<Real>& JxW_face = fe_face->get_JxW();
				const std::vector<Point >& qface_point = fe_face->get_xyz();
				fe_face->reinit(elem, side);

				// Loop over the face quadrature points for integration.
				for (unsigned int qp=0; qp<qface.n_points(); qp++)
				{
					const Real xf = qface_point[qp](0);
					const Real yf = qface_point[qp](1);
					const Real penalty = 1.e10;
					const Real value = boundary_value(xf, yf);

					for (unsigned int i=0; i<phi_face.size(); i++)
						for (unsigned int j=0; j<phi_face.size(); j++)
							Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

					for (unsigned int i=0; i<phi_face.size(); i++)
						Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
				}
			}

		// Assemble the local contribs
		system.matrix->add_matrix (Ke, dof_indices);
		system.rhs->add_vector    (Fe, dof_indices);
	}

	// All done!
}

void find_error(EquationSystems& es)
{
	const MeshBase& mesh = es.get_mesh();
	const int dim = mesh.mesh_dimension();

	LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Poisson");
	ExplicitSystem& err_system = es.get_system<ExplicitSystem>("Error");

	// indexing
	const DofMap& dof_map = system.get_dof_map();
	const DofMap& err_dof_map = err_system.get_dof_map();

	std::vector<dof_id_type> dof_indices;
	std::vector<dof_id_type> err_dof_indices;

	// finite element for quadrature
	FEType fe_type = dof_map.variable_type(0);

	UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, SIXTH);
	fe->attach_quadrature_rule (&qrule);

	const std::vector<Real>& JxW = fe->get_JxW();
	const std::vector<Point>& q_point = fe->get_xyz();

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
    	double error = 0;
    	const Elem* elem = *el;


    	// find the error
    	dof_map.dof_indices (elem, dof_indices);
    	fe->reinit (elem);
    	for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    	{
    		const Real x = q_point[qp](0);
    		const Real y = q_point[qp](1);
    		const Real u_exact = exact_solution(x,y);
    		const Real u_comp = system.point_value(0, q_point[qp], *elem);
    		const Real local_err = (u_exact - u_comp);

    		error += JxW[qp]*local_err*local_err;
    	}
		error = sqrt( error/elem->volume() );

    	// put the error in its place
    	err_dof_map.dof_indices (elem, err_dof_indices);
    	libmesh_assert(err_dof_indices.size() == 1);
    	err_system.solution->set(err_dof_indices[0], error);
    }

    err_system.solution->close();
    err_system.update();
}
