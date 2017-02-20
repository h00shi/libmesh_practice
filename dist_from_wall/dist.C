// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

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

// assemble the system
void assemble_poisson(EquationSystems& es, const std::string& system_name);
void post_process(EquationSystems& es);


extern Number right_hand_side_interior(const Point& ,const Parameters&,const std::string&,const std::string&);
extern Gradient right_hand_side_bdry(const Point& ,const Parameters& ,const std::string& ,const std::string&,const int);
extern Number exact_solution(const Point& ,const Parameters&,const std::string&,const std::string&);
extern Gradient exact_deriv(const Point& ,const Parameters&,const std::string&,const std::string&);


int main (int argc, char** argv)
{
	/*
	 * Init Libmesh
	 */
	LibMeshInit init (argc, argv);
	GetPot getpot(argc, argv);

	std::string mesh_name = getpot.follow("UNKNOWN.msh", "-m");
	std::string unk_name = getpot.follow("u", "-v");
	mesh_name = mesh_name.substr(0, mesh_name.rfind(".msh"));
	std::stringstream ss;

	std::cout << "mesh name is " << mesh_name << std::endl;

	/*
	 * Print the command line.
	 */
	std::cout << "Running " << argv[0];
	for (int i=1; i<argc; i++)
		std::cout << " " << argv[i];
	std::cout << std::endl << std::endl;
	if(argc < 2) libmesh_error_msg("Mesh file not given.");

	// create a mesh
	Mesh mesh(init.comm());
	GmshIO gmsh(mesh);
	ss.str(""); ss << mesh_name << ".msh";
	gmsh.read(ss.str());
	mesh.all_second_order();
	mesh.prepare_for_use();
	mesh.print_info();

	// define the system for u
	EquationSystems equation_systems (mesh);
	equation_systems.add_system<LinearImplicitSystem> ("Poisson");
	equation_systems.get_system("Poisson").add_variable(unk_name, SECOND, HIERARCHIC);
	equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
	// define error and distance from wall
	equation_systems.add_system<ExplicitSystem> ("PostProcess");
	equation_systems.get_system("PostProcess").add_variable("p1", CONSTANT, MONOMIAL);
	equation_systems.get_system("PostProcess").add_variable("p2", CONSTANT, MONOMIAL);

	// solve and find the error
	equation_systems.init();
	equation_systems.print_info();
	equation_systems.get_system("Poisson").solve();
	post_process(equation_systems);


	// find the exact solution
	if(unk_name == "u")
	{
		ExactSolution exact_sol(equation_systems);
		exact_sol.attach_exact_value(exact_solution);
		exact_sol.attach_exact_deriv(exact_deriv);
		exact_sol.extra_quadrature_order(1);
		exact_sol.compute_error("Poisson", "u");
		std::cout << "H^1 error: " << exact_sol.h1_error("Poisson", "u") << std::endl;
		std::cout << "L_2 error: " << exact_sol.l2_error("Poisson", "u") << std::endl;
	}

	// After solving the system write the solution
	// to a VTK-formatted plot file.
	ss.str(""); ss <<  mesh_name << ".exo";
	ExodusII_IO exo_io(mesh, /*single_precision=*/true);

	// First plot the displacement field using a nodal plot
	std::set<std::string> system_names;
	system_names.insert("Poisson");
	//system_names.insert("PostProcess");
	exo_io.write_equation_systems(ss.str(), equation_systems, &system_names);
	// then append element-based discontinuous plots of the stresses
	 exo_io.write_element_data(equation_systems);

	// All done.
	return 0;
}


void assemble_poisson(EquationSystems& es,
		const std::string& system_name)
{

	libmesh_assert_equal_to (system_name, "Poisson");

	const MeshBase& mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	LinearImplicitSystem& system = es.get_system<LinearImplicitSystem> ("Poisson");
	const std::string unk_name = system.variable_name(0);

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
			const Real fxy = right_hand_side_interior(Point(x,y), es.parameters, "Poisson", unk_name);
			for (unsigned int i=0; i<phi.size(); i++)
				Fe(i) += JxW[qp]*fxy*phi[i][qp];
		}

		// loop over sides for boundary conditions
		for (unsigned int side=0; side<elem->n_sides(); side++)
		{
			if ( elem->neighbor(side) == NULL )
			{
				const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
				const std::vector<Real>& JxW_face = fe_face->get_JxW();
				const std::vector<Point >& qface_point = fe_face->get_xyz();
				fe_face->reinit(elem, side);
				const int bid =mesh.boundary_info->boundary_id(elem,side);

				// Loop over the face quadrature points for integration.
				for (unsigned int qp=0; qp<qface.n_points(); qp++)
				{
					Gradient gr;
					gr = right_hand_side_bdry(qface_point[qp], es.parameters, "Poisson", unk_name, bid);

					for (unsigned int i=0; i<phi_face.size(); i++)
						for (unsigned int j=0; j<phi_face.size(); j++)
							Ke(i,j) += JxW_face[qp]*gr(0)*phi_face[i][qp]*phi_face[j][qp];

					for (unsigned int i=0; i<phi_face.size(); i++)
						Fe(i) += JxW_face[qp]*gr(1)*phi_face[i][qp];
				}
			}
		}

//		for (unsigned int side=0; side<elem->n_sides(); side++)
//			if (elem->neighbor(side) == NULL)
//			{
//				std::cout << mesh.boundary_info->boundary_id (elem,side) << std::endl;
//			}
		// Assemble the local contribs
		//dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices); //what is this!!??!
		//dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
		system.matrix->add_matrix (Ke, dof_indices);
		system.rhs->add_vector    (Fe, dof_indices);
	}

	// All done!
}

void post_process(EquationSystems& es)
{
	const MeshBase& mesh = es.get_mesh();
	const int dim = mesh.mesh_dimension();

	LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Poisson");
	ExplicitSystem& err_system = es.get_system<ExplicitSystem>("PostProcess");
	std::string unk_name = system.variable_name(0);

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
    	const Elem* elem = *el;

    	dof_map.dof_indices (elem, dof_indices);
    	fe->reinit (elem);

    	if(unk_name == "u")
    	{
    		double error = 0;

    		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    		{
    			const Real u_exact = exact_solution(q_point[qp],es.parameters, "Poisson", unk_name);
    			const Real u_comp = system.point_value(0, q_point[qp], *elem);
    			const Real local_err = (u_exact - u_comp);

    			error += JxW[qp]*local_err*local_err;
    		}
    		error = sqrt( error/elem->volume() );

    		// put the error in its place
    		err_dof_map.dof_indices (elem, err_dof_indices);
    		libmesh_assert(err_dof_indices.size() == 2);
    		err_system.solution->set(err_dof_indices[0], error);
      		err_system.solution->set(err_dof_indices[1], 0);
    	}
    	else
    	{
    		libmesh_assert_equal_to(unk_name, "phi");
    		double dist_ave1 = 0;
    		double dist_ave2 = 0;
    		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    		{
    			const Number phi      = system.point_value(0, q_point[qp], *elem);
    			const Gradient grad_phi = system.point_gradient(0, q_point[qp], *elem);
    			const Real grad_phi_mag2 = grad_phi.norm_sq();
    			const Real grad_phi_mag = sqrt(grad_phi_mag2);
    			const Real dist1 = -grad_phi_mag + sqrt(grad_phi_mag2 + 2*phi);
    			const Real dist2 = -grad_phi_mag - sqrt(grad_phi_mag2 + 2*phi);

    			dist_ave1 += JxW[qp]*dist1;
    			dist_ave2 += JxW[qp]*dist2;
    		}
    		dist_ave1 /= elem->volume() ;
    		dist_ave2 /= elem->volume() ;

    		err_dof_map.dof_indices (elem, err_dof_indices);
    		libmesh_assert(err_dof_indices.size() == 2);
    		err_system.solution->set(err_dof_indices[0], dist_ave1);
      		err_system.solution->set(err_dof_indices[1], dist_ave2);
    	}
    }

    err_system.solution->close();
    err_system.update();
}
