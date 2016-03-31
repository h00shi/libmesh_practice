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
#include "libmesh/gmsh_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

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
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Exact solution
#include "libmesh/exact_solution.h"
#include "libmesh/mesh_function.h"


// Get pot
#include "libmesh/getpot.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

/****************************************************************************
 *                         Enums and stuff                                  *
 ****************************************************************************/

/****************************************************************************
 *                         Functions                                        *
 ****************************************************************************/

//Assembly function
void assemble_steady(EquationSystems& es, const std::string& system_name);

/****************************************************************************
 *                             Main                                         *
 ****************************************************************************/

int main (int argc, char** argv)
{

	// Initialize libraries, like in example 2.
	LibMeshInit init (argc, argv);
	GetPot getpot(argc, argv);

	// important data
	const int n_elem = getpot("NX", 100);
	std::string mesh_name = getpot("mesh","square");
	const bool exo_io = getpot("exo",false);
	Order fe_order=THIRD;
	FEFamily fe_family= BERNSTEIN;

	// create a mesh
	Mesh mesh(init.comm());
	if (mesh_name == "line")
	{
		MeshTools::Generation::build_line(mesh,n_elem,0.,M_PI,EDGE2);
	}
	else if (mesh_name == "square")
	{
		MeshTools::Generation::build_square(mesh,n_elem,n_elem,
				0.,M_PI,
				0.,M_PI,
				QUAD9, true);
	}
	else
	{
		GmshIO gmsh(mesh);
		gmsh.read(mesh_name);
		mesh.prepare_for_use();
	}
	mesh.print_info();

	// define the system for u
	EquationSystems equation_systems (mesh);

	// define the steady state solution
	LinearImplicitSystem &system_steady =
	equation_systems.add_system<LinearImplicitSystem> ("Steady");
	system_steady.add_variable("sig", fe_order, fe_family);
	system_steady.add_variable("tau", fe_order, fe_family);
	system_steady.attach_assemble_function (assemble_steady);

	// init system
	equation_systems.init ();
	equation_systems.print_info();

	// solve the steady state
	system_steady.solve();

	//write the first time step
	std::string exodus_filename = "transient_ex1.e";
	if(exo_io)
	{
		std::set<std::string> system_names;
		system_names.insert("Steady");
		std::string exodus_filename = "transient_ex1.e";
		ExodusII_IO(mesh).write_equation_systems (exodus_filename, equation_systems, &system_names);
	}

	  // All done.
	  return 0;
}


/****************************************************************************
 *                          Steady State                                    *
 ****************************************************************************/

void assemble_steady(EquationSystems& es, const std::string&)
{

	// mesh stuff
	const MeshBase& mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	// system and variable numbers
	LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Steady");
	const uint sig_var = system.variable_number ("sig");
	const uint tau_var = system.variable_number ("tau");

	// Finite element and quadrature
	FEType fe_type = system.variable_type(sig_var);
	UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));
	UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qrule (dim,   TENTH);
	QGauss qface (dim-1, TENTH);
	fe->attach_quadrature_rule      (&qrule);
	fe_face->attach_quadrature_rule (&qface);

	// Value of shape function and its derivatives
	const std::vector<Real>& JxW      = fe->get_JxW();
	const std::vector<Real>& JxW_face = fe_face->get_JxW();
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<Real> >& psi = fe_face->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	// dof map and its indices
	const DofMap & dof_map = system.get_dof_map();
	std::vector<dof_id_type> dof_indices;
	std::vector<dof_id_type> dof_indices_sig;
	std::vector<dof_id_type> dof_indices_tau;

	// Vectors and subvectors
	DenseMatrix<Number> Ke;
	DenseVector<Number> Fe;
	DenseSubMatrix<Number>  Kss(Ke), Kst(Ke), Kts(Ke), Ktt(Ke);
	DenseSubVector<Number>  Fs(Fe), Ft(Fe);

	/*
	 * Assembly loop over elements
	 */
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el; ++el)
	{
		// store pointer to element
		const Elem* elem = *el;

		// Find the dof maps
		dof_map.dof_indices (elem, dof_indices);
		dof_map.dof_indices (elem, dof_indices_sig, sig_var);
		dof_map.dof_indices (elem, dof_indices_tau, tau_var);

		// find number of dof's
		const unsigned int n_dofs   = dof_indices.size();
		const unsigned int n_sig_dofs = dof_indices_sig.size();

		// Compute the element specific stuff
		fe->reinit  (elem);

		// Zero the matrices
		Ke.resize (n_dofs, n_dofs);
		Fe.resize (n_dofs);

		// Reposition the submatrices...  The idea is this:
		//        -        -          -  -
		//        | Kss Kst |        | Fs |
		//   Ke = | Kts Ktt |;  Fe = | Ft |
		//        -         -         -  -
		Kss.reposition (sig_var*n_sig_dofs, sig_var*n_sig_dofs, n_sig_dofs, n_sig_dofs);
		Kst.reposition (sig_var*n_sig_dofs, tau_var*n_sig_dofs, n_sig_dofs, n_sig_dofs);

		Kts.reposition (tau_var*n_sig_dofs, sig_var*n_sig_dofs, n_sig_dofs, n_sig_dofs);
		Ktt.reposition (tau_var*n_sig_dofs, tau_var*n_sig_dofs, n_sig_dofs, n_sig_dofs);

		Fs.reposition (sig_var*n_sig_dofs, n_sig_dofs);
		Ft.reposition (tau_var*n_sig_dofs, n_sig_dofs);

		// Sum over all quadrature points
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{

			for (unsigned int i=0; i<n_sig_dofs; i++)
				Fs(i) +=  JxW[qp]*phi[i][qp];

			for (unsigned int i=0; i<n_sig_dofs; i++)
				for (unsigned int j=0; j<n_sig_dofs; j++)
				{
					// sig sig
					Kss(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];

					// sig tau
					Kst(i,j) += -JxW[qp]*phi[i][qp]*phi[j][qp];

					// tau sig
					Kts(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];

					// tau tau
					Ktt(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];
				}
		} // end of the quadrature point qp-loop

		/*
		 * Loop over sides four bnd conditions
		 */
		for (unsigned int s=0; s<elem->n_sides(); s++)
		{
			if (elem->neighbor(s) == NULL)
			{
				const Real penalty = 1.e10;

				fe_face->reinit(elem,s);

				for (unsigned int qp=0; qp<qface.n_points(); qp++)
				{
					// Dirichlet 0 right hand size
					const Number value = 0;

					// RHS contribution
					for (unsigned int i=0; i<psi.size(); i++)
					{
						Fs(i) += penalty*JxW_face[qp]*value*psi[i][qp];
						Ft(i) += penalty*JxW_face[qp]*value*psi[i][qp];
					}

					// Matrix contribution
					for (unsigned int i=0; i<psi.size(); i++)
						for (unsigned int j=0; j<psi.size(); j++)
						{
							Kss(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
							Ktt(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
						}
				}
			}
		} // end of bnd condition

		system.matrix->add_matrix (Ke, dof_indices);
		system.rhs->add_vector    (Fe, dof_indices);
	} // end of element loop

	return;
}

