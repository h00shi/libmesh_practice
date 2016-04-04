#include "linear_elasticity.hxx"


/****************************************************************************
 *                                Test Case                                 *
 ****************************************************************************/

/*
 * MMS functions
 */
double TestCaseMMS::exact_solution(const Point& p, const uint n_unknown)
{

	double ans;

	switch(n_unknown)
	{
	case U:
		ans = cos(M_PI*p(X))*cos(M_PI*p(Y))*cosh(M_PI*p(Z));
		break;
	case V:
		ans = sin(M_PI*p(X))*sin(2*M_PI*p(Y))*sinh(M_PI*p(Z));
		break;
	case W:
		ans = sin(2*M_PI*p(X))*sin(M_PI*p(Y))*sinh(2*M_PI*p(Z));
		break;
	default:
		libmesh_assert(0);
		break;
	}

	return ans;
}

Gradient TestCaseMMS::exact_derivative(const Point& p, const uint n_unknown)
{

	Gradient g;
	switch(n_unknown)
	{
	case U:
		g(X) = -sin(M_PI * p(X)) * M_PI * cos(M_PI * p(Y)) * cosh(M_PI * p(Z));
		g(Y) = -cos(M_PI * p(X)) * sin(M_PI * p(Y)) * M_PI * cosh(M_PI * p(Z));
		g(Z) = cos(M_PI * p(X)) * cos(M_PI * p(Y)) * sinh(M_PI * p(Z)) * M_PI;
		break;
	case V:
		g(X) = cos(M_PI * p(X)) * M_PI * sin(2. * M_PI * p(Y)) * sinh(M_PI * p(Z));
		g(Y) = 2. * sin(M_PI * p(X)) * cos(2. * M_PI * p(Y)) * M_PI * sinh(M_PI * p(Z));
		g(Z) = sin(M_PI * p(X)) * sin(2. * M_PI * p(Y)) * cosh(M_PI * p(Z)) * M_PI;
		break;
	case W:
		g(X) = 2. * cos(2. * M_PI * p(X)) * M_PI * sin(M_PI * p(Y)) * sinh(2 * M_PI * p(Z));
		g(Y) = sin(2. * M_PI * p(X)) * cos(M_PI * p(Y)) * M_PI * sinh(2 * M_PI * p(Z));
		g(Z) = 2. * sin(2. * M_PI * p(X)) * sin(M_PI * p(Y)) * cosh(2. * M_PI * p(Z)) * M_PI;
		break;
	default:
		libmesh_assert(0);
	}

	return g;
}

double TestCaseMMS::right_hand_side(const Point& p, const uint n_unknown)
{
	double f;
	const double x=p(X), y=p(Y), z=p(Z);
	const double pi=M_PI, mu=_material.mu, lambda=_material.lambda;

	switch(n_unknown)
	{
	case U:
		f = (double) (2 * mu + lambda) * cos(pi * x) * pi * pi * cos(pi * y) * cosh(pi * z) - 0.2e1 * (double) lambda * cos(pi * x) * pi * pi * cos(0.2e1 * pi * y) * sinh(pi * z) - 0.4e1 * (double) lambda * cos(0.2e1 * pi * x) * pi * pi * sin(pi * y) * cosh(0.2e1 * pi * z) - (double) mu * (-cos(pi * x) * cos(pi * y) * pi * pi * cosh(pi * z) + 0.2e1 * cos(pi * x) * pi * pi * cos(0.2e1 * pi * y) * sinh(pi * z)) - (double) mu * (cos(pi * x) * cos(pi * y) * pi * pi * cosh(pi * z) + 0.4e1 * cos(0.2e1 * pi * x) * pi * pi * sin(pi * y) * cosh(0.2e1 * pi * z));
		break;
	case V:
		f = -mu * (sin(pi * x) * pi * pi * sin(pi * y) * cosh(pi * z) - sin(pi * x) * pi * pi * sin(0.2e1 * pi * y) * sinh(pi * z)) - lambda * sin(pi * x) * pi * pi * sin(pi * y) * cosh(pi * z) + 0.4e1 * (0.2e1 * mu + lambda) * sin(pi * x) * sin(0.2e1 * pi * y) * pi * pi * sinh(pi * z) - 0.2e1 * lambda * sin(0.2e1 * pi * x) * cos(pi * y) * pi * pi * cosh(0.2e1 * pi * z) - mu * (sin(pi * x) * pi * pi * sin(0.2e1 * pi * y) * sinh(pi * z) + 0.2e1 * sin(0.2e1 * pi * x) * cos(pi * y) * pi * pi * cosh(0.2e1 * pi * z));
		break;
	case W:
		f = -mu * (-sin(pi * x) * pi * pi * cos(pi * y) * sinh(pi * z) - 0.4e1 * sin(0.2e1 * pi * x) * pi * pi * sin(pi * y) * sinh(0.2e1 * pi * z)) - mu * (0.2e1 * sin(pi * x) * cos(0.2e1 * pi * y) * pi * pi * cosh(pi * z) - sin(0.2e1 * pi * x) * pi * pi * sin(pi * y) * sinh(0.2e1 * pi * z)) + lambda * sin(pi * x) * pi * pi * cos(pi * y) * sinh(pi * z) - 0.2e1 * lambda * sin(pi * x) * cos(0.2e1 * pi * y) * pi * pi * cosh(pi * z) - 0.4e1 * (0.2e1 * mu + lambda) * sin(0.2e1 * pi * x) * sin(pi * y) * sinh(0.2e1 * pi * z) * pi * pi;
		break;
	default:
		libmesh_assert(0);
	}

	return f;
}

void TestCaseMMS::boundary_penalty(const Point& p, const uint , const uint n_unknown, double &t, double &q)
{
	const double penalty = 1.e10;
	q = penalty;
	t = penalty*exact_solution(p,n_unknown);
}

void TestCaseSphere::boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q)
{
	const double penalty = 1.e10;

	switch(id)
	{
	// Perpendicular to x
	case 1:
	{
		// No U displacement
		if (n_unknown == U)
		{
			q = penalty;
			t = 0;
		}
		// Others are free to move
		else
		{
			q = t = 0;
		}
		break;
	}

	// Perpendicular to y
	case 2:
	{
		// No V displacement
		if (n_unknown == V)
		{
			q = penalty;
			t = 0;
		}
		// Others are free to move
		else
		{
			q = t = 0;
		}
		break;
	}

	// Perpendicular to z
	case 3:
	{
		// No W displacement
		if (n_unknown == W)
		{
			q = penalty;
			t = 0;
		}
		// Others are free to move
		else
		{
			q = t = 0;
		}
		break;
	}

	// Outer surface
	case 4:
	{
		// Stop everyone from moving
		q = penalty;
		t = 0;
		break;
	}

	// The inner surface that must be curved
	case 5:
	{
		// Find the new coordinates
		Point p_new = _r_in / p.size() * p;

		// Force the displacement
		q = penalty;

		switch(n_unknown)
		{
		case U: t = (p_new(X) - p(X)) * penalty; break;
		case V: t = (p_new(Y) - p(Y)) * penalty; break;
		case W: t = (p_new(Z) - p(Z)) * penalty; break;
		}

		break;
	}

	//end of case
	}

}


/****************************************************************************
 *                             Assembler                                    *
 ****************************************************************************/

void LinearElasticity::post_process()
{
	// find the errors
	if(_test_case.has_exact_solution())
	{
		// Wrappers for the exact solution functions
		ExactValue eval(_test_case);
		ExactDeriv eder(_test_case);

		// Exact solution object
		ExactSolution esol(es);
		esol.attach_exact_value(0, &eval);
		esol.attach_exact_deriv(0, &eder);
		esol.extra_quadrature_order(1);

		esol.compute_error("Elasticity", "u");
		esol.compute_error("Elasticity", "v");
		esol.compute_error("Elasticity", "w");

		std::cout << esol.l2_error("Elasticity", "u") << " " <<esol.l2_error("Elasticity", "v") << " " << esol.l2_error("Elasticity", "w") << std::endl;
		std::cout << esol.h1_error("Elasticity", "u") << " " <<esol.h1_error("Elasticity", "v") << " " << esol.h1_error("Elasticity", "w") << std::endl;
	}
}


/**
 * Assemble the system matrix and right-hand side vector.
 */
void LinearElasticity::assemble()
{
	// mesh and dimensions
	const MeshBase& mesh = es.get_mesh();
	const uint dim = mesh.mesh_dimension();

	// the system and var numbers
	LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");
	const uint u_var = system.variable_number ("u");

	// Dof Map
	const DofMap& dof_map = system.get_dof_map();

	// FE for volume integrals
	FEType fe_type = dof_map.variable_type(u_var);
	UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, fe_type.default_quadrature_order());
	fe->attach_quadrature_rule (&qrule);

	// shape functions and integration weights
	const std::vector<Real>& JxW = fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	// FE for boundary integrals
	UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	QGauss qface(dim-1, fe_type.default_quadrature_order());
	fe_face->attach_quadrature_rule (&qface);

	// Shape function and integration weigths
	const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
	const std::vector<Real>& JxW_face = fe_face->get_JxW();

	// Location of quadrature points
	const std::vector<Point>& q_point = fe->get_xyz();
	const std::vector<Point >& qface_point = fe_face->get_xyz();

	// Right hand side mats and vectors
	DenseMatrix<Number> Ke;
	DenseSubMatrix<Number> Ke_var[3][3] =
	{
			{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
			{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
			{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
	};

	DenseVector<Number> Fe;
	DenseSubVector<Number> Fe_var[3] =
	{DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)};

	// Dof indices
	std::vector<dof_id_type> dof_indices;
	std::vector< std::vector<dof_id_type> > dof_indices_var(3);

	/*
	 * For over all elements
	 */
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el; ++el)
	{
		// pointer to element
		const Elem* elem = *el;

		// get dof indices
		dof_map.dof_indices (elem, dof_indices);
		for(unsigned int var=0; var<3; var++)
		{
			dof_map.dof_indices (elem, dof_indices_var[var], var);
		}

		const unsigned int n_dofs   = dof_indices.size();
		const unsigned int n_var_dofs = dof_indices_var[0].size();

		// reinit element
		fe->reinit (elem);

		// position the vector and matrices
		Ke.resize (n_dofs,n_dofs);
		for(unsigned int var_i=0; var_i<3; var_i++)
			for(unsigned int var_j=0; var_j<3; var_j++)
			{
				Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);
			}

		Fe.resize (n_dofs);
		for(unsigned int var=0; var<3; var++)
		{
			Fe_var[var].reposition (var*n_var_dofs, n_var_dofs);
		}

		/*
		 * Volume quadrature points
		 */
		 for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		 {

			 // Assemble the Jacobian
			 for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
				 for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
				 {
					 for(unsigned int i=0; i<3; i++)
						 for(unsigned int j=0; j<3; j++)
							 for(unsigned int k=0; k<3; k++)
								 for(unsigned int l=0; l<3; l++)
								 {
									 Ke_var[i][j](dof_i,dof_j) += JxW[qp] *
											 elasticity_tensor(i,j,k,l) *
											 dphi[dof_i][qp](k)* dphi[dof_j][qp](l);
								 }
				 }

			 // Assemble the right hand side
			 DenseVector<Number> f_vec(3);
			 if (_test_case.has_nonzero_rhs())
			 {
				 f_vec(0) = _test_case.right_hand_side(q_point[qp],0);
				 f_vec(1) = _test_case.right_hand_side(q_point[qp],1);
				 f_vec(2) = _test_case.right_hand_side(q_point[qp],2);
			 }
			 else
			 {
				 f_vec(0) = f_vec(1) = f_vec(2) = 0;
			 }
			 for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
			 {
				 for(unsigned int i=0; i<3; i++)
				 {
					 Fe_var[i](dof_i) += JxW[qp] *
							 ( f_vec(i) * phi[dof_i][qp] );
				 }
			 }
		 } // End of volume quadrature points


		 /*
		  * Boundary quadrature points for the penalty method
		  */
		 DenseVector<Number> t_vec(3), q_vec(3);

		 // find a boundary side
		 for (unsigned int side=0; side<elem->n_sides(); side++)
			 if (elem->neighbor(side) == NULL)
			 {
				 fe_face->reinit(elem, side);

				 for (unsigned int qp=0; qp<qface.n_points(); qp++)
				 {
					 // Find the boundary id
					 const uint id = mesh.get_boundary_info().boundary_id(elem,side);

					 // Get the penalty parameters
					 for (uint i=0; i<3 ; i++)
					 {
						 _test_case.boundary_penalty(qface_point[qp],id,i,t_vec(i),q_vec(i));
					 }

					 // Add contribution to RHS
					 for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
					 {
						 for(unsigned int i=0; i<3; i++)
						 {

							 Fe_var[i](dof_i) += JxW_face[qp] *
									 ( t_vec(i) * phi_face[dof_i][qp] );
						 }
					 }

					 // Add contribution to LHS
					 for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
						 for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
						 {
							 for(unsigned int i=0; i<3; i++)
							 {
								 Ke_var[i][i](dof_i,dof_j) += JxW_face[qp] *
										 q_vec(i) * phi_face[dof_j][qp]*phi_face[dof_i][qp];
							 }
						 }

				 } // end of quadrature points
			 } // end of boundary side

		 system.matrix->add_matrix (Ke, dof_indices);
		 system.rhs->add_vector    (Fe, dof_indices);
	} // end of elements
}
