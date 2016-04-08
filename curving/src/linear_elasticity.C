#include "linear_elasticity.hxx"

#include "mpi.h"


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

	// The upper surface is freee!!!
	case 6:
		t=q=0;
		break;

	default:
		libmesh_assert(0);

	//end of case
	}

}

TestCaseBump::TestCaseBump(Material &material, const double r, const double h, const double c):
			  TestCase(material, false, false),
			  _r(r),
			  _h(h),
			  _c(c)
{
	PetscErrorCode ierr;
	const uint n_unknown = 3;

	// create the snes for solving projection problems
	ierr = SNESCreate(PETSC_COMM_SELF,&_snes);libmesh_assert(ierr==0);

	// create the solution and residual norms
	ierr = VecCreateSeq(PETSC_COMM_SELF, n_unknown, &_vec_x);libmesh_assert(ierr==0);
	ierr = VecDuplicate(_vec_x,&_vec_r);libmesh_assert(ierr==0);

	// create the matrix
	ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, 3, NULL, &_mat_j);libmesh_assert(ierr==0);
	{
		int idx[20]; double ones[100];
		for (uint i=0; i<n_unknown; i++) idx[i]=i;
		for (uint i=0; i<n_unknown*n_unknown; i++) ones[i]=1;
		ierr=MatSetValues(_mat_j, n_unknown, idx, n_unknown, idx, ones, INSERT_VALUES);libmesh_assert(ierr==0);
	}
	ierr = MatAssemblyBegin(_mat_j, MAT_FINAL_ASSEMBLY); libmesh_assert(ierr==0);
	ierr = MatAssemblyEnd(_mat_j, MAT_FINAL_ASSEMBLY); libmesh_assert(ierr==0);

	// set the snes residual function
	ierr = SNESSetFunction(_snes,_vec_r,TestCaseBump::snes_residual, this);libmesh_assert(ierr==0);

	// set the snes jacobian function (normal finite difference)
	ierr = SNESSetJacobian(_snes,_mat_j,_mat_j,SNESComputeJacobianDefault, this);libmesh_assert(ierr==0);

	// other snes options
	KSP ksp;
	PC pc;

	ierr = SNESGetKSP(_snes,&ksp);libmesh_assert(ierr==0);
	ierr = KSPGetPC(ksp,&pc);libmesh_assert(ierr==0);
	ierr = KSPSetType(ksp, KSPPREONLY);libmesh_assert(ierr==0);
	ierr = PCSetType(pc,PCLU);libmesh_assert(ierr==0);

	// Strice tolerances
	ierr = SNESSetTolerances(_snes,1e-16,1e-20,1e-20,100,1000);libmesh_assert(ierr==0);
	ierr = SNESSetFromOptions(_snes);libmesh_assert(ierr==0);


}

TestCaseBump::~TestCaseBump()
{
	PetscErrorCode ierr;

	// destroy all the petsc stuff
	ierr = VecDestroy(&_vec_x);libmesh_assert(ierr==0);
	ierr = VecDestroy(&_vec_r);libmesh_assert(ierr==0);
	ierr = MatDestroy(&_mat_j);libmesh_assert(ierr==0);
	ierr = SNESDestroy(&_snes);libmesh_assert(ierr==0);
}

PetscErrorCode TestCaseBump::snes_residual(SNES snes,Vec xvec,Vec f,void *ctx)
{
  PetscErrorCode    ierr;
  const PetscScalar *xx;
  PetscScalar       *ff;

  // get pointer to data
  ierr = VecGetArrayRead(xvec,&xx);libmesh_assert(ierr==0);
  ierr = VecGetArray(f,&ff);libmesh_assert(ierr==0);

  /* extract data */
  const TestCaseBump *g = (TestCaseBump *) ctx;

  // xs, ys and zs
  const double xs = g->_xs;
  const double ys = g->_ys;
  const double zs = g->_zs;

  // u, v and d
  const double u = xx[0];
  const double v = xx[1];
  const double d = xx[2];

  // r, h and c
  const double r = g->_r;
  const double h = g->_h;
  const double c = g->_c;

  // find n
  const double denom_n = sqrt(c*c + sin(v)*sin(v) );
  const double n_x = c * sin(v) / denom_n;
  const double n_y = -sin(v) / denom_n;
  const double n_z = c * cos(v) / denom_n;

  // find x, y and z
  const double x = u + r * sin(v);
  const double y = c * u;
  const double z = r * cos(v) - h;

  // evaluate residuals
  ff[0] = xs - x - d * n_x;
  ff[1] = ys - y - d * n_y;
  ff[2] = zs - z - d * n_z;

  /* Restore vectors */
  ierr = VecRestoreArrayRead(xvec,&xx);libmesh_assert(ierr==0);
  ierr = VecRestoreArray(f,&ff);libmesh_assert(ierr==0);
  return 0;
}

Point TestCaseBump::project_on_bump(const Point& p)
{

	PetscErrorCode ierr;

	// set the xs,ys,zs to be used by snes
	_xs = p(X); _ys=p(Y); _zs=p(Z);

	// Initial guess for u, v and d
	double *xm;
    ierr  = VecGetArray(_vec_x,&xm);libmesh_assert(ierr==0);
	xm[0]= _ys / _c;
	xm[1]=  asin( (_xs - xm[0]) / _r );
	xm[2]=  0;
    ierr  = VecRestoreArray(_vec_x,&xm);libmesh_assert(ierr==0);

    // solve the system to get u, v, d
    ierr = SNESSolve(_snes,NULL,_vec_x);libmesh_assert(ierr==0);

    // create the point and return it
    const double *xr;
    Point p_proj;
    ierr = VecGetArrayRead(_vec_x,&xr);libmesh_assert(ierr==0);

    const double u = xr[0];
    const double v = xr[1];
    const double d = xr[2];

    p_proj(X) = u + _r * sin(v);
    p_proj(Y) = _c * u;
    p_proj(Z) = _r * cos(v) - _h;

    ierr = VecRestoreArrayRead(_vec_x,&xr);libmesh_assert(ierr==0);

    // check d to be consistent with the returned point
    libmesh_assert( fabs( fabs(d) - (p_proj-p).size() ) < 1e-13 );

    // return the point
    return p_proj;
}

void TestCaseBump::project_on_bump_test()
{
	// Hard code some values for c,r and h
	const double c = _c;
	const double r = _r;
	const double h = _h;
	_c = 1;
	_r = 1;
	_h = 0.5;

	// the exact solution in u and v
	const double ue = 0.5;
	const double ve = 0.7 * M_PI/3;
	const double de = 0.25;

	// exact solution in x and y
	const double denom_n = sqrt(_c * _c + sin(ve) * sin(ve) );
	const double n_x = _c * sin(ve) / denom_n;
	const double n_y = -sin(ve) / denom_n;
	const double n_z = _c * cos(ve) / denom_n;

	const double xe = ue + _r * sin(ve);
	const double ye = _c * ue;
	const double ze = _r * cos(ve) - _h;

	// point to be projected
	Point p;
	p(X) = xe + de*n_x;
	p(Y) = ye + de*n_y;
	p(Z) = ze + de*n_z;

	// project the point
	Point p_proj = project_on_bump(p);

    // print the results
	printf("Projection results: \n");
	printf(" init_x: %10.8lf %10.8lf %10.8lf\n", p(X), p(Y), p(Z));
	printf("exact_x: %10.8lf %10.8lf %10.8lf\n", xe, ye, ze);
	printf("  num_x: %10.8lf %10.8lf %10.8lf\n", p_proj(X) , p_proj(Y) , p_proj(Z) );
	printf("   err: %10e %10e %10e\n", p_proj(X)-xe , p_proj(Y)-ye , p_proj(Z)-ze );
	printf("\n");

	// Restore the values for r, c and h
	_c = c;
	_r = r;
	_h = h;

	// Exit after execution
	//libmesh_assert_msg(0,"This function causes assertion and exit!!");
}

void TestCaseBump::boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q)
{
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
