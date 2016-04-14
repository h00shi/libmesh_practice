#include "linear_elasticity.hxx"

#include "mpi.h"
#include <iomanip>


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

TestCaseBump::TestCaseBump(Material &material, const double r, const double v_max, const double c, const double b, Parallel::Communicator& comm):
			  TestCase(material, false, false),
			  _r(r),
			  _h(r*cos(v_max)),
			  _c(c),
			  _l(5 * r * sin(v_max)),
			  _b(b),
			  _last_projectee((Point *) NULL),
			  _comm(comm)
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

	// tolerances
	ierr = SNESSetTolerances(_snes,1e-16,1e-20,1e-20,100,1000);libmesh_assert(ierr==0);

	// do not set from options - when I say gmres I mean the global solver not this!!!!
	// ierr = SNESSetFromOptions(_snes);libmesh_assert(ierr==0);


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

PetscErrorCode TestCaseBump::snes_residual(SNES,Vec xvec,Vec f,void *ctx)
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

  // find n
  Point pn = g->normal(u,v);

  // find x, y and z
  Point pxyz = g->xyz(u,v);

  // evaluate residuals
  ff[0] = xs - pxyz(X) - d * pn(X);
  ff[1] = ys - pxyz(Y) - d * pn(Y);
  ff[2] = zs - pxyz(Z) - d * pn(Z);

  /* Restore vectors */
  ierr = VecRestoreArrayRead(xvec,&xx);libmesh_assert(ierr==0);
  ierr = VecRestoreArray(f,&ff);libmesh_assert(ierr==0);
  return 0;
}

Point TestCaseBump::project_on_bump(const Point& p)
{

	PetscErrorCode ierr;

	// check if we just projected this guy
	if (&p == _last_projectee)
	{
		//if(!_comm.rank()) printf("used_previous \n");
		return _last_projection;
	}
	else
	{

		// set the xs,ys,zs to be used by snes
		_xs = p(X); _ys=p(Y); _zs=p(Z);

		// Initial guess for u, v and d
		double *xm;
		ierr  = VecGetArray(_vec_x,&xm);libmesh_assert(ierr==0);
		u0v0(p,xm[0],xm[1]);
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

		p_proj = xyz(u,v);

		ierr = VecRestoreArrayRead(_vec_x,&xr);libmesh_assert(ierr==0);

		// check d to be consistent with the returned point
		libmesh_assert( fabs( fabs(d) - (p_proj-p).size() ) < 1e-12 );

		// save the projection
		_last_projectee = &p;
		_last_projection = p_proj;

		// return the projection
		return p_proj;
	}
}

void TestCaseBump::project_on_bump_test()
{

	// do this only for rank zero
	if (_comm.rank() > 0 ) return;

	// the exact solution in u and v
	const double ue = 0.3 * _l/2./_c;
	const double ve = 0.7 * M_PI/6;
	const double de = 0.1;

	// exact solution in x and y
	Point ne = normal(ue, ve);
	Point xe = xyz(ue, ve);

	// point to be projected
	Point p;
	p(X) = xe(X) + de*ne(X);
	p(Y) = xe(Y) + de*ne(Y);
	p(Z) = xe(Z) + de*ne(Z);

	// find the initial guess
	double u0, v0;
	u0v0(p,u0,v0);
	printf("Initial guess: %10.8lf %10.8lf \n", u0, v0);

	// project the point
	Point p_proj = project_on_bump(p);

    // print the results
	printf("Projection results: \n");
	printf(" init_x: %10.8lf %10.8lf %10.8lf\n", p(X), p(Y), p(Z));
	printf("exact_x: %10.8lf %10.8lf %10.8lf\n", xe(X), xe(Y), xe(Z));
	printf("  num_x: %10.8lf %10.8lf %10.8lf\n", p_proj(X) , p_proj(Y) , p_proj(Z) );
	printf("   err: %10e %10e %10e\n", p_proj(X)-xe(X) , p_proj(Y)-xe(Y) , p_proj(Z)-xe(Z) );
	printf("\n");

	// Exit after execution
	//libmesh_assert_msg(0,"This function causes assertion and exit!!");
}


void TestCaseBump::boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q)
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
			q = 0; // let u move
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

	// The surface that must be curved
	case 4:
	{
		// Find the new coordinates
		Point p_new = project_on_bump(p);

		// For debugging:
		// check correct convergence
		//printf("%.8lf \n", p_new(Z) - p(Z));

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

	default:
		libmesh_assert(0);

	//end of case
	}

}

/****************************************************************************
 *                             Assembler                                    *
 ****************************************************************************/

void LinearElasticity::post_process(std::ostream& outstream)
{


	// save old setting
	std::ios::fmtflags old_settings;
	if (es.comm().rank() == 0)
	{
		old_settings = outstream.flags();
		outstream.precision(10);
		outstream.setf(std::ios::scientific, std::ios::floatfield);
	}

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

		if (es.comm().rank() == 0)
		{
			outstream << esol.l2_error("Elasticity", "u") << " " <<esol.l2_error("Elasticity", "v") << " " << esol.l2_error("Elasticity", "w") << std::endl;
			outstream << esol.h1_error("Elasticity", "u") << " " <<esol.h1_error("Elasticity", "v") << " " << esol.h1_error("Elasticity", "w") << std::endl;
		}
	}
	// print the convergence history
	else
	{
		// get the object to the solver
		LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");

		PetscLinearSolver<Number>* petsc_linear_solver =
		    libmesh_cast_ptr<PetscLinearSolver<Number>*>(system.get_linear_solver());
		libmesh_assert(petsc_linear_solver);

		// get the convergence history
		std::vector<Real> history;
		petsc_linear_solver->get_residual_history(history);

		// print the histroy
		if (es.comm().rank() == 0)
		{
			for (uint i = 0 ; i < history.size() ; i++)
			{
				outstream << std::setw(10) << i+1  << " "
						<< std::setw(20) << history[i] << std::endl;
			}
		}

	}

	// restore flags
	if (es.comm().rank() == 0)
		outstream.flags(old_settings);
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
