// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/exact_solution.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

/****************************************************************************
 *                       Petsc Configurations                               *
 ****************************************************************************/
class PetscSolverConfiguration : public SolverConfiguration
{
public:

  PetscSolverConfiguration(PetscLinearSolver<Number>& petsc_linear_solver)
  :
  _petsc_linear_solver(petsc_linear_solver)
  {
  }

  virtual void configure_solver()
  {
    PetscErrorCode ierr = 0;
    ierr = KSPSetType (_petsc_linear_solver.ksp(), const_cast<KSPType>(KSPCG));
    libmesh_assert(ierr == 0);

    ierr = PCSetType (_petsc_linear_solver.pc(), const_cast<PCType>(PCBJACOBI));
    libmesh_assert(ierr == 0);
  }

  // The linear solver object that we are configuring
  PetscLinearSolver<Number>& _petsc_linear_solver;

};

/****************************************************************************
 *                           Exact solutions                                *
 ****************************************************************************/
enum Variables{U=0, V, W};
enum Coordinates{X=0, Y, Z};

class TestCase
{
	  double _lambda, _mu;
public:

	  TestCase(const double lambda, const double mu):
		  _lambda(lambda),
		  _mu(mu)
		{}

	  /**
	    * Exact solution.
	    */
	  double exact_solution(const Point& p, const uint n_unknown)
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

	  /**
	    * Exact derivative.
	    */
	  Gradient exact_derivative(const Point& p, const uint n_unknown)
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


	  /**
	   * The required parameters for boundary conditions.
	   */
	  double right_hand_side(const Point& p, const uint n_unknown)
	  {
		  double f;
		  const double x=p(X), y=p(Y), z=p(Z);
		  const double pi=M_PI, mu=_mu, lambda=_lambda;

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

};

/**
 * Finding the errors
 */
class ExactValue: public FunctionBase<Number>
{
public:

	// constructor and destructor
	ExactValue(TestCase &test_case):
		_test_case(&test_case)
	{
		this->_initialized=true;
	}
	~ExactValue(){}


	// libmesh musts
	UniquePtr<FunctionBase<Number> > clone() const libmesh_override
	{
		return UniquePtr<FunctionBase<Number> >  (new ExactValue(*_test_case));
	}

	Number operator() (const Point & p, const Real = 0) libmesh_override
			{
		return _test_case->exact_solution(p,0);
			}

	virtual void operator() (const Point & p, const Real, DenseVector<Number> & output) libmesh_override
			{
		unsigned int size = output.size();
		for (unsigned int i=0; i != size; ++i)
			output(i) = _test_case->exact_solution(p,i);
			}

private:
	TestCase *_test_case;
};

class ExactDeriv: public FunctionBase<Gradient>
{
public:

	// constructor and destructor
	ExactDeriv(TestCase &test_case):
		_test_case(&test_case)
	{
		this->_initialized=true;
	}
	~ExactDeriv(){}


	// libmesh musts
	UniquePtr<FunctionBase<Gradient> > clone() const libmesh_override
	{
		return UniquePtr<FunctionBase<Gradient> >  (new ExactDeriv(*_test_case));
	}

	Gradient operator() (const Point & p, const Real = 0) libmesh_override
			{
		return _test_case->exact_derivative(p,0);
			}

	virtual void operator() (const Point & p, const Real, DenseVector<Gradient> & output) libmesh_override
			{
		unsigned int size = output.size();
		for (unsigned int i=0; i != size; ++i)
			output(i) = _test_case->exact_derivative(p,i);
			}

private:
	TestCase *_test_case;
};
/****************************************************************************
 *                             Assembler                                    *
 ****************************************************************************/

class LinearElasticity : public System::Assembly
{
private:
  EquationSystems &es;
  double _mu, _lambda;

public:

  LinearElasticity (EquationSystems &es_in) :
    es(es_in),
    _mu(0),
    _lambda(0)
  {
	    // Hard code material parameters for the sake of simplicity
	    const Real poisson_ratio = 0.3;
	    const Real young_modulus = 1.;

	    // Define the Lame constants
	    _lambda = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
	    _mu = young_modulus/(2.*(1.+poisson_ratio));
  }

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(
    unsigned int i,
    unsigned int j)
  {
    return i == j ? 1. : 0.;
  }

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(
    unsigned int i,
    unsigned int j,
    unsigned int k,
    unsigned int l)
  {

    return _lambda * kronecker_delta(i,k) * kronecker_delta(j,l) +
           _mu * (kronecker_delta(i,j) * kronecker_delta(k,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
  }

  /**
    * Exact solution.
    */
  double exact_solution(const Point& p, const uint n_unknown)
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

  /**
    * Exact derivative.
    */
  Gradient exact_derivative(const Point& p, const uint n_unknown)
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


  /**
   * The required parameters for boundary conditions.
   */
  double right_hand_side(const Point& p, const uint n_unknown)
  {
	  double f;
	  const double x=p(X), y=p(Y), z=p(Z);
	  const double pi=M_PI, mu=_mu, lambda=_lambda;

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

  void evaluate_errors()
  {
	  TestCase tc(_lambda, _mu);
	  ExactValue eval(tc);
	  ExactDeriv eder(tc);

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

  /**
   * The required parameters for boundary conditions.
   */
  void boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q)
  {
	  double penalty = 1.e10;
	  q = penalty;
	  t = penalty*exact_solution(p,n_unknown);
  }


  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble()
  {
	  const MeshBase& mesh = es.get_mesh();

	  const unsigned int dim = mesh.mesh_dimension();

	  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");

	  const unsigned int u_var = system.variable_number ("u");

	  // Dof Map
	  const DofMap& dof_map = system.get_dof_map();
	  FEType fe_type = dof_map.variable_type(u_var);
	  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
	  QGauss qrule (dim, fe_type.default_quadrature_order());
	  fe->attach_quadrature_rule (&qrule);

	  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	  QGauss qface(dim-1, fe_type.default_quadrature_order());
	  fe_face->attach_quadrature_rule (&qface);

	  const std::vector<Real>& JxW = fe->get_JxW();
	  const std::vector<std::vector<Real> >& phi = fe->get_phi();
	  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	  const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
	  const std::vector<Real>& JxW_face = fe_face->get_JxW();

	  const std::vector<Point>& q_point = fe->get_xyz();
	  const std::vector<Point >& qface_point = fe_face->get_xyz();

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

	  std::vector<dof_id_type> dof_indices;
	  std::vector< std::vector<dof_id_type> > dof_indices_var(3);

	  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	  for ( ; el != end_el; ++el)
	  {
		  const Elem* elem = *el;

		  dof_map.dof_indices (elem, dof_indices);
		  for(unsigned int var=0; var<3; var++)
		  {
			  dof_map.dof_indices (elem, dof_indices_var[var], var);
		  }

		  const unsigned int n_dofs   = dof_indices.size();
		  const unsigned int n_var_dofs = dof_indices_var[0].size();

		  fe->reinit (elem);

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
        	f_vec(0) = right_hand_side(q_point[qp],0);
        	f_vec(1) = right_hand_side(q_point[qp],1);
        	f_vec(2) = right_hand_side(q_point[qp],2);
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
        				boundary_penalty(qface_point[qp],id,i,t_vec(i),q_vec(i));
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
      }
  }

};


// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     40,
                                     40,
                                     40,
                                     0., 1.,
                                     0., 1.,
                                     0., 1.,
                                     HEX27);


  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem& system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Attach a SolverConfiguration object to system.linear_solver
  PetscLinearSolver<Number>* petsc_linear_solver =
    libmesh_cast_ptr<PetscLinearSolver<Number>*>(system.get_linear_solver());
  libmesh_assert(petsc_linear_solver);
  PetscSolverConfiguration petsc_solver_config(*petsc_linear_solver);
  petsc_linear_solver->set_solver_configuration(petsc_solver_config);

  // Add three displacement variables, u and v, to the system
  system.add_variable("u", FIRST, HIERARCHIC);
  system.add_variable("v", FIRST, HIERARCHIC);
  system.add_variable("w", FIRST, HIERARCHIC);

  LinearElasticity le(equation_systems);
  system.attach_assemble_object(le);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Print the errors
  le.evaluate_errors();

  // Plot the solution
  // Use single precision in this case (reduces the size of the exodus file)
  ExodusII_IO exo_io(mesh, /*single_precision=*/true);

  // First plot the displacement field using a nodal plot
  std::set<std::string> system_names;
  system_names.insert("Elasticity");
  exo_io.write_equation_systems("displacement_and_stress.exo",equation_systems,&system_names);

  // All done.
  return 0;
}
