

#ifndef TEST_CASE_HXX
#define TEST_CASE_HXX

/****************************************************************************
 *                                 Includes                                 *
 ****************************************************************************/
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>

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

//Petsc stuff
#include "petscsnes.h"
#include "petscsys.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


/****************************************************************************
 *                                    Enum                                  *
 ****************************************************************************/
enum Variables{U=0, V, W};
enum Coordinates{X=0, Y, Z};


/****************************************************************************
 *                                 Material                                 *
 ****************************************************************************/

class Material
{
	Material(){}
public:
	double lambda, mu;

	Material(const double lambda_val, const double mu_val):
		lambda(lambda_val),
		mu(mu_val)
	{}
};

/****************************************************************************
 *                                Test Case                                *
 ****************************************************************************/

/**
 * A class containing data for each test case (rhs, bdry condition).
 */
class TestCase
{

protected:
	  Material &_material;
	  bool _has_exact_solution;
	  bool _has_nonzero_rhs;

	  // constructor
	  TestCase(Material &material, bool has_exact_solution, bool has_nonzero_rhs):
		  _material(material),
		  _has_exact_solution(has_exact_solution),
		  _has_nonzero_rhs(has_nonzero_rhs)
		{}

public:

	  // destructor
	  virtual ~TestCase(){}

	  /**
	    * Exact solution if available.
	    */
	  virtual double exact_solution(const Point& , const uint )
	  { libmesh_assert(0); return 0;}

	  /**
	    * Exact derivative if available.
	    */
	  virtual Gradient exact_derivative(const Point& , const uint )
	  {libmesh_assert(0); Gradient g; return g;}

	  /**
	   * Right hand side of the equations.
	   */
	  virtual double right_hand_side(const Point& , const uint )
	  {libmesh_assert(0); return 0;}

	  /**
	   * Dirichlet boundary conditions parameters.
	   * Must always be overriden, so pure virtual.
	   */
	  virtual void boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q) = 0;

	  /**
	   * Query functions
	   */
	  bool has_exact_solution()
	  {return _has_exact_solution;}

	  bool has_nonzero_rhs()
	  {return _has_nonzero_rhs;}
};

/**
 * Manufactured solution test case.
 */
class TestCaseMMS: public TestCase
{
public:

	  // constructor
	  TestCaseMMS(Material &material):
		  TestCase(material, true, true)
		{}

	  // destructor
	  ~TestCaseMMS(){}

	  double exact_solution(const Point& p, const uint n_unknown);
	  Gradient exact_derivative(const Point& p, const uint n_unknown);
	  double right_hand_side(const Point& p, const uint n_unknown);
	  void boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q);
};

/**
 * Curving the sphere test case.
 */
class TestCaseSphere: public TestCase
{
private:
	double _r_in, _r_out;
public:

	  // constructor
	  TestCaseSphere(Material &material, const double r_in, const double r_out):
		  TestCase(material, false, false),
		  _r_in(r_in),
		  _r_out(r_out)
		{}

	  // destructor
	  ~TestCaseSphere(){}

	  // Boundary surfaces should be numbered like this:
	  // 1-> perpendicular to x
	  // 2-> perpendicular to y
	  // 3-> perpendicular to z
	  // 4-> inner surface
	  // 5-> outer surface
	  // 6-> upper surface
	  void boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q);
};

/**
 * Extruded bump test case.
 */
class TestCaseBump: public TestCase
{
private:
	// geometry parameters
	double _r, _h, _c, _l, _b;

	// petsc stuff for projecting
	SNES _snes;
	Vec _vec_x, _vec_r;
	Mat _mat_j;

	double _xs, _ys, _zs;

	// avoid projecting one point twice
	const Point *_last_projectee;
	Point _last_projection;

	// parallel stuff
	Parallel::Communicator& _comm;

public:

	  // constructor
	  TestCaseBump(Material &material, const double r, const double v_max, const double c, const double b, Parallel::Communicator& comm);

	  // destructor
	  ~TestCaseBump();

	  // x, y and z as functions of u and v
	  Point xyz(const double u, const double v) const
	  {
		  Point p;

		  const double r = _r - _b * u;
		  const double h = _h * (1 - _b / _r * u);
		  p(X) = u +  _l/2. + r * sin(v);
		  p(Y) = _c * u;
		  p(Z) = r * cos(v) - h;

		  return p;
	  }
	  // n as function of u and v
	  Point normal(const double, const double v) const
	  {
		  Point n;

		  const double rurv_x = _c*sin(v);
		  const double rurv_y = _b*(1-_h/_r*cos(v))-sin(v);
		  const double rurv_z = _c*cos(v);
		  const double denom = sqrt( _c*_c + rurv_y*rurv_y);

		  n(X) = rurv_x / denom;
		  n(Y) = rurv_y / denom;
		  n(Z) = rurv_z / denom;

		  return n;
	  }
	  // Initial Guess
	  void u0v0(const Point &p, double& u, double& v) const
	  {

		  u= p(Y) / _c;

		  const double r = _r - _b * u;

		  const double sin_candidate  = (p(X) - _l/2 - u) / r;
		  if(abs(sin_candidate) < 1)
		  {
			  v=  asin(sin_candidate);
		  }
		  else
		  {
			  const double v_max = acos(_h/_r);
			  v = v_max * .9;
		  }
	  }

	  // residual for projection on bump
	  static PetscErrorCode snes_residual(SNES snes,Vec xvec,Vec f,void *ctx);
	  // project on the bump surface
	  Point project_on_bump(const Point& p);
	  // test the projection code
	  void project_on_bump_test();

	  // Boundary surfaces should be numbered like this:
	  // 1 -> perpendicular to x
	  // 2 -> perpendicular to y
	  // 3 -> perpendicular to z
	  // 4 -> to be curved
	  void boundary_penalty(const Point& p, const uint id, const uint n_unknown, double &t, double &q);
};

/****************************************************************************
 *                    Exact solution function wrappers                      *
 ****************************************************************************/

/**
 * Function wrapper for computing the exact solution
 */
class ExactValue: public FunctionBase<Number>
{
private:
	TestCase &_test_case;

public:

	// constructor and destructor
	ExactValue(TestCase &test_case):
		_test_case(test_case)
	{
		libmesh_assert(_test_case.has_exact_solution());
		this->_initialized=true;
	}
	~ExactValue(){}


	// libmesh musts
	UniquePtr<FunctionBase<Number> > clone() const libmesh_override
	{
		return UniquePtr<FunctionBase<Number> >  (new ExactValue(_test_case));
	}

	Number operator() (const Point & p, const Real = 0) libmesh_override
			{
		return _test_case.exact_solution(p,0);
			}

	virtual void operator() (const Point & p, const Real, DenseVector<Number> & output) libmesh_override
			{
		unsigned int size = output.size();
		for (unsigned int i=0; i != size; ++i)
			output(i) = _test_case.exact_solution(p,i);
			}
};

/**
 * Function wrapper for computing the exact derivative
 */
class ExactDeriv: public FunctionBase<Gradient>
{
private:
	TestCase &_test_case;

public:

	// constructor and destructor
	ExactDeriv(TestCase &test_case):
		_test_case(test_case)
	{
		libmesh_assert(_test_case.has_exact_solution());
		this->_initialized=true;
	}
	~ExactDeriv(){}


	// libmesh musts
	UniquePtr<FunctionBase<Gradient> > clone() const libmesh_override
	{
		return UniquePtr<FunctionBase<Gradient> >  (new ExactDeriv(_test_case));
	}

	Gradient operator() (const Point & p, const Real = 0) libmesh_override
			{
		return _test_case.exact_derivative(p,0);
			}

	virtual void operator() (const Point & p, const Real, DenseVector<Gradient> & output) libmesh_override
			{
		unsigned int size = output.size();
		for (unsigned int i=0; i != size; ++i)
			output(i) = _test_case.exact_derivative(p,i);
			}
};

/****************************************************************************
 *                             Assembler                                    *
 ****************************************************************************/

class LinearElasticity : public System::Assembly
{
private:
  EquationSystems &es;
  Material &_material;
  TestCase &_test_case;

public:

  LinearElasticity (EquationSystems &es_in, Material &material, TestCase &tc) :
    es(es_in),
    _material(material),
    _test_case(tc)
  { }

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(const uint i, const uint j)
  {
    return i == j ? 1. : 0.;
  }

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(const uint i, const uint j, const uint k, const uint l)
  {

    return _material.lambda * kronecker_delta(i,k) * kronecker_delta(j,l) +
           _material.mu * (kronecker_delta(i,j) * kronecker_delta(k,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
  }

  /**
   * Find the errors if we have a manufactured solution.
   */
  void post_process(std::ostream &);


  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble();

 };

#endif /* TEST_CASE_HXX */
