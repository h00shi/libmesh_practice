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




Number exact_solution
(const Point& p,
const Parameters&,
const std::string& system_name,
const std::string& unk_name)
{
	libmesh_assert_equal_to(system_name, "Poisson");

	if(unk_name == "u")
	{
		const Real x = p(0);
		const Real y = p(1);
		const Real z = 0;
		return cos(.5*M_PI*x)*sin(.5*M_PI*y)*cos(.5*M_PI*z);
	}
	else
	{
		libmesh_assert_msg(0, "Variable " << unk_name << " has no exact solution.");
		return 0;
	}
}

Gradient exact_deriv
(const Point& p,
const Parameters& ,
const std::string& system_name,
const std::string& unk_name)
{
	libmesh_assert_equal_to(system_name, "Poisson");

	Gradient g(0);

	if(unk_name == "u")
	{
		const Real x = p(0);
		const Real y = p(1);
		const Real z = 0;

		g(0) = -0.5*M_PI*sin(.5*M_PI*x)*sin(.5*M_PI*y)*cos(.5*M_PI*z);
		g(1) = 0.5*M_PI *cos(.5*M_PI*x)*cos(.5*M_PI*y)*cos(.5*M_PI*z);
		g(2) = -0.5*M_PI*cos(.5*M_PI*x)*sin(.5*M_PI*y)*sin(.5*M_PI*z);
	}
	else
	{
		libmesh_assert_msg(0, "Variable " << unk_name << " has no exact derivative.");
	}

	return g;
}

Number right_hand_side_interior
(const Point& p,
const Parameters& params,
const std::string& system_name,
const std::string& unk_name)
{
	libmesh_assert_equal_to(system_name, "Poisson");

	if(unk_name == "phi")
	{
		return 1;
	}
	else if(unk_name =="u")
	{
		const Real eps = 1.e-6;
		const double u_u = exact_solution(Point(p(0),p(1)+eps), params, system_name, unk_name);
		const double u_d = exact_solution(Point(p(0),p(1)-eps), params, system_name, unk_name);
		const double u_r = exact_solution(Point(p(0)+eps,p(1)), params, system_name, unk_name);
		const double u_l = exact_solution(Point(p(0)-eps,p(1)), params, system_name, unk_name);
		const double u_c = exact_solution(Point(p(0),p(1)), params, system_name, unk_name);
		return -( u_d + u_u + u_l + u_r - 4. * u_c ) / (eps * eps) ;
    }
	else
	{
		libmesh_assert_msg(0, "Variable name" << unk_name << "not supported");
		return 0;
	}
}



Gradient right_hand_side_bdry
(const Point& p,
const Parameters& params,
const std::string& system_name,
const std::string& unk_name,
const int id)
{
	libmesh_assert_equal_to(system_name, "Poisson");
	const double penalty = 1.e10;
	Gradient g;

	if(unk_name == "phi")
	{
		switch(id)
		{
		case 1:
			g(0) = penalty;
			g(1) = 0;
			break;
		case 5:
			g(0) = 0;
			g(1) = 0;
			break;
		default:
			libmesh_assert_msg(0, "boundary id " << id << " not supported. ");
		}
	}
	else if(unk_name =="u")
	{
		g(0) = penalty;
		g(1) = penalty * exact_solution(p, params, system_name, unk_name);
	}
	else
	{
		libmesh_assert_msg(0, "Variable name" << unk_name << "not supported");
	}

	return g;
}
