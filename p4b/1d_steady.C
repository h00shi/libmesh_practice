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

Real exact_solution_1d(const Real x, const Real t)
{
	const double tau = -0.1e1 + (exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * exp(-sqrt(0.2e1) * x / 0.2e1) * cos(sqrt(0.2e1) * x / 0.2e1) - 0.1e1 / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * (0.2e1 * cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) * exp(-sqrt(0.2e1) * x / 0.2e1) * sin(sqrt(0.2e1) * x / 0.2e1) + (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * exp(sqrt(0.2e1) * x / 0.2e1) * cos(sqrt(0.2e1) * x / 0.2e1) - 0.1e1 / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * (0.2e1 * cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) * exp(sqrt(0.2e1) * x / 0.2e1) * sin(sqrt(0.2e1) * x / 0.2e1);
	const double sig = (exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * exp(-sqrt(0.2e1) * x / 0.2e1) * sin(sqrt(0.2e1) * x / 0.2e1) + 0.1e1 / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * (0.2e1 * cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) * exp(-sqrt(0.2e1) * x / 0.2e1) * cos(sqrt(0.2e1) * x / 0.2e1) - (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * exp(sqrt(0.2e1) * x / 0.2e1) * sin(sqrt(0.2e1) * x / 0.2e1) - 0.1e1 / (pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) - 0.2e1 * pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + pow(exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) + 0.2e1 * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) + pow(sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1) * pow(exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1), 0.2e1)) * sin(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * (0.2e1 * cos(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) * exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(-sqrt(0.2e1) * 0.3141592654e1 / 0.2e1) - exp(sqrt(0.2e1) * 0.3141592654e1 / 0.2e1)) * exp(sqrt(0.2e1) * x / 0.2e1) * cos(sqrt(0.2e1) * x / 0.2e1);

	return tau * cos(t) + sig * sin(t);
}

//Assembly function
void assemble_heat(EquationSystems& es, const std::string& system_name);
void assemble_err(EquationSystems& es);
void assemble_steady(EquationSystems& es, const std::string& system_name);

void init_heat (EquationSystems& es, const std::string& system_name);

/****************************************************************************
 *                             Main                                         *
 ****************************************************************************/

int main (int argc, char** argv)
{

	// Initialize libraries, like in example 2.
	LibMeshInit init (argc, argv);
	GetPot getpot(argc, argv);

	// important data
	const uint n_timestep = getpot("NT", 50);
	const Real dt = getpot("dt", 0.025);
	const int n_elem = getpot("NX", 100);
	const double time_write = getpot("WT",10);
	std::string mesh_name = getpot("mesh","square");
	const bool exact_time = getpot("exact", false);
	const bool exo_io = getpot("exo",false);

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
				QUAD4);
		if(exact_time)
		{
			std::cerr << "No exact time" << std::endl;
			return 0;
		}
	}
	else
	{
		GmshIO gmsh(mesh);
		gmsh.read(mesh_name);
		mesh.all_first_order();
		mesh.prepare_for_use();
		if(exact_time)
		{
			std::cerr << "No exact time" << std::endl;
			return 0;
		}
	}
	mesh.print_info();

	// define the system for u
	EquationSystems equation_systems (mesh);

	// define system for transient problem
	TransientLinearImplicitSystem & system_heat =
			equation_systems.add_system<TransientLinearImplicitSystem> ("Heat");
	system_heat.add_variable("u",FIRST);
	system_heat.attach_assemble_function (assemble_heat);
	system_heat.attach_init_function (init_heat);

	// define the steady state solution
	ExplicitSystem &system_err = equation_systems.add_system<ExplicitSystem> ("Error");
	system_err.add_variable("e", FIRST);

	// define the steady state solution
	LinearImplicitSystem &system_steady =
	equation_systems.add_system<LinearImplicitSystem> ("Steady");
	system_steady.add_variable("sig",FIRST);
	system_steady.add_variable("tau",FIRST);
	system_steady.attach_assemble_function (assemble_steady);

	// init system
	equation_systems.init ();
	equation_systems.print_info();

	// solve the steady state
	system_steady.solve();

	// set global options
	equation_systems.parameters.set<Real>("time_write")= time_write;
	equation_systems.parameters.set<bool>("exact_time")= exact_time;


	//write the first time step
	std::string exodus_filename = "transient_ex1.e";
	if(exo_io)
	{
		std::set<std::string> system_names;
		system_names.insert("Error");
		system_names.insert("Heat");
		system_names.insert("Steady");
		std::string exodus_filename = "transient_ex1.e";
		ExodusII_IO(mesh).write_equation_systems (exodus_filename, equation_systems, &system_names);
	}

	/*
	 * start solving
	 */
	system_heat.time   = 0.;
	for (unsigned int t_step = 0; t_step < n_timestep; t_step++)
	{
		// update the time
		equation_systems.parameters.set<Real> ("time_prev") = system_heat.time;
		system_heat.time += dt;
		equation_systems.parameters.set<Real> ("time") = system_heat.time;
		equation_systems.parameters.set<Real> ("dt")   = dt;

		// A pretty update message
		std::cout << " Solving time step "
				<< std::setw(2)
		<< std::right
		<< t_step
		<< ", time="
		<< std::fixed
		<< std::setw(6)
		<< std::setprecision(3)
		<< std::setfill('0')
		<< std::left
		<< system_heat.time
		<<  "..."
		<< std::endl;

		// update the old solution
		*system_heat.old_local_solution = *system_heat.current_local_solution;

		// solve
		system_heat.solve();
		assemble_err(equation_systems);

		// Output evey 10 timesteps to file.
		if ((t_step+1)%10 == 0)
		{
			if(exo_io)
			{
				ExodusII_IO exo(mesh);
				exo.append(true);
				exo.write_timestep(exodus_filename, equation_systems, t_step+1, system_heat.time);
			}
			else
			{
				VTKIO vtk(mesh);
				char fname[300];
				sprintf(fname, "out/transient_%d.pvtu",(t_step+1)/10);
				std::string fname2(fname);
				vtk.write_equation_systems(fname2,equation_systems);
			}
		}
	}

	  // All done.
	  return 0;
}

/****************************************************************************
 *                          Assemble Heat                                   *
 ****************************************************************************/
Number exact_value_heat (const Point& p,const Parameters& param,
                    const std::string&, const std::string&)
{
	Number ans = 0;
 	ans =  exact_solution_1d(p(0), param.get<Real>("time"));
	return ans;
}

Number init_value_heat (const Point&,const Parameters& ,
                    const std::string&, const std::string&)
{
	return 0;
}

void assemble_heat(EquationSystems& es, const std::string& system_name)
{

	libmesh_assert_equal_to (system_name, "Heat"); //check system
	const MeshBase& mesh = es.get_mesh();  //mesh
	const unsigned int dim = mesh.mesh_dimension(); //dimension
	TransientLinearImplicitSystem & system =
			es.get_system<TransientLinearImplicitSystem> ("Heat"); //system
	FEType fe_type = system.variable_type(0); // finite element type

	// finite element objects
	UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));
	UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

	// DofMap
	const DofMap& dof_map = system.get_dof_map();
	// Dof indices
	std::vector<dof_id_type> dof_indices;

	// quadrature objects
	QGauss qrule (dim,   fe_type.default_quadrature_order());
	QGauss qface (dim-1, fe_type.default_quadrature_order());
	fe->attach_quadrature_rule      (&qrule);
	fe_face->attach_quadrature_rule (&qface);

	// Jacobian times quadrature weight
	const std::vector<Real>& JxW      = fe->get_JxW();
	const std::vector<Real>& JxW_face = fe_face->get_JxW();

	// Value of shape function and its derivatives
	const std::vector<std::vector<Real> >& phi = fe->get_phi();
	const std::vector<std::vector<Real> >& psi = fe_face->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

	// location of quadrature points
	//const std::vector<Point>& qface_points = fe_face->get_xyz();

	// Local LHS ans RHS
	DenseMatrix<Number> Ke;
	DenseVector<Number> Fe;

	// time step
	const Real dt = es.parameters.get<Real>   ("dt");
	const Real t = es.parameters.get<Real>   ("time");
	const Real to = es.parameters.get<Real>   ("time_prev");

	/*
	 * Loop over all elements
	 */
	MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.local_elements_end();

	for ( ; el != end_el; ++el)
	{
		// pointer
		const Elem* elem = *el;

		// element related data
		dof_map.dof_indices (elem, dof_indices);
		fe->reinit (elem);

		// Zero LHS and RHS
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fe.resize (dof_indices.size());

		/*
		 * Loop over element quadrature points
		 */
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{
			// Values to hold the old solution & its gradient.
			Number   u_old = 0.;
			Number st=sin(t), sto=sin(to);
			Gradient grad_u_old;

			// Compute the old solution & its gradient.
			for (unsigned int l=0; l<phi.size(); l++)
			{
				u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
				grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
			}

			// Now compute the element matrix and RHS contributions.
			for (unsigned int i=0; i<phi.size(); i++)
			{
				Fe(i) += JxW[qp]* (
						u_old*phi[i][qp]
						- .5*dt*( grad_u_old*dphi[i][qp] )
						+ .5*dt*phi[i][qp]*(st+sto)
				);

				for (unsigned int j=0; j<phi.size(); j++)
				{
					// The matrix contribution
					Ke(i,j) += JxW[qp]*(
							phi[i][qp]*phi[j][qp]
							+ .5*dt*( dphi[i][qp]*dphi[j][qp])
					);
				}
			}
		}

		/*
		 * Loop over sides
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
						Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];

					// Matrix contribution
					for (unsigned int i=0; i<psi.size(); i++)
						for (unsigned int j=0; j<psi.size(); j++)
							Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
				}
			}
		}

		// add contribution to the system
		system.matrix->add_matrix (Ke, dof_indices);
		system.rhs->add_vector    (Fe, dof_indices);
	}

	// All done.
}

void init_heat (EquationSystems& es, const std::string& system_name)
{
	// check system name
	libmesh_assert_equal_to (system_name, "Heat");

	// reference
	TransientLinearImplicitSystem & system =
			es.get_system<TransientLinearImplicitSystem>("Heat");

	// set time
	es.parameters.set<Real> ("time") = system.time = 0;
	es.parameters.set<Real> ("time_prev") = 0;

	// project solution
	system.project_solution(init_value_heat, NULL, es.parameters);
}


/****************************************************************************
 *                          Assemble Error                                  *
 ****************************************************************************/

void assemble_err(EquationSystems& es)
{
	// mesh
	const MeshBase& mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	// systems
	TransientLinearImplicitSystem & system_heat =
			es.get_system<TransientLinearImplicitSystem>("Heat");
	ExplicitSystem & system_err =
			es.get_system<ExplicitSystem> ("Error");
	LinearImplicitSystem & system_ss =
			es.get_system<LinearImplicitSystem> ("Steady");

	// indexing
	const DofMap& dof_map = system_heat.get_dof_map();
	const DofMap& err_dof_map = system_err.get_dof_map();
	const DofMap& ss_dof_map = system_ss.get_dof_map();
	std::vector<dof_id_type> dof_indices;
	std::vector<dof_id_type> err_dof_indices;
	std::vector<dof_id_type> ss_dof_indices;
	std::vector<dof_id_type> sig_dof_indices;
	std::vector<dof_id_type> tau_dof_indices;

	//numbers
	const uint sig_var = system_ss.variable_number ("sig");
	const uint tau_var = system_ss.variable_number ("tau");

	// finite element for quadrature
	FEType fe_type = dof_map.variable_type(0);
	UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
	QGauss qrule (dim, SIXTH);
	fe->attach_quadrature_rule (&qrule);

	//get time and if we are exact time solution
	double t = es.parameters.get<Real>("time");
	const bool exact_time = es.parameters.get<bool>("exact_time");

	//project the steady state solution - only for 1d exact space exact time
	if(exact_time) system_err.project_solution(exact_value_heat, NULL, es.parameters);

	// find the error
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double error = 0;
	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;

		// update element data and dof numbers
		fe->reinit (elem);
		dof_map.dof_indices (elem, dof_indices);
		err_dof_map.dof_indices (elem, err_dof_indices);
		ss_dof_map.dof_indices (elem, ss_dof_indices);
		ss_dof_map.dof_indices (elem, sig_dof_indices,sig_var);
		ss_dof_map.dof_indices (elem, tau_dof_indices,tau_var);

		//get solution at dof's
		std::vector<Number> sig,tau,u;
		system_ss.solution->get(sig_dof_indices, sig);
		system_ss.solution->get(tau_dof_indices, tau);
		system_heat.solution->get(dof_indices, u);

		// for dof's
		for (uint i = 0 ; i < dof_indices.size() ; i++)
		{

			// state solution at dof
			double steady_solution;

			// exact time and space - only 1d case
			if(exact_time)
			{
				steady_solution = (*system_err.solution)(err_dof_indices[i]);
			}
			// exact time, numerical space
			else
			{
				steady_solution = tau[i]*cos(t) + sig[i]*sin(t);
				// set the error at this dof
				system_err.solution->set(err_dof_indices[i], steady_solution);
			}

			// update the l2 norm
			error += elem->volume()/dof_indices.size() * pow( steady_solution - u[i] , 2);
		}
	}
	error = sqrt(error);

	system_err.solution->close();
	system_err.update();

	if(es.parameters.get<Real>("time") > es.parameters.get<Real>("time_write"))
		fprintf(stderr, "%e %e \n", es.parameters.get<Real>("time"), error );

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
	QGauss qrule (dim,   fe_type.default_quadrature_order());
	QGauss qface (dim-1, fe_type.default_quadrature_order());
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

