#include "linear_elasticity.hxx"
#include "libmesh/gmsh_io.h"

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

    ierr = KSPSetTolerances(_petsc_linear_solver.ksp(), 1e-20, 1e-11, 1e30, 3000);
    libmesh_assert(ierr == 0);
  }

  // The linear solver object that we are configuring
  PetscLinearSolver<Number>& _petsc_linear_solver;

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

  // Read the mesh
  {
	  GmshIO gmsh(mesh);
	  gmsh.read("../mesh/sphere.msh");
	  mesh.all_second_order(true);
	  mesh.prepare_for_use();
  }

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
  system.add_variable("u", THIRD, HIERARCHIC);
  system.add_variable("v", THIRD, HIERARCHIC);
  system.add_variable("w", THIRD, HIERARCHIC);

  // add the assembler
  // Hard code material parameters for the sake of simplicity
  const Real poisson_ratio = 0.3;
  const Real young_modulus = 1.;

  // Define the Lame constants
  const double lambda = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
  const double mu = young_modulus/(2.*(1.+poisson_ratio));

  Material material(lambda, mu);
  //TestCaseMMS test_case(material);
  TestCaseSphere test_case(material, 1, 2);
  LinearElasticity le(equation_systems, material, test_case);
  system.attach_assemble_object(le);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Print the errors
  le.post_process();

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
