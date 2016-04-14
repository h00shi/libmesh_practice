#include "linear_elasticity.hxx"
#include "libmesh/gmsh_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel.h"
#include "libmesh/metis_partitioner.h"

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

    ierr = KSPSetTolerances(_petsc_linear_solver.ksp(), PETSC_DEFAULT , 1e-11, 1e60 , 3000);
    libmesh_assert(ierr == 0);

    // override previous stuff if need be
    ierr = KSPSetFromOptions(_petsc_linear_solver.ksp());
    libmesh_assert(ierr == 0);

    // print solver options
    ierr =  KSPView(_petsc_linear_solver.ksp(), PETSC_VIEWER_STDOUT_WORLD);
    libmesh_assert(ierr == 0);

  }

  // The linear solver object that we are configuring
  PetscLinearSolver<Number>& _petsc_linear_solver;

};

/****************************************************************************
 *                                   MAIN                                   *
 ****************************************************************************/
// Begin the main program.
int main (int argc, char** argv)
{
  // Getpot
  GetPot getpot(argc, argv);

  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);

  // Read the mesh
  GmshIO gmsh(mesh);
  if(init.comm().rank() == 0)
  {
	  std::string mesh_name;
	  mesh_name = getpot("mesh_name", "../mesh/simple_box_5.msh");
	  gmsh.read(mesh_name);
  }
  // broadcast the mesh
  MeshCommunication mc;
  mc.broadcast(mesh);
  // paratition the mesh
  MetisPartitioner metis;
  metis.partition(mesh, init.comm().size());
  // initialize the mesh
  mesh.prepare_for_use();
  mesh.all_second_order(true);

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
  const int order= getpot("order", (int)SECOND);
  const int fe_family= getpot("fe_family",(int)HIERARCHIC);
  system.add_variable("u", (Order)order, (FEFamily)fe_family);
  system.add_variable("v", (Order)order, (FEFamily)fe_family);
  system.add_variable("w", (Order)order, (FEFamily)fe_family);

  // add the assembler
  // Hard code material parameters for the sake of simplicity
  const Real poisson_ratio = 0.3;
  const Real young_modulus = 1.;

  // Define the Lame constants
  const double lambda = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
  const double mu = young_modulus/(2.*(1.+poisson_ratio));

  Material material(lambda, mu);
  UniquePtr<TestCase> test_case;
  std::string test_case_str = getpot("test_case", "mms");
  if(test_case_str == "mms")
	  test_case.reset(new TestCaseMMS(material));
  else if(test_case_str == "sphere")
	  test_case.reset(new TestCaseSphere(material, 1, 2));
  else if(test_case_str == "bump")
	  test_case.reset(new TestCaseBump(material, 1 /*r*/, M_PI/6. /*vmax*/, 4 /*c*/, 2 /*b*/, init.comm()));
  else
	  libmesh_assert(0);

  LinearElasticity le(equation_systems, material, *test_case);
  system.attach_assemble_object(le);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system and print reason
  system.solve();
  petsc_linear_solver->print_converged_reason();

  // output file name
  std::string out_name = getpot("out_name", "displacement_and_stress");
  std::stringstream string_stream;

  /*
   *   post process
   */

  // open the file
  std::ofstream out_stream;
  if (init.comm().rank() == 0)
  {
	  string_stream.str("");
	  string_stream << out_name << ".pp";
	  out_stream.open(string_stream.str().c_str(), std::ios::out);
	  libmesh_assert(out_stream.good());
  }

  // post process and print in file
  le.post_process(out_stream);

  // close the file
  if (init.comm().rank() == 0)
  {
	  out_stream.close();
  }

  // Plot the solution
  // Use single precision in this case (reduces the size of the exodus file)
  ExodusII_IO exo_io(mesh, /*single_precision=*/true);

  std::set<std::string> system_names;
  system_names.insert("Elasticity");

  string_stream.str("");
  string_stream << out_name << ".exo";

  exo_io.write_equation_systems(string_stream.str(),equation_systems,&system_names);
  std::cout << string_stream.str()  << std::endl;

  // All done.
  return 0;
}
