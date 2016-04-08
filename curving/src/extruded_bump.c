
static char help[] = "Newton's method for a two-variable system, sequential.\n\n";

/**
 * Here we have a surface with the parameterization:
 * (x,y,z) = r(u,v)
 * n(u,v) = r_u \times r_v / \| r_u \times r_v \|
 *
 * We wish to project a point xs, ys, zs into this surface, i.e.,
 * solve for u, v and d such that:
 * xs - x(u,v) - d n_x(u,v) = 0
 * ys - y(u,v) - d n_y(u,v) = 0
 * zs - z(u,v) - d n_z(u,v) = 0
 *
 * We will use PETSc's SNES framework.
 */
#include "petscsnes.h"
#include <math.h>

/*
   User-defined routines
*/

#define n_unknown 3

struct Geometry
{
	// Geometry of the problem
	double h,c,r;

	// The coordinates in the question
	double xs,ys,zs;   //point to be projected
	double x, y, z;    //numerical solution
	double xe, ye, ze; //exact solution

	// The generalized coordinates
	double u0, v0, d0; //initial guess
	double u, v, d;    //numerical solution
	double ue, ve, de; //exact solution
	
};
typedef struct Geometry Geometry;

extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode SetUpGeometry(Geometry*);
extern PetscErrorCode CheckGeometry(Vec, Geometry*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  SNES           snes;         /* nonlinear solver context */
  KSP            ksp;          /* linear solver context */
  PC             pc;           /* preconditioner context */
  Vec            x,r;          /* solution, residual vectors */
  Mat            J;            /* Jacobian matrix */
  PetscErrorCode ierr;
  PetscInt       its;
  PetscMPIInt    size;
  PetscScalar    *xx;
  int            i;
  Geometry       geom;

  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Example is only for sequential runs");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the geometry
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetUpGeometry(&geom);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix and vector data structures; set corresponding routines
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create vectors for solution and nonlinear function
  */
  ierr = VecCreateSeq(PETSC_COMM_SELF, n_unknown, &x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

  /*
     Create Jacobian matrix data structure
  */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, n_unknown, n_unknown, NULL, &J);CHKERRQ(ierr);

  {
	  int idx[20]; double ones[100];
	  for (i=0; i<n_unknown; i++) idx[i]=i;
	  for (i=0; i<n_unknown*n_unknown; i++) ones[i]=1;
	  ierr=MatSetValues(J, n_unknown, idx, n_unknown, idx, ones, INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  

  // set the function
  ierr = SNESSetFunction(snes,r,FormFunction, &geom);CHKERRQ(ierr);

  // set the snes function
  ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault, &geom);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Set linear solver defaults for this problem. By extracting the
     KSP and PC contexts from the SNES context, we can then
     directly call any KSP and PC routines to set various options.
  */
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPPREONLY);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);

  /*
     Set SNES/KSP/KSP/PC runtime options, e.g.,
         -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
     These options will override those specified above as long as
     SNESSetFromOptions() is called _after_ any other customization
     routines.
  */
  ierr=SNESSetTolerances(snes,1e-16,1e-20,1e-20,100,1000);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr  = VecGetArray(x,&xx);CHKERRQ(ierr);
    xx[0] = geom.u0;
	xx[1] = geom.v0;
	xx[2] = geom.d0;
    ierr  = VecRestoreArray(x,&xx);CHKERRQ(ierr);
	
  /*
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */

  ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  CheckGeometry(x, &geom);
  
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of SNES iterations = %D\n",its);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr); ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/*
   FormFunction1 - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x    - input vector
.  ctx  - optional user-defined context

   Output Parameter:
.  f - function vector
 */
PetscErrorCode FormFunction(SNES snes,Vec xvec,Vec f,void *ctx)
{
  PetscErrorCode    ierr;
  const PetscScalar *xx;
  PetscScalar       *ff;

  /*
   Get pointers to vector data.
      - For default PETSc vectors, VecGetArray() returns a pointer to
        the data array.  Otherwise, the routine is implementation dependent.
      - You MUST call VecRestoreArray() when you no longer need access to
        the array.
   */
  ierr = VecGetArrayRead(xvec,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

  /* extract data */
  const Geometry *g = (Geometry *) ctx;

  // xs, ys and zs
  const double xs = g->xs;
  const double ys = g->ys;
  const double zs = g->zs;

  // u, v and d
  const double u = xx[0];
  const double v = xx[1];
  const double d = xx[2];

  // r, h and c
  const double r = g->r;
  const double h = g->h;
  const double c = g->c;

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
  ierr = VecRestoreArrayRead(xvec,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SetUpGeometry(Geometry* g)
{
	// set the geometry
	g->c = 1;
	g->r = 1;
	g->h = 0.5;

	// set the point for which surface friend has to be found
	// lets create a manufactured solution
	g->ue = 0.5;
	g->ve = 0.7 * M_PI/3;
	g->de = 0.25;
		
	const double denom_n = sqrt(g->c * g->c + sin(g->ve) * sin(g->ve) );
	const double n_x = g->c * sin(g->ve) / denom_n;
	const double n_y = -sin(g->ve) / denom_n;
	const double n_z = g->c * cos(g->ve) / denom_n;
	
    // Exact solution
	g->xe = g->ue + g->r * sin(g->ve);
	g->ye = g->c * g->ue;
	g->ze = g->r * cos(g->ve) - g->h;

	// The point to be projected
	g->xs= g->xe + g->de*n_x;
	g->ys= g->ye + g->de*n_y;
	g->zs= g->ze + g->de*n_z;
	  
	// Initial guess
	g->u0= g->ys / g->c;
	g->v0= asin( (g->xs - g->u0) / g->r );
	g->d0 = 0;

	// set the rest to zero
	g->x = g->y = g->z = 0;
	g->u = g->v = g->d = 0;
	
	return 0;
}

PetscErrorCode CheckGeometry(Vec xvec, Geometry* g)
{
  PetscErrorCode    ierr;
  const PetscScalar *xx;

  // get vector
  ierr = VecGetArrayRead(xvec,&xx);CHKERRQ(ierr);

  // set the answer
  g->u = xx[0];
  g->v = xx[1];
  g->d = xx[2];

  // set the new coordinates
  g->x = g->u + g->r * sin(g->v);
  g->y = g->c * g->u;
  g->z = g->r * cos(g->v) - g->h;

  // check the difference
  printf(" init_x: %10.8lf %10.8lf %10.8lf\n", g->xs, g->ys, g->zs);  
  printf("exact_x: %10.8lf %10.8lf %10.8lf\n", g->xe, g->ye, g->ze);
  printf("  num_x: %10.8lf %10.8lf %10.8lf\n", g->x , g->y , g->z );  

  printf(" init_u: %10.8lf %10.8lf %10.8lf\n", g->u0, g->v0, g->d0);  
  printf("exact_u: %10.8lf %10.8lf %10.8lf\n", g->ue, g->ve, g->de);
  printf("  num_u: %10.8lf %10.8lf %10.8lf\n", g->u , g->v , g->d );  

  // restore vector
  ierr = VecRestoreArrayRead(xvec,&xx);CHKERRQ(ierr);

  return 0;	
}



