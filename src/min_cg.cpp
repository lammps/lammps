/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "mpi.h"
#include "min_cg.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

MinCG::MinCG(LAMMPS *lmp) : Min(lmp) {}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

int MinCG::iterate(int n)
{
  int i,fail,ntimestep;
  double beta,gg,dot[2],dotall[2];

  double *x = NULL;
  double *f = NULL;

  int ndoftotal;
  MPI_Allreduce(&ndof,&ndoftotal,1,MPI_INT,MPI_SUM,world);
  ndoftotal += nextra;

  if (ndof) f = atom->f[0];
  for (i = 0; i < ndof; i++) h[i] = g[i] = f[i];
  if (nextra)
    for (i = 0; i < nextra; i++)
      hextra[i] = gextra[i] = fextra[i];
  
  dot[0] = 0.0;
  for (i = 0; i < ndof; i++) dot[0] += f[i]*f[i];
  MPI_Allreduce(dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);
  if (nextra)
    for (i = 0; i < nextra; i++) gg += fextra[i]*fextra[i];

  neval = 0;

  for (niter = 0; niter < n; niter++) {

    ntimestep = ++update->ntimestep;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    if (ndof) x = atom->x[0];
    fail = (this->*linemin)(ndof,x,h,x0,ecurrent,dmax,alpha_final,neval);
    if (fail) return fail;

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) < 
	update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    if (ndof) f = atom->f[0];
    dot[0] = dot[1] = 0.0;
    for (i = 0; i < ndof; i++) {
      dot[0] += f[i]*f[i];
      dot[1] += f[i]*g[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    if (nextra)
      for (i = 0; i < nextra; i++) {
	dotall[0] += fextra[i]*fextra[i];
	dotall[1] += fextra[i]*gextra[i];
      }

    if (dotall[0] < update->ftol * update->ftol) return FTOL;

    // update new search direction h from new f = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be Fletcher-Reeves
    // reinitialize CG every ndof iterations by setting beta = 0.0

    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    if ((niter+1) % ndoftotal == 0) beta = 0.0;
    gg = dotall[0];

    for (i = 0; i < ndof; i++) {
      g[i] = f[i];
      h[i] = g[i] + beta*h[i];
    }
    if (nextra)
      for (i = 0; i < nextra; i++) {
	gextra[i] = fextra[i];
	hextra[i] = gextra[i] + beta*hextra[i];
      }

    // reinitialize CG if new search direction h is not downhill

    dot[0] = 0.0;
    for (i = 0; i < ndof; i++) dot[0] += g[i]*h[i];
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);

    if (nextra)
      for (i = 0; i < nextra; i++)
	dotall[0] += gextra[i]*hextra[i];

    if (dotall[0] <= 0.0) {
      for (i = 0; i < ndof; i++) h[i] = g[i];
      if (nextra)
	for (i = 0; i < nextra; i++)
	  hextra[i] = gextra[i];
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  return MAXITER;
}

