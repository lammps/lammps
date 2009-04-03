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
#include "mpi.h"
#include "min_sd.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinSD::MinSD(LAMMPS *lmp) : Min(lmp) {}

/* ----------------------------------------------------------------------
   minimization via steepest descent
------------------------------------------------------------------------- */

int MinSD::iterate(int n)
{
  int i,fail,ntimestep;
  double dot,dotall;

  double *x = NULL;
  double *f = NULL;

  if (ndof) f = atom->f[0];
  for (i = 0; i < ndof; i++) h[i] = f[i];
  if (nextra)
    for (i = 0; i < nextra; i++)
      hextra[i] = fextra[i];

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
    dot = 0.0;
    for (i = 0; i < ndof; i++) dot += f[i]*f[i];
    MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
    if (nextra)
      for (i = 0; i < nextra; i++)
	dotall += fextra[i]*fextra[i];

    if (dotall < update->ftol * update->ftol) return FTOL;

    // set new search direction h to f = -Grad(x)

    for (i = 0; i < ndof; i++) h[i] = f[i];
    if (nextra)
      for (i = 0; i < nextra; i++)
	hextra[i] = fextra[i];

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
  
  return MAXITER;
}
