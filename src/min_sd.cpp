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

#define EPS       1.0e-6

/* ---------------------------------------------------------------------- */

MinSD::MinSD(LAMMPS *lmp) : MinCG(lmp) {}

/* ----------------------------------------------------------------------
   minimization via steepest descent
------------------------------------------------------------------------- */

void MinSD::iterate(int n)
{
  int i,fail;
  double alpha,dot,dotall;
  double *f;

  f = atom->f[0];
  for (int i = 0; i < ndof; i++) h[i] = f[i];

  neval = 0;

  for (niter = 0; niter < n; niter++) {

    update->ntimestep++;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    fail = (this->*linemin)(ndof,atom->x[0],h,ecurrent,dmin,dmax,alpha,neval);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently, all done

    if (neval >= update->max_eval) break;

    if (fail || fabs(ecurrent-eprevious) <= 
    	update->tolerance * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS))
      break;

    // set h to new f = -Grad(x)
    // done if size sq of grad vector < EPS

    f = atom->f[0];
    dot = 0.0;
    for (i = 0; i < ndof; i++) dot += f[i]*f[i];
    MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
    if (dotall < EPS) break;

    for (i = 0; i < ndof; i++) h[i] = f[i];

    // output for thermo, dump, restart files

    if (output->next == update->ntimestep) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}
