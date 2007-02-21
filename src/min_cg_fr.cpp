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
#include "min_cg_fr.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define EPS       1.0e-6

/* ---------------------------------------------------------------------- */

MinCGFR::MinCGFR(LAMMPS *lmp) : MinCG(lmp) {}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
   Fletcher-Reeves formulation
------------------------------------------------------------------------- */

void MinCGFR::iterate(int n)
{
  int i,gradsearch,fail;
  double alpha,beta,gg,dot,dotall;

  for (int i = 0; i < ndof; i++) h[i] = g[i] = f[i];
  for (i = 0; i < ndof_extra; i++) hextra[i] = gextra[i] = fextra[i];

  dot = 0.0;
  for (i = 0; i < ndof; i++) dot += f[i]*f[i];
  MPI_Allreduce(&dot,&gg,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < ndof_extra; i++) gg += fextra[i]*fextra[i];

  neval = 0;
  gradsearch = 1;

  for (niter = 0; niter < n; niter++) {

    update->ntimestep++;

    // line minimization along direction h from current atom->x

    eprevious = energy;
    fail = (this->*linemin)(neval);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently:
    //   all done if searched in grad direction
    //   else force next search to be in grad direction (CG restart)

    if (neval >= update->max_eval) break;

    if (fail || fabs(energy-eprevious) <= 
    	update->tolerance * 0.5*(fabs(energy) + fabs(eprevious) + EPS)) {
      if (gradsearch == 1) break;
      gradsearch = -1;
    }

    // update h from new f = -Grad(x) and old g
    // old g,h must have migrated with atoms to do this correctly
    // done if size sq of grad vector < EPS
    // force new search dir to be grad dir if need to restart CG
    // set gradsesarch to 1 if will search in grad dir on next iteration

    dot = 0.0;
    for (i = 0; i < ndof; i++) dot += f[i]*f[i];
    MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < ndof_extra; i++) dotall += fextra[i]*fextra[i];

    beta = dotall/gg;
    gg = dotall;
    if (gg < EPS) break;

    if (gradsearch == -1) beta = 0.0;
    if (beta == 0.0) gradsearch = 1;
    else gradsearch = 0;

    for (i = 0; i < ndof; i++) {
      g[i] = f[i];
      h[i] = g[i] + beta*h[i];
    }
    for (i = 0; i < ndof_extra; i++) {
      gextra[i] = fextra[i];
      hextra[i] = gextra[i] + beta*hextra[i];
    }

    // output for thermo, dump, restart files

    if (output->next == update->ntimestep) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}
