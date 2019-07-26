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

#include <mpi.h>
#include <cmath>
#include <cstring>
#include "min_cg.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

/* ---------------------------------------------------------------------- */

MinCG::MinCG(LAMMPS *lmp) : MinLineSearch(lmp) {}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

int MinCG::iterate(int maxiter)
{
  int i,m,n,fail,ntimestep;
  double beta,gg,dot[2],dotall[2],fmax,fmaxall;
  double *fatom,*gatom,*hatom;

  // nlimit = max # of CG iterations before restarting
  // set to ndoftotal unless too big

  int nlimit = static_cast<int> (MIN(MAXSMALLINT,ndoftotal));

  // initialize working vectors

  for (i = 0; i < nvec; i++) h[i] = g[i] = fvec[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      gatom = gextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hatom[i] = gatom[i] = fatom[i];
    }
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) hextra[i] = gextra[i] = fextra[i];

  gg = fnorm_sqr();

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    fail = (this->*linemin)(ecurrent,alpha_final);
    if (fail) return fail;

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) <
        update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    fmax = fmaxall = 0.0;
    dot[0] = dot[1] = 0.0;
    for (i = 0; i < nvec; i++) {
      dot[0] += fvec[i]*fvec[i];
      dot[1] += fvec[i]*g[i];
      fmax = MAX(fmax,fvec[i]*fvec[i]);
    }
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
        gatom = gextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) {
          dot[0] += fatom[i]*fatom[i];
          dot[1] += fatom[i]*gatom[i];
	  fmax = MAX(fmax,fatom[i]*fatom[i]);
        }
      }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&fmax,&fmaxall,2,MPI_DOUBLE,MPI_MAX,world);
    if (nextra_global)
      for (i = 0; i < nextra_global; i++) {
        dotall[0] += fextra[i]*fextra[i];
        dotall[1] += fextra[i]*gextra[i];
      }

    if (normstyle == 1) {	// max force norm
      if (fmax < update->ftol*update->ftol) return FTOL;
    } else {			// Euclidean force norm
      if (dotall[0] < update->ftol*update->ftol) return FTOL;
    }

    // update new search direction h from new f = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be Fletcher-Reeves
    // reinitialize CG every ndof iterations by setting beta = 0.0

    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    if ((niter+1) % nlimit == 0) beta = 0.0;
    gg = dotall[0];

    for (i = 0; i < nvec; i++) {
      g[i] = fvec[i];
      h[i] = g[i] + beta*h[i];
    }
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
        gatom = gextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) {
          gatom[i] = fatom[i];
          hatom[i] = gatom[i] + beta*hatom[i];
        }
      }
    if (nextra_global)
      for (i = 0; i < nextra_global; i++) {
        gextra[i] = fextra[i];
        hextra[i] = gextra[i] + beta*hextra[i];
      }

    // reinitialize CG if new search direction h is not downhill

    dot[0] = 0.0;
    for (i = 0; i < nvec; i++) dot[0] += g[i]*h[i];
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        gatom = gextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) dot[0] += gatom[i]*hatom[i];
      }
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);
    if (nextra_global)
      for (i = 0; i < nextra_global; i++)
        dotall[0] += gextra[i]*hextra[i];

    if (dotall[0] <= 0.0) {
      for (i = 0; i < nvec; i++) h[i] = g[i];
      if (nextra_atom)
        for (m = 0; m < nextra_atom; m++) {
          gatom = gextra_atom[m];
          hatom = hextra_atom[m];
          n = extra_nlen[m];
          for (i = 0; i < n; i++) hatom[i] = gatom[i];
        }
      if (nextra_global)
        for (i = 0; i < nextra_global; i++) hextra[i] = gextra[i];
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
