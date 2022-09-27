// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "min_cg_kokkos.h"
#include <mpi.h>
#include <cmath>
#include "update.h"
#include "output.h"
#include "timer.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_minimize_kokkos.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

/* ---------------------------------------------------------------------- */

MinCGKokkos::MinCGKokkos(LAMMPS *lmp) : MinLineSearchKokkos(lmp)
{
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

int MinCGKokkos::iterate(int maxiter)
{
  int fail,ntimestep;
  double beta,gg,dot[2],dotall[2],fdotf;

  fix_minimize_kk->k_vectors.sync<LMPDeviceType>();
  fix_minimize_kk->k_vectors.modify<LMPDeviceType>();

  atomKK->sync(Device,F_MASK);

  // nlimit = max # of CG iterations before restarting
  // set to ndoftotal unless too big

  int nlimit = static_cast<int> (MIN(MAXSMALLINT,ndoftotal));

  // initialize working vectors

  {
    // local variables for lambda capture

    auto l_h = h;
    auto l_g = g;
    auto l_fvec = fvec;

    Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
      l_h[i] = l_fvec[i];
      l_g[i] = l_fvec[i];
    });
  }
  if (nextra_global)
    for (int i = 0; i < nextra_global; i++) hextra[i] = gextra[i] = fextra[i];

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

    s_double2 sdot;
    {
      // local variables for lambda capture

      auto l_g = g;
      auto l_fvec = fvec;

      Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, s_double2& sdot) {
        sdot.d0 += l_fvec[i]*l_fvec[i];
        sdot.d1 += l_fvec[i]*l_g[i];
      },sdot);
    }
    dot[0] = sdot.d0;
    dot[1] = sdot.d1;
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);
    if (nextra_global)
      for (int i = 0; i < nextra_global; i++) {
        dotall[0] += fextra[i]*fextra[i];
        dotall[1] += fextra[i]*gextra[i];
      }

    fdotf = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotf = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf = fnorm_inf();   // infinite force norm
      else if (normstyle == TWO) fdotf = dotall[0];     // same as fnorm_sqr(), Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      if (fdotf < update->ftol*update->ftol) return FTOL;
    }

    // update new search direction h from new f = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be Fletcher-Reeves
    // reinitialize CG every ndof iterations by setting beta = 0.0

    beta = MAX(0.0,(dotall[0] - dotall[1])/gg);
    if ((niter+1) % nlimit == 0) beta = 0.0;
    gg = dotall[0];

    {
      // local variables for lambda capture

      auto l_h = h;
      auto l_g = g;
      auto l_fvec = fvec;

      Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
        l_g[i] = l_fvec[i];
        l_h[i] = l_g[i] + beta*l_h[i];
      });
    }
    if (nextra_global)
      for (int i = 0; i < nextra_global; i++) {
        gextra[i] = fextra[i];
        hextra[i] = gextra[i] + beta*hextra[i];
      }

    // reinitialize CG if new search direction h is not downhill

    double dot_0 = 0.0;

    {
      // local variables for lambda capture

      auto l_h = h;
      auto l_g = g;

      Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(const int& i, double& dot_0) {
        dot_0 += l_g[i]*l_h[i];
      },dot_0);
    }
    dot[0] = dot_0;
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);
    if (nextra_global)
      for (int i = 0; i < nextra_global; i++)
        dotall[0] += gextra[i]*hextra[i];

    if (dotall[0] <= 0.0) {
      // local variables for lambda capture

      auto l_h = h;
      auto l_g = g;

      Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
         l_h[i] = l_g[i];
      });
      if (nextra_global)
        for (int i = 0; i < nextra_global; i++) hextra[i] = gextra[i];
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      atomKK->sync(Host,ALL_MASK);

      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
