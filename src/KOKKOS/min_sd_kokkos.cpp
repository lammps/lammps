// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "min_sd_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_minimize_kokkos.h"
#include "output.h"
#include "timer.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

static constexpr double EPS_ENERGY = 1.0e-8;

/* ---------------------------------------------------------------------- */

MinSDKokkos::MinSDKokkos(LAMMPS *lmp) : MinLineSearchKokkos(lmp)
{
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;
}

/* ----------------------------------------------------------------------
   set the search direction h to the downhill gradient f = -Grad(x)
------------------------------------------------------------------------- */

void MinSDKokkos::set_search_direction()
{
  // local variables for lambda capture

  auto l_h = h;

  if constexpr (F_LAYOUTRIGHT) {
    auto l_fvec = fvec;
    Kokkos::parallel_for(nvec, LAMMPS_LAMBDA(const int& i) {
      l_h[i] = static_cast<KK_FLOAT>(l_fvec[i]);
    });
  } else {
    auto l_f = atomKK->k_f.view_device();
    Kokkos::parallel_for(atom->nlocal, LAMMPS_LAMBDA(const int& i) {
      const int j = i*3;
      l_h[j] = static_cast<KK_FLOAT>(l_f(i,0));
      l_h[j+1] = static_cast<KK_FLOAT>(l_f(i,1));
      l_h[j+2] = static_cast<KK_FLOAT>(l_f(i,2));
    });
  }

  if (nextra_global)
    for (int i = 0; i < nextra_global; i++) hextra[i] = fextra[i];
}

/* ----------------------------------------------------------------------
   minimization via steepest descent
------------------------------------------------------------------------- */

int MinSDKokkos::iterate(int maxiter)
{
  int fail,ntimestep;
  double fdotf;

  fix_minimize_kk->k_vectors.sync_device();
  fix_minimize_kk->k_vectors.modify_device();

  atomKK->sync(Device,F_MASK);

  // initialize working vectors

  set_search_direction();

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along h from current position x
    // h = downhill gradient direction

    eprevious = ecurrent;
    fail = (this->*linemin)(ecurrent,alpha_final);
    if (fail) return fail;

    // the line search ends with a force evaluation, and with
    // modify->min_reset_ref() when there are extra global dof.  the styles and
    // fixes it runs may be non-KOKKOS ones, which claim the host side of f and
    // leave the device view fvec read by set_search_direction() below stale

    atomKK->sync(Device,F_MASK);

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) <
        update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    fdotf = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotf = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf = fnorm_inf();   // infinite force norm
      else if (normstyle == TWO) fdotf = fnorm_sqr();   // Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      if (fdotf < update->ftol*update->ftol) return FTOL;
    }

    // set new search direction h to f = -Grad(x)

    set_search_direction();

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
