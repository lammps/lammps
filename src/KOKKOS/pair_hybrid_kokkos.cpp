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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include "pair_hybrid_kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "respa.h"
#include "atom_masks.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridKokkos::PairHybridKokkos(LAMMPS *lmp) : PairHybrid(lmp)
{
  atomKK = (AtomKokkos *) atom;

 // prevent overlapping host/device computation, which isn't
 //  yet supported by pair_hybrid_kokkos
 execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

PairHybridKokkos::~PairHybridKokkos()
{

}

/* ----------------------------------------------------------------------
  call each sub-style's compute() or compute_outer() function
  accumulate sub-style global/peratom energy/virial in hybrid
  for global vflag = 1:
    each sub-style computes own virial[6]
    sum sub-style virial[6] to hybrid's virial[6]
  for global vflag = 2:
    call sub-style with adjusted vflag to prevent it calling
      virial_fdotr_compute()
    hybrid calls virial_fdotr_compute() on final accumulated f
------------------------------------------------------------------------- */

void PairHybridKokkos::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = 2, then
  // reset vflag as if global component were 1
  // necessary since one or more sub-styles cannot compute virial as F dot r

  int neighflag = lmp->kokkos->neighflag;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  if (no_virial_fdotr_compute && vflag % 4 == 2) vflag = 1 + vflag/4 * 4;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // check if global component of incoming vflag = 2
  // if so, reset vflag passed to substyle as if it were 0
  // necessary so substyle will not invoke virial_fdotr_compute()

  int vflag_substyle;
  if (vflag % 4 == 2) vflag_substyle = vflag/4 * 4;
  else vflag_substyle = vflag;

  double *saved_special = save_special();

  // check if we are running with r-RESPA using the hybrid keyword

  Respa *respa = NULL;
  respaflag = 0;
  if (strstr(update->integrate_style,"respa")) {
    respa = (Respa *) update->integrate;
    if (respa->nhybrid_styles > 0) respaflag = 1;
  }

  for (m = 0; m < nstyles; m++) {

    set_special(m);

    if (!respaflag || (respaflag && respa->hybrid_compute[m])) {

      // invoke compute() unless compute flag is turned off or
      // outerflag is set and sub-style has a compute_outer() method

      if (styles[m]->compute_flag == 0) continue;
      atomKK->sync(styles[m]->execution_space,styles[m]->datamask_read);
      if (outerflag && styles[m]->respa_enable)
        styles[m]->compute_outer(eflag,vflag_substyle);
      else styles[m]->compute(eflag,vflag_substyle);
      atomKK->modified(styles[m]->execution_space,styles[m]->datamask_modify);
    }

    restore_special(saved_special);

    // jump to next sub-style if r-RESPA does not want global accumulated data

    if (respaflag && !respa->tally_global) continue;

    if (eflag_global) {
      eng_vdwl += styles[m]->eng_vdwl;
      eng_coul += styles[m]->eng_coul;
    }
    if (vflag_global) {
      for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
    }
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++)
          vatom[i][j] += vatom_substyle[i][j];
    }
  }

  delete [] saved_special;

  // perform virial_fdotr on device

  atomKK->sync(Device,X_MASK|F_MASK);
  x = atomKK->k_x.view<LMPDeviceType>();
  f = atomKK->k_f.view<LMPDeviceType>();

  if (vflag_fdotr)
    pair_virial_fdotr_compute(this);
}
