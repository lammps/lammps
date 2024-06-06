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

#include "pair_hybrid_kokkos.h"

#include "atom_kokkos.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "atom_masks.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridKokkos::PairHybridKokkos(LAMMPS *lmp) : PairHybrid(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;

  execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHybridKokkos::init_style()
{
  PairHybrid::init_style();

  for (int m = 0; m < nstyles; m++)
    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
}

/* ----------------------------------------------------------------------
  call each sub-style's compute() or compute_outer() function
  accumulate sub-style global/peratom energy/virial in hybrid
  for global vflag = VIRIAL_PAIR:
    each sub-style computes own virial[6]
    sum sub-style virial[6] to hybrid's virial[6]
  for global vflag = VIRIAL_FDOTR:
    call sub-style with adjusted vflag to prevent it calling
      virial_fdotr_compute()
    hybrid calls virial_fdotr_compute() on final accumulated f
------------------------------------------------------------------------- */

void PairHybridKokkos::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // check if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag as if global component were VIRIAL_PAIR
  // necessary since one or more sub-styles cannot compute virial as F dot r

  if (lmp->kokkos->neighflag == FULL) no_virial_fdotr_compute = 1;

  if (no_virial_fdotr_compute && (vflag & VIRIAL_FDOTR))
    vflag = VIRIAL_PAIR | (vflag & ~VIRIAL_FDOTR);

  ev_init(eflag,vflag);

  // check if global component of incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag passed to substyle so VIRIAL_FDOTR is turned off
  // necessary so substyle will not invoke virial_fdotr_compute()

  int vflag_substyle;
  if (vflag & VIRIAL_FDOTR) vflag_substyle = vflag & ~VIRIAL_FDOTR;
  else vflag_substyle = vflag;

  double *saved_special = save_special();

  // check if we are running with r-RESPA using the hybrid keyword

  Respa *respa = nullptr;
  respaflag = 0;
  if (utils::strmatch(update->integrate_style,"^respa")) {
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

    // substyles may be CENTROID_SAME or CENTROID_AVAIL

    if (cvflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      if (styles[m]->centroidstressflag == CENTROID_AVAIL) {
        double **cvatom_substyle = styles[m]->cvatom;
        for (i = 0; i < n; i++)
          for (j = 0; j < 9; j++)
            cvatom[i][j] += cvatom_substyle[i][j];
      } else {
        double **vatom_substyle = styles[m]->vatom;
        for (i = 0; i < n; i++) {
          for (j = 0; j < 6; j++) {
            cvatom[i][j] += vatom_substyle[i][j];
          }
          for (j = 6; j < 9; j++) {
            cvatom[i][j] += vatom_substyle[i][j-3];
          }
        }
      }
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
