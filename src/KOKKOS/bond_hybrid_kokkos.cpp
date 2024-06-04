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

#include "bond_hybrid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cstring>

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

BondHybridKokkos::BondHybridKokkos(LAMMPS *lmp) : BondHybrid(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;

  execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

BondHybridKokkos::~BondHybridKokkos()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void BondHybridKokkos::compute(int eflag, int vflag)
{
  // save ptrs to original bondlist

  int nbondlist_orig = neighbor->nbondlist;
  neighborKK->k_bondlist.sync_device();
  auto k_bondlist_orig = neighborKK->k_bondlist;
  auto d_bondlist_orig = k_bondlist_orig.d_view;
  auto d_nbondlist = k_nbondlist.d_view;
  auto h_nbondlist = k_nbondlist.h_view;

  // if this is re-neighbor step, create sub-style bondlists
  // nbondlist[] = length of each sub-style list
  // realloc sub-style bondlist if necessary
  // load sub-style bondlist with 3 values from original bondlist

  if (neighbor->ago == 0) {
    Kokkos::deep_copy(d_nbondlist,0);

    k_map.sync_device();
    auto d_map = k_map.d_view;

    Kokkos::parallel_for(nbondlist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_bondlist_orig(i,2)];
      if (m >= 0) Kokkos::atomic_increment(&d_nbondlist[m]);
    });

    k_nbondlist.modify_device();
    k_nbondlist.sync_host();

    maxbond_all = 0;
    for (int m = 0; m < nstyles; m++)
      if (h_nbondlist[m] > maxbond_all)
        maxbond_all = h_nbondlist[m] + EXTRA;

    if (k_bondlist.d_view.extent(1) < maxbond_all)
      MemKK::realloc_kokkos(k_bondlist, "bond_hybrid:bondlist", nstyles, maxbond_all, 3);
    auto d_bondlist = k_bondlist.d_view;

    Kokkos::deep_copy(d_nbondlist,0);

    Kokkos::parallel_for(nbondlist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_bondlist_orig(i,2)];
      if (m < 0) return;
      const int n = Kokkos::atomic_fetch_add(&d_nbondlist[m],1);
      d_bondlist(m,n,0) = d_bondlist_orig(i,0);
      d_bondlist(m,n,1) = d_bondlist_orig(i,1);
      d_bondlist(m,n,2) = d_bondlist_orig(i,2);
    });
  }

  // call each sub-style's compute function
  // set neighbor->bondlist to sub-style bondlist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag, vflag);

  k_nbondlist.modify_device();
  k_nbondlist.sync_host();

  for (int m = 0; m < nstyles; m++) {
    neighbor->nbondlist = h_nbondlist[m];
    auto k_bondlist_m = Kokkos::subview(k_bondlist,m,Kokkos::ALL,Kokkos::ALL);
    k_bondlist_m.modify_device();
    neighborKK->k_bondlist = k_bondlist_m;

    auto style = styles[m];
    atomKK->sync(style->execution_space,style->datamask_read);
    style->compute(eflag, vflag);
    atomKK->modified(style->execution_space,style->datamask_modify);

    if (eflag_global) energy += style->energy;
    if (vflag_global)
      for (int n = 0; n < 6; n++) virial[n] += style->virial[n];

    if (eflag_atom) {
      int n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (int i = 0; i < n; i++) eatom[i] += eatom_substyle[i];
    }
    if (vflag_atom) {
      int n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < 6; j++) vatom[i][j] += vatom_substyle[i][j];
    }
  }

  // restore ptrs to original bondlist

  neighbor->nbondlist = nbondlist_orig;
  neighborKK->k_bondlist = k_bondlist_orig;
}

/* ---------------------------------------------------------------------- */

void BondHybridKokkos::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memoryKK->create_kokkos(k_map, map, n + 1, "bond:map");
  memory->create(setflag, n + 1, "bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  k_nbondlist = DAT::tdual_int_1d("bond:nbondlist", nstyles);
}

/* ---------------------------------------------------------------------- */

void BondHybridKokkos::deallocate()
{
  if (!allocated) return;

  allocated = 0;

  memory->destroy(setflag);
  memoryKK->destroy_kokkos(k_map,map);
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void BondHybridKokkos::coeff(int narg, char **arg)
{
  BondHybrid::coeff(narg,arg);

  k_map.modify_host();
}

/* ---------------------------------------------------------------------- */

void BondHybridKokkos::init_style()
{
  BondHybrid::init_style();

  for (int m = 0; m < nstyles; m++) {
    if (!styles[m]->kokkosable)
      error->all(FLERR,"Must use only Kokkos-enabled bond styles with bond_style hybrid/kk");

    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double BondHybridKokkos::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += (double) maxbond_all * 3 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
