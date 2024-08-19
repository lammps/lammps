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

#include "angle_hybrid_kokkos.h"

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

AngleHybridKokkos::AngleHybridKokkos(LAMMPS *lmp) : AngleHybrid(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;

  execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

AngleHybridKokkos::~AngleHybridKokkos()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void AngleHybridKokkos::compute(int eflag, int vflag)
{
  // save ptrs to original anglelist

  int nanglelist_orig = neighbor->nanglelist;
  neighborKK->k_anglelist.sync_device();
  auto k_anglelist_orig = neighborKK->k_anglelist;
  auto d_anglelist_orig = k_anglelist_orig.d_view;
  auto d_nanglelist = k_nanglelist.d_view;
  auto h_nanglelist = k_nanglelist.h_view;

  // if this is re-neighbor step, create sub-style anglelists
  // nanglelist[] = length of each sub-style list
  // realloc sub-style anglelist if necessary
  // load sub-style anglelist with 3 values from original anglelist

  if (neighbor->ago == 0) {
    Kokkos::deep_copy(d_nanglelist,0);

    k_map.sync_device();
    auto d_map = k_map.d_view;

    Kokkos::parallel_for(nanglelist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_anglelist_orig(i,3)];
      if (m >= 0) Kokkos::atomic_increment(&d_nanglelist[m]);
    });

    k_nanglelist.modify_device();
    k_nanglelist.sync_host();

    maxangle_all = 0;
    for (int m = 0; m < nstyles; m++)
      if (h_nanglelist[m] > maxangle_all)
        maxangle_all = h_nanglelist[m] + EXTRA;

    if (k_anglelist.d_view.extent(1) < maxangle_all)
      MemKK::realloc_kokkos(k_anglelist, "angle_hybrid:anglelist", nstyles, maxangle_all, 4);
    auto d_anglelist = k_anglelist.d_view;

    Kokkos::deep_copy(d_nanglelist,0);

    Kokkos::parallel_for(nanglelist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_anglelist_orig(i,3)];
      if (m < 0) return;
      const int n = Kokkos::atomic_fetch_add(&d_nanglelist[m],1);
      d_anglelist(m,n,0) = d_anglelist_orig(i,0);
      d_anglelist(m,n,1) = d_anglelist_orig(i,1);
      d_anglelist(m,n,2) = d_anglelist_orig(i,2);
      d_anglelist(m,n,3) = d_anglelist_orig(i,3);
    });
  }

  // call each sub-style's compute function
  // set neighbor->anglelist to sub-style anglelist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag, vflag);

  k_nanglelist.modify_device();
  k_nanglelist.sync_host();

  for (int m = 0; m < nstyles; m++) {
    neighbor->nanglelist = h_nanglelist[m];
    auto k_anglelist_m = Kokkos::subview(k_anglelist,m,Kokkos::ALL,Kokkos::ALL);
    k_anglelist_m.modify_device();
    neighborKK->k_anglelist = k_anglelist_m;

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
    if (cvflag_atom) {
      int n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double **cvatom_substyle = styles[m]->cvatom;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < 9; j++) cvatom[i][j] += cvatom_substyle[i][j];
    }
  }

  // restore ptrs to original anglelist

  neighbor->nanglelist = nanglelist_orig;
  neighborKK->k_anglelist = k_anglelist_orig;
}

/* ---------------------------------------------------------------------- */

void AngleHybridKokkos::allocate()
{
  allocated = 1;
  int np1 = atom->nangletypes + 1;

  memoryKK->create_kokkos(k_map, map, np1, "angle:map");
  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;

  k_nanglelist = DAT::tdual_int_1d("angle:nanglelist", nstyles);
}

/* ---------------------------------------------------------------------- */

void AngleHybridKokkos::deallocate()
{
  if (!allocated) return;

  allocated = 0;

  memory->destroy(setflag);
  memoryKK->destroy_kokkos(k_map,map);
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void AngleHybridKokkos::coeff(int narg, char **arg)
{
  AngleHybrid::coeff(narg,arg);

  k_map.modify_host();
}

/* ---------------------------------------------------------------------- */

void AngleHybridKokkos::init_style()
{
  AngleHybrid::init_style();

  for (int m = 0; m < nstyles; m++) {
    if (!styles[m]->kokkosable)
      error->all(FLERR,"Must use only Kokkos-enabled angle styles with angle_style hybrid/kk");

    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double AngleHybridKokkos::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  bytes += (double) maxcvatom * 9 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += (double) maxangle_all * 4 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
