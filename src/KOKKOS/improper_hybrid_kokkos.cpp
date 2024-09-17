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

#include "improper_hybrid_kokkos.h"

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

ImproperHybridKokkos::ImproperHybridKokkos(LAMMPS *lmp) : ImproperHybrid(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;

  execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

ImproperHybridKokkos::~ImproperHybridKokkos()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void ImproperHybridKokkos::compute(int eflag, int vflag)
{

  // save ptrs to original improperlist

  int nimproperlist_orig = neighbor->nimproperlist;
  neighborKK->k_improperlist.sync_device();
  auto k_improperlist_orig = neighborKK->k_improperlist;
  auto d_improperlist_orig = k_improperlist_orig.d_view;
  auto d_nimproperlist = k_nimproperlist.d_view;
  auto h_nimproperlist = k_nimproperlist.h_view;

  // if this is re-neighbor step, create sub-style improperlists
  // nimproperlist[] = length of each sub-style list
  // realloc sub-style improperlist if necessary
  // load sub-style improperlist with 3 values from original improperlist

  if (neighbor->ago == 0) {
    Kokkos::deep_copy(d_nimproperlist,0);

    k_map.sync_device();
    auto d_map = k_map.d_view;

    Kokkos::parallel_for(nimproperlist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_improperlist_orig(i,4)];
      if (m >= 0) Kokkos::atomic_increment(&d_nimproperlist[m]);
    });

    k_nimproperlist.modify_device();
    k_nimproperlist.sync_host();

    maximproper_all = 0;
    for (int m = 0; m < nstyles; m++)
      if (h_nimproperlist[m] > maximproper_all)
        maximproper_all = h_nimproperlist[m] + EXTRA;

    if (k_improperlist.d_view.extent(1) < maximproper_all)
      MemKK::realloc_kokkos(k_improperlist, "improper_hybrid:improperlist", nstyles, maximproper_all, 5);
    auto d_improperlist = k_improperlist.d_view;

    Kokkos::deep_copy(d_nimproperlist,0);

    Kokkos::parallel_for(nimproperlist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_improperlist_orig(i,4)];
      if (m < 0) return;
      const int n = Kokkos::atomic_fetch_add(&d_nimproperlist[m],1);
      d_improperlist(m,n,0) = d_improperlist_orig(i,0);
      d_improperlist(m,n,1) = d_improperlist_orig(i,1);
      d_improperlist(m,n,2) = d_improperlist_orig(i,2);
      d_improperlist(m,n,3) = d_improperlist_orig(i,3);
      d_improperlist(m,n,4) = d_improperlist_orig(i,4);
    });
  }

  // call each sub-style's compute function
  // set neighbor->improperlist to sub-style improperlist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag, vflag);

  k_nimproperlist.modify_device();
  k_nimproperlist.sync_host();

  for (int m = 0; m < nstyles; m++) {
    neighbor->nimproperlist = h_nimproperlist[m];
    auto k_improperlist_m = Kokkos::subview(k_improperlist,m,Kokkos::ALL,Kokkos::ALL);
    k_improperlist_m.modify_device();
    neighborKK->k_improperlist = k_improperlist_m;

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

  // restore ptrs to original improperlist

  neighbor->nimproperlist = nimproperlist_orig;
  neighborKK->k_improperlist = k_improperlist_orig;
}

/* ---------------------------------------------------------------------- */

void ImproperHybridKokkos::allocate()
{
  allocated = 1;
  int np1 = atom->nimpropertypes + 1;

  memoryKK->create_kokkos(k_map, map, np1, "improper:map");
  memory->create(setflag, np1, "improper:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;

  k_nimproperlist = DAT::tdual_int_1d("improper:nimproperlist", nstyles);
}

/* ---------------------------------------------------------------------- */

void ImproperHybridKokkos::deallocate()
{
  if (!allocated) return;

  allocated = 0;

  memory->destroy(setflag);
  memoryKK->destroy_kokkos(k_map,map);
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void ImproperHybridKokkos::coeff(int narg, char **arg)
{
  ImproperHybrid::coeff(narg,arg);

  k_map.modify_host();
}

/* ---------------------------------------------------------------------- */

void ImproperHybridKokkos::init_style()
{
  ImproperHybrid::init_style();

  for (int m = 0; m < nstyles; m++) {
    if (!styles[m]->kokkosable)
      error->all(FLERR,"Must use only Kokkos-enabled improper styles with improper_style hybrid/kk");

    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ImproperHybridKokkos::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  bytes += (double) maxcvatom * 9 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += (double) maximproper_all * 5 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
