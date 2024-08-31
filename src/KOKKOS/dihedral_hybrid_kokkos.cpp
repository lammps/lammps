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

#include "dihedral_hybrid_kokkos.h"

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

DihedralHybridKokkos::DihedralHybridKokkos(LAMMPS *lmp) : DihedralHybrid(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;

  execution_space = Device;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

DihedralHybridKokkos::~DihedralHybridKokkos()
{
  deallocate();
}

/* ---------------------------------------------------------------------- */

void DihedralHybridKokkos::compute(int eflag, int vflag)
{
  // save ptrs to original dihedrallist

  int ndihedrallist_orig = neighbor->ndihedrallist;
  neighborKK->k_dihedrallist.sync_device();
  auto k_dihedrallist_orig = neighborKK->k_dihedrallist;
  auto d_dihedrallist_orig = k_dihedrallist_orig.d_view;
  auto d_ndihedrallist = k_ndihedrallist.d_view;
  auto h_ndihedrallist = k_ndihedrallist.h_view;

  // if this is re-neighbor step, create sub-style dihedrallists
  // ndihedrallist[] = length of each sub-style list
  // realloc sub-style dihedrallist if necessary
  // load sub-style dihedrallist with 3 values from original dihedrallist

  if (neighbor->ago == 0) {
    Kokkos::deep_copy(d_ndihedrallist,0);

    k_map.sync_device();
    auto d_map = k_map.d_view;

    Kokkos::parallel_for(ndihedrallist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_dihedrallist_orig(i,4)];
      if (m >= 0) Kokkos::atomic_increment(&d_ndihedrallist[m]);
    });

    k_ndihedrallist.modify_device();
    k_ndihedrallist.sync_host();

    maxdihedral_all = 0;
    for (int m = 0; m < nstyles; m++)
      if (h_ndihedrallist[m] > maxdihedral_all)
        maxdihedral_all = h_ndihedrallist[m] + EXTRA;

    if (k_dihedrallist.d_view.extent(1) < maxdihedral_all)
      MemKK::realloc_kokkos(k_dihedrallist, "dihedral_hybrid:dihedrallist", nstyles, maxdihedral_all, 5);
    auto d_dihedrallist = k_dihedrallist.d_view;

    Kokkos::deep_copy(d_ndihedrallist,0);

    Kokkos::parallel_for(ndihedrallist_orig,LAMMPS_LAMBDA(int i) {
      const int m = d_map[d_dihedrallist_orig(i,4)];
      if (m < 0) return;
      const int n = Kokkos::atomic_fetch_add(&d_ndihedrallist[m],1);
      d_dihedrallist(m,n,0) = d_dihedrallist_orig(i,0);
      d_dihedrallist(m,n,1) = d_dihedrallist_orig(i,1);
      d_dihedrallist(m,n,2) = d_dihedrallist_orig(i,2);
      d_dihedrallist(m,n,3) = d_dihedrallist_orig(i,3);
      d_dihedrallist(m,n,4) = d_dihedrallist_orig(i,4);
    });
  }

  // call each sub-style's compute function
  // set neighbor->dihedrallist to sub-style dihedrallist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag, vflag);

  k_ndihedrallist.modify_device();
  k_ndihedrallist.sync_host();

  for (int m = 0; m < nstyles; m++) {
    neighbor->ndihedrallist = h_ndihedrallist[m];
    auto k_dihedrallist_m = Kokkos::subview(k_dihedrallist,m,Kokkos::ALL,Kokkos::ALL);
    k_dihedrallist_m.modify_device();
    neighborKK->k_dihedrallist = k_dihedrallist_m;

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

  // restore ptrs to original dihedrallist

  neighbor->ndihedrallist = ndihedrallist_orig;
  neighborKK->k_dihedrallist = k_dihedrallist_orig;
}

/* ---------------------------------------------------------------------- */

void DihedralHybridKokkos::allocate()
{
  allocated = 1;
  int np1 = atom->ndihedraltypes + 1;

  memoryKK->create_kokkos(k_map, map, np1, "dihedral:map");
  memory->create(setflag, np1, "dihedral:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;

  k_ndihedrallist = DAT::tdual_int_1d("dihedral:ndihedrallist", nstyles);
}

/* ---------------------------------------------------------------------- */

void DihedralHybridKokkos::deallocate()
{
  if (!allocated) return;

  allocated = 0;

  memory->destroy(setflag);
  memoryKK->destroy_kokkos(k_map,map);
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void DihedralHybridKokkos::coeff(int narg, char **arg)
{
  DihedralHybrid::coeff(narg,arg);

  k_map.modify_host();
}

/* ---------------------------------------------------------------------- */

void DihedralHybridKokkos::init_style()
{
  DihedralHybrid::init_style();

  for (int m = 0; m < nstyles; m++) {
    if (!styles[m]->kokkosable)
      error->all(FLERR,"Must use only Kokkos-enabled dihedral styles with dihedral_style hybrid/kk");

    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double DihedralHybridKokkos::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  bytes += (double) maxcvatom * 9 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += (double) maxdihedral_all * 5 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
