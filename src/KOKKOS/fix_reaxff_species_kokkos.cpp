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

/* ----------------------------------------------------------------------
   Contributing authors: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "fix_reaxff_species_kokkos.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include "fix_ave_atom.h"
#include "pair_reaxff_kokkos.h"
#include "reaxff_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxFFSpeciesKokkos::FixReaxFFSpeciesKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReaxFFSpecies(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;

  // NOTE: Could improve performance if a Kokkos version of ComputeSpecAtom is added

  datamask_read = X_MASK | V_MASK | Q_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::init()
{
  Pair* pair_kk = force->pair_match("^reax../kk",0);
  if (pair_kk == nullptr) error->all(FLERR,"Cannot use fix reaxff/species/kk without "
                  "pair_style reaxff/kk");

  FixReaxFFSpecies::init();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::FindMolecule()
{
  int i,j,ii,jj,inum,itype,jtype,loop,looptot;
  int change,done,anychange;
  int *mask = atom->mask;
  double bo_tmp,bo_cut;
  double **spec_atom = f_SPECBOND->array_atom;

  inum = reaxff->list->inum;
  typename ArrayTypes<LMPHostType>::t_int_1d ilist;
  if (reaxff->execution_space == Host) {
    NeighListKokkos<LMPHostType>* k_list = static_cast<NeighListKokkos<LMPHostType>*>(reaxff->list);
    k_list->k_ilist.sync<LMPHostType>();
    ilist = k_list->k_ilist.h_view;
  } else {
    NeighListKokkos<LMPDeviceType>* k_list = static_cast<NeighListKokkos<LMPDeviceType>*>(reaxff->list);
    k_list->k_ilist.sync<LMPHostType>();
    ilist = k_list->k_ilist.h_view;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      clusterID[i] = atom->tag[i];
      x0[i].x = spec_atom[i][1];
      x0[i].y = spec_atom[i][2];
      x0[i].z = spec_atom[i][3];
    }
    else clusterID[i] = 0.0;
  }

  loop = 0;
  while (true) {
    comm->forward_comm(this);
    loop ++;

    change = 0;
    while (true) {
      done = 1;

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        itype = atom->type[i];

        for (jj = 0; jj < MAXSPECBOND; jj++) {
          j = reaxff->tmpid[i][jj];

          if ((j == 0) && (j < i)) continue;
          if (!(mask[j] & groupbit)) continue;

          if (clusterID[i] == clusterID[j]
            && x0[i].x == x0[j].x && x0[i].y == x0[j].y && x0[i].z == x0[j].z) continue;

          jtype = atom->type[j];
          bo_cut = BOCut[itype][jtype];
          bo_tmp = spec_atom[i][jj+7];

          if (bo_tmp > bo_cut) {
            clusterID[i] = clusterID[j] = MIN(clusterID[i], clusterID[j]);
            x0[i] = x0[j] = chAnchor(x0[i], x0[j]);
            done = 0;
          }
        }
      }
      if (!done) change = 1;
      if (done) break;
    }
    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;

    MPI_Allreduce(&loop,&looptot,1,MPI_INT,MPI_SUM,world);
    if (looptot >= 400*nprocs) break;

  }
}
