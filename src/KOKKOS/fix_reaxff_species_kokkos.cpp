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

/* ----------------------------------------------------------------------
   Contributing authors: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "fix_reaxff_species_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "pair_reaxff_kokkos.h"
#include "update.h"

#include "fix_ave_atom.h"
#include "pair_reaxff.h"
#include "reaxff_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxFFSpeciesKokkos::FixReaxFFSpeciesKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReaxFFSpecies(lmp, narg, arg)
{
  // not all functions in FixReaxFFSpecies are ported to KOKKOS
  //  so set kokkosable flag to zero

  kokkosable = 0;

  atomKK = (AtomKokkos *) atom;

  // the per-atom data and the bonds are only copied to the host on the steps
  // when they are needed, see post_integrate()

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::init()
{
  Pair* pair_kk = force->pair_match("^reax../kk",0);
  if (pair_kk == nullptr)
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot use fix reaxff/species/kk without pair_style reaxff/kk");

  FixReaxFFSpecies::init();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::setup(int vflag)
{
  sync_host_data();
  FixReaxFFSpecies::setup(vflag);
  atomKK->modified(Host, DVECTOR_MASK);
}

/* ----------------------------------------------------------------------
   the averaging of the bond data samples on the steps of the averaging
   windows that end on multiples of Nfreq, which also include the steps of
   the analysis.  on all other steps, no per-atom data is used on the host.
------------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::post_integrate()
{
  if (sample_step()) {
    sync_host_data();
    FixReaxFFSpecies::post_integrate();
    atomKK->modified(Host, DVECTOR_MASK);
  } else {
    FixReaxFFSpecies::post_integrate();
  }
}

/* ----------------------------------------------------------------------
   true if the current step is within an averaging window on a step that
   is a multiple of Nevery.  this includes all steps on which the averaging
   samples the bond data, see FixAveAtom::nextvalid()
------------------------------------------------------------------------- */

bool FixReaxFFSpeciesKokkos::sample_step() const
{
  const bigint ntimestep = update->ntimestep;
  if (ntimestep % nevery) return false;
  const bigint ahead = (nfreq - (ntimestep % nfreq)) % nfreq;
  return ahead <= static_cast<bigint>(nrepeat - 1) * nevery;
}

/* ----------------------------------------------------------------------
   make the per-atom data and the bonds of pair reaxff/kk current on the host
------------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::sync_host_data()
{
  atomKK->sync(Host, X_MASK | V_MASK | Q_MASK | MASK_MASK | DVECTOR_MASK | TYPE_MASK |
                         RMASS_MASK | TAG_MASK);
  if (reaxff->execution_space == HostKK)
    static_cast<PairReaxFFKokkos<LMPHostType> *>(reaxff)->sync_species_host();
  else
    static_cast<PairReaxFFKokkos<LMPDeviceType> *>(reaxff)->sync_species_host();
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
  HAT::t_int_1d ilist;
  if (reaxff->execution_space == HostKK) {
    NeighListKokkos<LMPHostType>* k_list = static_cast<NeighListKokkos<LMPHostType>*>(reaxff->list);
    k_list->k_ilist.sync_host();
    ilist = k_list->k_ilist.view_host();
  } else {
    NeighListKokkos<LMPDeviceType>* k_list = static_cast<NeighListKokkos<LMPDeviceType>*>(reaxff->list);
    k_list->k_ilist.sync_host();
    ilist = k_list->k_ilist.view_host();
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

          if ((j == 0) || (j < i)) continue;
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
    if (looptot >= 400*comm->nprocs) break;

  }
}

/* ----------------------------------------------------------------------
   deleting molecules removes atoms from the host copy of all per-atom
   data with AtomVec::copy(), so all of it must be current on the host
   before, and the device copy must be updated after atoms were removed
------------------------------------------------------------------------- */

void FixReaxFFSpeciesKokkos::DeleteSpecies(int Nmole, int Nspec)
{
  atomKK->sync(Host, ALL_MASK);
  const int nlocal_old = atom->nlocal;

  FixReaxFFSpecies::DeleteSpecies(Nmole, Nspec);

  if (atom->nlocal != nlocal_old) atomKK->modified(Host, ALL_MASK);
}
