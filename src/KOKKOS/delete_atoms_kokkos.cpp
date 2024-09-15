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
   Contributing author: Mitch Murphy, alphataubio at gmail
------------------------------------------------------------------------- */

#include "delete_atoms_kokkos.h"

#include "angle.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kokkos.h"
#include "kspace.h"
#include "neighbor_kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DeleteAtomsKokkos<DeviceType>::DeleteAtomsKokkos(LAMMPS *lmp) : DeleteAtoms(lmp)
{
  atomKK = (AtomKokkos *) atom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DeleteAtomsKokkos<DeviceType>::command(int narg, char **arg)
{
  atomKK->sync(Host, X_MASK|RMASS_MASK|TYPE_MASK);
  DeleteAtoms::command(narg, arg);
}

/* ----------------------------------------------------------------------
   delete atoms so there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

template<class DeviceType>
void DeleteAtomsKokkos<DeviceType>::delete_overlap(int narg, char **arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "delete_atoms overlap", error);

  // read args

  const double cut = utils::numeric(FLERR, arg[1], false, lmp);
  const double cutsq = cut * cut;

  int igroup1 = group->find(arg[2]);
  if (igroup1 < 0)
    error->all(FLERR, "Could not find delete_atoms overlap first group ID {}", arg[2]);
  int igroup2 = group->find(arg[3]);
  if (igroup2 < 0)
    error->all(FLERR, "Could not find delete_atoms overlap second group ID {}", arg[3]);
  options(narg - 4, &arg[4]);

  const int group1bit = group->bitmask[igroup1];
  const int group2bit = group->bitmask[igroup2];

  if (comm->me == 0) utils::logmesg(lmp, "System init for delete_atoms/kk ...\n");

  // request a full neighbor list for use by this command

  neighbor->add_request(this, "delete_atoms/kk", NeighConst::REQ_FULL);

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  lmp->init();

  // error check on cutoff
  // if no pair style, neighbor list will be empty

  if (force->pair == nullptr) error->all(FLERR, "Delete_atoms requires a pair style be defined");
  if (cut > neighbor->cutneighmax) error->all(FLERR, "Delete_atoms cutoff > max neighbor cutoff");
  if (cut > neighbor->cutneighmin && comm->me == 0)
    error->warning(FLERR, "Delete_atoms cutoff > minimum neighbor cutoff");

  // setup domain, communication and neighboring
  // acquire ghosts and build standard neighbor lists

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  neighbor->build(1);

  // build neighbor list this command needs based on the earlier request

  auto list = neighbor->find_list(this);
  neighbor->build_one(list);

  auto inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  auto d_numneigh = k_list->d_numneigh;
  auto d_neighbors = k_list->d_neighbors;
  auto d_ilist = k_list->d_ilist;

  // allocate and initialize deletion list
  // must be after exchange potentially changes nlocal

  int nlocal = atom->nlocal;
  memoryKK->create_kokkos(k_dlist, dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;
  k_dlist.template sync<DeviceType>();



  // double loop over owned atoms and their full neighbor list
  // at end of loop, there are no more overlaps
  // only ever delete owned atom I in I loop iteration, never J even if owned

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_tag = atomKK->k_tag.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int i, j, ii, jj, jnum;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double factor_lj, factor_coul;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  copymode = 1;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(d_mask[i] & (group1bit | group2bit))) continue;
    double xtmp = d_x(i,0);
    double ytmp = d_x(i,1);
    double ztmp = d_x(i,2);
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      if (!(d_mask[j] & (group1bit | group2bit))) continue;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      // only consider deletion if I,J distance < cutoff
      // compute rsq identically on both I,J loop iterations
      // ignoring possibility that I,J tags are equal

      double delx, dely, delz;

      if (d_tag(i) < d_tag(j)) {
        delx = xtmp - d_x(j,0);
        dely = ytmp - d_x(j,1);
        delz = ztmp - d_x(j,2);
      } else {
        delx = d_x(j,0) - xtmp;
        dely = d_x(j,1) - ytmp;
        delz = d_x(j,2) - ztmp;
      }
      double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq) continue;

      // only consider deletion if I,J are in groups 1,2 respectively
      // true whether J is owned or ghost atom

      if (!(d_mask[i] & group1bit)) continue;
      if (!(d_mask[j] & group2bit)) continue;

      // J is owned atom:
      //   delete atom I if atom J has not already been deleted
      // J is ghost atom:
      //   delete atom I if J,I is not a candidate deletion pair
      //     due to being in groups 1,2 respectively
      //   if they are candidate pair, then either:
      //      another proc owns J and could delete J
      //      J is a ghost of another of my owned atoms, and I could delete J
      //   test on tags of I,J ensures that only I or J is deleted

      if (j < nlocal) {
        if (dlist[j]) continue;
      } else if ((d_mask[i] & group2bit) && (d_mask[j] & group1bit)) {
        if (d_tag(i) > d_tag(j)) continue;
      }

      dlist[i] = 1;
      break;
    }
  }
  neighbor->init();
  k_dlist.template modify<DeviceType>();
}

namespace LAMMPS_NS {
template class DeleteAtomsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class DeleteAtomsKokkos<LMPHostType>;
#endif
}
