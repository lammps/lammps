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

#include "fix_update_special_bonds.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixUpdateSpecialBonds::FixUpdateSpecialBonds(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix update/special/bonds command");

  restart_global = 1;
}

/* ---------------------------------------------------------------------- */

int FixUpdateSpecialBonds::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::setup(int /*vflag*/)
{
  // error if more than one fix update/special/bonds

  if (modify->get_fix_by_style("UPDATE_SPECIAL_BONDS").size() > 1)
    error->all(FLERR, "More than one fix update/special/bonds");

  // Require atoms know about all of their bonds and if they break
  if (force->newton_bond) error->all(FLERR, "Fix update/special/bonds requires Newton bond off");

  if (!atom->avec->bonds_allow)
    error->all(FLERR, "Fix update/special/bonds requires an atom style supporting bonds");

  // special lj must be 0 1 1 to censor pair forces between bonded particles
  // special coulomb must be 1 1 1 to ensure all pairs are included in the
  //   neighbor list and 1-3 and 1-4 special bond lists are skipped
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR, "Fix update/special/bonds requires special LJ weights = 0,1,1");
  if (force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0)
    error->all(FLERR, "Fix update/special/bonds requires special Coulomb weights = 1,1,1");
  // Implies neighbor->special_flag = [X, 2, 1, 1]
}

/* ----------------------------------------------------------------------
  Update special bond list and atom bond arrays, empty broken/created lists
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_exchange()
{
  int i, j, m, n1;
  tagint tagi, tagj;
  int nlocal = atom->nlocal;

  tagint *slist;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  for (auto const &it : broken_pairs) {
    tagi = it.first;
    tagj = it.second;
    i = atom->map(tagi);
    j = atom->map(tagj);

    // remove i from special bond list for atom j and vice versa
    // ignore n2, n3 since 1-3, 1-4 special factors required to be 1.0
    if (i < nlocal) {
      slist = special[i];
      n1 = nspecial[i][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagj) break;
      if (m == n1) error->one(FLERR, "Special bond {} {} not found", tagi, tagj);
      for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
      nspecial[i][0]--;
      nspecial[i][1] = nspecial[i][2] = nspecial[i][0];
    }

    if (j < nlocal) {
      slist = special[j];
      n1 = nspecial[j][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagi) break;
      if (m == n1) error->one(FLERR, "Special bond {} {} not found", tagi, tagj);
      for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
      nspecial[j][0]--;
      nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
    }
  }

  for (auto const &it : created_pairs) {
    tagi = it.first;
    tagj = it.second;
    i = atom->map(tagi);
    j = atom->map(tagj);

    // add i to special bond list for atom j and vice versa
    // ignore n2, n3 since 1-3, 1-4 special factors required to be 1.0
    n1 = nspecial[i][0];
    if (n1 >= atom->maxspecial)
      error->one(FLERR, "Special list size exceeded in fix update/special/bond");
    special[i][n1] = tagj;
    nspecial[i][0] += 1;
    nspecial[i][1] = nspecial[i][2] = nspecial[i][0];

    n1 = nspecial[j][0];
    if (n1 >= atom->maxspecial)
      error->one(FLERR, "Special list size exceeded in fix update/special/bond");
    special[j][n1] = tagi;
    nspecial[j][0] += 1;
    nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
  }

  broken_pairs.clear();
  created_pairs.clear();
}

/* ----------------------------------------------------------------------
  Update special lists for recently broken/created bonds
  Assumes appropriate atom/bond arrays were updated, e.g. had called
      neighbor->add_temporary_bond(i1, i2, btype);
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_force(int /*vflag*/)
{
  int i1, i2, j, jj, jnum;
  int *jlist, *numneigh, **firstneigh;
  tagint tag1, tag2;
  NeighList *list;

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  // In theory could communicate a list of broken bonds to neighboring processors here
  // to remove restriction that users use Newton bond off

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (auto const &it : new_broken_pairs) {
      tag1 = it.first;
      tag2 = it.second;
      i1 = atom->map(tag1);
      i2 = atom->map(tag2);

      // Loop through atoms of owned atoms i j
      if (i1 < nlocal) {
        jlist = firstneigh[i1];
        jnum = numneigh[i1];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag2) jlist[jj] = j;
        }
      }

      if (i2 < nlocal) {
        jlist = firstneigh[i2];
        jnum = numneigh[i2];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag1) jlist[jj] = j;
        }
      }
    }
  }

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    for (auto const &it : new_created_pairs) {
      tag1 = it.first;
      tag2 = it.second;
      i1 = atom->map(tag1);
      i2 = atom->map(tag2);

      // Loop through atoms of owned atoms i j and update SB bits
      if (i1 < nlocal) {
        jlist = firstneigh[i1];
        jnum = numneigh[i1];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag2) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }

      if (i2 < nlocal) {
        jlist = firstneigh[i2];
        jnum = numneigh[i2];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag1) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }
    }
  }

  new_broken_pairs.clear();
  new_created_pairs.clear();
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::add_broken_bond(int i, int j)
{
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_broken_pairs.push_back(tag_pair);
  broken_pairs.push_back(tag_pair);
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::add_created_bond(int i, int j)
{
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_created_pairs.push_back(tag_pair);
  created_pairs.push_back(tag_pair);
}

/* ----------------------------------------------------------------------
   Use write_restart to invoke pre_exchange
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::write_restart(FILE *fp)
{
  // Call pre-exchange to process any broken/created bonds

  pre_exchange();
  if (comm->me == 0) {
    int size = 0;
    fwrite(&size, sizeof(int), 1, fp);
  }
}
