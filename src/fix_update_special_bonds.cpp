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

#include "fix_update_special_bonds.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixUpdateSpecialBonds::FixUpdateSpecialBonds(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal fix update/special/bonds command");
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
  // Require atoms know about all of their bonds and if they break
  if (force->newton_bond) error->all(FLERR, "Fix update/special/bonds requires Newton bond off");

  if (!atom->avec->bonds_allow) error->all(FLERR, "Fix update/special/bonds requires atom bonds");

  // special lj must be 0 1 1 to censor pair forces between bonded particles
  // special coulomb must be 1 1 1 to ensure all pairs are included in the
  //   neighbor list and 1-3 and 1-4 special bond lists are skipped
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR, "Fix update/special/bonds requires special LJ weights = 0,1,1");
  if (force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0)
    error->all(FLERR, "Fix update/special/bonds requires special Coulomb weights = 1,1,1");

  new_broken_pairs.clear();
  broken_pairs.clear();
}

/* ----------------------------------------------------------------------
  Update special bond list and atom bond arrays, empty broken bond list
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_exchange()
{
  int i, j, m, n1, n3;
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
    if (i < nlocal) {
      slist = special[i];
      n1 = nspecial[i][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagj) break;
      n3 = nspecial[i][2];
      for (; m < n3 - 1; m++) slist[m] = slist[m + 1];
      nspecial[i][0]--;
      nspecial[i][1]--;
      nspecial[i][2]--;
    }

    if (j < nlocal) {
      slist = special[j];
      n1 = nspecial[j][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == tagi) break;
      n3 = nspecial[j][2];
      for (; m < n3 - 1; m++) slist[m] = slist[m + 1];
      nspecial[j][0]--;
      nspecial[j][1]--;
      nspecial[j][2]--;
    }
  }

  broken_pairs.clear();
}

/* ----------------------------------------------------------------------
  Loop neighbor list and update special bond lists for recently broken bonds
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_force(int /*vflag*/)
{
  int i1, i2, j, jj, jnum;
  int *jlist, *numneigh, **firstneigh;
  tagint tag1, tag2;

  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  NeighList *list = force->pair->list;    // may need to be generalized to work with pair hybrid*
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // In theory could communicate a list of broken bonds to neighboring processors here
  // to remove restriction that users use Newton bond off

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
  new_broken_pairs.clear();
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::add_broken_bond(int i, int j)
{
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_broken_pairs.push_back(tag_pair);
  broken_pairs.push_back(tag_pair);
}
