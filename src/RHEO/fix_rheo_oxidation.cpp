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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "fix_rheo_oxidation.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
enum {NONE, CONSTANT};

/* ---------------------------------------------------------------------- */

FixRHEOOxidation::FixRHEOOxidation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) //, fix_bond_history(nullptr)
{
  if (narg != 5) error->all(FLERR,"Illegal fix command");

  cut = utils::numeric(FLERR, arg[3], false, lmp);
  if (cut <= 0.0) error->all(FLERR, "Illegal bond cutoff {} in fix rheo/oxidation", cut);

  btype = utils::inumeric(FLERR, arg[4], false, lmp);
  if (btype < 1 || btype > atom->nbondtypes) error->all(FLERR, "Illegal value {} for bond type in fix rheo/oxidation", btype);

  cutsq = cut * cut;
}

/* ---------------------------------------------------------------------- */

FixRHEOOxidation::~FixRHEOOxidation()
{
}

/* ---------------------------------------------------------------------- */

int FixRHEOOxidation::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/oxidation");
  class FixRHEO *fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);
  double cut_kernel = fix_rheo->h;

  if (cut > cut_kernel)
    error->all(FLERR, "Bonding length exceeds kernel cutoff");

  if (!force->bond) error->all(FLERR, "Must define a bond style with fix rheo/oxidation");
  if (!atom->avec->bonds_allow) error->all(FLERR, "Fix rheo/oxidation requires atom bonds");

  //// find instances of bond history to delete data
  //histories = modify->get_fix_by_style("BOND_HISTORY");
  //for (auto &ihistory: histories)
  //  if (strcmp(histories[i]->id, "HISTORY_RHEO_SHELL") == 0)
  //    fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
//
  //if (!fix_bond_history)
  //  error->all(FLERR, "Must define bond style rheo/shell to use fix rheo/oxidation");

  // need a half neighbor list
  auto req = neighbor->add_request(this, NeighConst::REQ_DEFAULT);
  req->set_cutoff(cut);
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::post_integrate()
{
  int i, j, n, ii, jj, inum, jnum, bflag;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  tagint tagi, tagj;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int *status = atom->status;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(status[i] & STATUS_SURFACE)) continue;

    tagi = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(status[j] & STATUS_SURFACE)) continue;

      tagj = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq) continue;

      // Check if already have an oxide bond
      bflag = 0;
      for (n = 0; n < num_bond[i]; n++) {
        if (bond_type[i][n] == btype && bond_atom[i][n] == tagj) {
          bflag = 1;
          break;
        }
      }
      if (bflag) continue;

      for (n = 0; n < num_bond[j]; n++) {
        if (bond_type[j][n] == btype && bond_atom[j][n] == tagi) {
          bflag = 1;
          break;
        }
      }
      if (bflag) continue;

      // Add bonds to owned atoms
      // If newton bond, add to both, otherwise add to whichever has a smaller tag
      if (i < nlocal && (!newton_bond || tagi < tagj)) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR,"New bond exceeded bonds per atom in fix rheo/oxidation for atom {}", tagi);
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tagj;
        num_bond[i]++;
      }

      if (j < nlocal && (!newton_bond || tagj < tagi)) {
        if (num_bond[j] == atom->bond_per_atom)
          error->one(FLERR,"New bond exceeded bonds per atom in fix rheo/oxidation for atom {}", tagj);
        bond_type[j][num_bond[j]] = btype;
        bond_atom[j][num_bond[j]] = tagi;
        num_bond[j]++;
      }
    }
  }
}
