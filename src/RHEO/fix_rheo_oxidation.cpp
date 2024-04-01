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
#include "compute_rheo_surface.h"
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
  Fix(lmp, narg, arg), compute_surface(nullptr), fix_rheo(nullptr)
{
  if (narg != 6) error->all(FLERR,"Illegal fix command");

  cut = utils::numeric(FLERR, arg[3], false, lmp);
  if (cut <= 0.0) error->all(FLERR, "Illegal bond cutoff {} in fix rheo/oxidation", cut);

  btype = utils::inumeric(FLERR, arg[4], false, lmp);
  if (btype < 1 || btype > atom->nbondtypes) error->all(FLERR, "Illegal value {} for bond type in fix rheo/oxidation", btype);

  rsurf = utils::numeric(FLERR, arg[5], false, lmp);
  if (rsurf <= 0.0) error->all(FLERR, "Illegal surface distance {} in fix rheo/oxidation", cut);

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
  mask |= PRE_FORCE;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/oxidation");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (cut > fix_rheo->h)
    error->all(FLERR, "Bonding length exceeds kernel cutoff");

  if (rsurf >= fix_rheo->h)
    error->all(FLERR, "Surface distance must be less than kernel cutoff");

  if (!force->bond) error->all(FLERR, "Must define a bond style with fix rheo/oxidation");
  if (!atom->avec->bonds_allow) error->all(FLERR, "Fix rheo/oxidation requires atom bonds");

  int tmp1, tmp2;
  index_nb = atom->find_custom("shell_nbond", tmp1, tmp2);
  if (index_nb == -1)
    error->all(FLERR, "Must use bond style rheo/shell to use fix rheo/oxidation");
  nbond = atom->ivector[index_nb];

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

void FixRHEOOxidation::setup_pre_force(int /*vflag*/)
{
  // Not strictly required that this fix be after FixRHEO,
  // but enforce to be consistent with other RHEO fixes
  fix_rheo->oxidation_fix_defined = 1;

  if (!fix_rheo->surface_flag) error->all(FLERR,
      "fix rheo/oxidation requires surface calculation in fix rheo");
  compute_surface = fix_rheo->compute_surface;

  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEOOxidation::pre_force(int /*vflag*/)
{
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
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  double *rsurface = compute_surface->rsurface;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (rsurface[i] > rsurf) continue;

    tagi = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (rsurface[j] > rsurf) continue;

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
