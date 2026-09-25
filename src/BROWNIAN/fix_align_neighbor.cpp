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

#include "fix_align_neighbor.h"

#include "atom.h"
#include "error.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { DIPOLE };
enum { POLAR, NEMATIC };

/* ---------------------------------------------------------------------- */

FixAlignNeighbor::FixAlignNeighbor(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), symmetry(POLAR), ilevel_respa(0), list(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix align/neighbor", error);

  respa_level_support = 1;

  if (strcmp(arg[3], "dipole") == 0) {
    mode = DIPOLE;
  } else {
    error->all(FLERR, 3, "Unknown fix align/neighbor mode {}", arg[3]);
  }

  magnitude = utils::numeric(FLERR, arg[4], false, lmp);

  cutoff = utils::numeric(FLERR, arg[5], false, lmp);
  if (cutoff <= 0.0) error->all(FLERR, 5, "Fix align/neighbor cutoff must be > 0");
  cutsq = cutoff * cutoff;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "symmetry") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix align/neighbor symmetry", error);
      if (strcmp(arg[iarg + 1], "polar") == 0) {
        symmetry = POLAR;
      } else if (strcmp(arg[iarg + 1], "nematic") == 0) {
        symmetry = NEMATIC;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix align/neighbor symmetry {}", arg[iarg + 1]);
      }
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix align/neighbor keyword {}", arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixAlignNeighbor::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAlignNeighbor::init()
{
  if ((mode == DIPOLE) && (!atom->mu_flag || !atom->torque_flag))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix align/neighbor with mode dipole requires atom attributes mu and torque");

  // full neighbor list with the same fixed cutoff for all atom type pairs.
  // the cutoff may be larger than the pair style cutoff; the neighbor list
  // code checks that it fits within the communication cutoff.

  neighbor->add_request(this, NeighConst::REQ_FULL)->set_cutoff_fixed(cutoff);

  if (utils::strmatch(update->integrate_style, "^respa")) {
    int max_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    ilevel_respa = max_respa;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAlignNeighbor::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixAlignNeighbor::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    post_force(vflag);
  } else {
    auto *respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    respa->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAlignNeighbor::post_force(int /*vflag*/)
{
  if (mode == DIPOLE) post_force_dipole();
}

/* ---------------------------------------------------------------------- */

void FixAlignNeighbor::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ----------------------------------------------------------------------
   torque on atom i from all atoms j of the group within the cutoff:
     polar:   K * sum_j (e_i x e_j)
     nematic: K * sum_j (e_i . e_j) (e_i x e_j)
   with e = dipole vector.  the torques of a pair are equal and opposite.
------------------------------------------------------------------------- */

void FixAlignNeighbor::post_force_dipole()
{
  double **x = atom->x;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *mask = atom->mask;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    const double *ui = mu[i];
    double tx = 0.0, ty = 0.0, tz = 0.0;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;

      const double delx = x[i][0] - x[j][0];
      const double dely = x[i][1] - x[j][1];
      const double delz = x[i][2] - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq) continue;

      const double *uj = mu[j];
      double prefactor = magnitude;
      if (symmetry == NEMATIC) prefactor *= ui[0] * uj[0] + ui[1] * uj[1] + ui[2] * uj[2];

      tx += prefactor * (ui[1] * uj[2] - ui[2] * uj[1]);
      ty += prefactor * (ui[2] * uj[0] - ui[0] * uj[2]);
      tz += prefactor * (ui[0] * uj[1] - ui[1] * uj[0]);
    }

    torque[i][0] += tx;
    torque[i][1] += ty;
    torque[i][2] += tz;
  }
}
