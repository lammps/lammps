/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
  Contributed by: Jeremy Fersula @ Sorbonne University
----------------------------------------------------------------------- */

#include "fix_align_self.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { DIPOLE, QUAT };

static constexpr double SMALL = 1.0e-14;

/* ---------------------------------------------------------------------- */

FixAlignSelf::FixAlignSelf(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), avec(nullptr)
{

  if (narg != 5 && narg != 9) error->all(FLERR, "Incorrect number of fix align/self arguments");

  if (strcmp(arg[3], "dipole") == 0) {
    mode = DIPOLE;
  } else if (strcmp(arg[3], "quat") == 0) {
    mode = QUAT;
  } else {
    error->all(FLERR, 3, "Unknown fix align/self keyword {}", arg[3]);
  }

  magnitude = utils::numeric(FLERR, arg[4], false, lmp);

  // check for keyword

  if (narg == 9) {
    if (mode != QUAT)
      error->all(FLERR, 3, "Incorrect number of arguments for 'quat' mode of fix align/self");
    if (strcmp(arg[5], "qvector") == 0) {
      sx = utils::numeric(FLERR, arg[6], false, lmp);
      sy = utils::numeric(FLERR, arg[7], false, lmp);
      sz = utils::numeric(FLERR, arg[8], false, lmp);
      double snorm = sqrt(sx * sx + sy * sy + sz * sz);
      if (snorm < SMALL)
        error->all(FLERR, 5, "Fix align/self qvector magnitude {} is too small", snorm);
      sx = sx / snorm;
      sy = sy / snorm;
      sz = sz / snorm;
    } else {
      error->all(FLERR, 5, "Unknown fix align/self keyword {}", arg[5]);
    }
  } else {
    sx = 1.0;
    sy = 0.0;
    sz = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

int FixAlignSelf::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAlignSelf::init()
{
  if (mode == DIPOLE && (!atom->mu_flag || !atom->torque_flag))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix align/self with option dipole requires atom attributes mu + torque");

  if (mode == QUAT) {
    avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
    if (!avec)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix align/self with option quat requires atom style ellipsoid");

    // check that all particles are finite-size ellipsoids
    // no point particles allowed, spherical is OK

    int *ellipsoid = atom->ellipsoid;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (ellipsoid[i] < 0)
          error->one(FLERR, Error::NOLASTLINE,
                     "Fix align/self with option quat requires extended particles");
  }
}

/* ---------------------------------------------------------------------- */

void FixAlignSelf::post_force(int vflag)
{
  if (mode == DIPOLE)
    post_force_dipole(vflag);
  else if (mode == QUAT)
    post_force_quaternion(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAlignSelf::post_force_dipole(int /*vflag*/)
{
  double **torque = atom->torque;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **mu = atom->mu;
  double selfTorque[3];

  // Add the active torque to the atom torques:
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      MathExtra::cross3(mu[i], v[i], selfTorque);

      torque[i][0] += selfTorque[0] * magnitude;
      torque[i][1] += selfTorque[1] * magnitude;
      torque[i][2] += selfTorque[2] * magnitude;
    }
}

/* ---------------------------------------------------------------------- */

void FixAlignSelf::post_force_quaternion(int /*vflag*/)
{
  double **torque = atom->torque;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *ellipsoid = atom->ellipsoid;

  // ellipsoidal properties
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double f_act[3] = {sx, sy, sz};
  double f_rot[3];
  double *quat;
  double Q[3][3];
  double selfTorque[3];

  // Add the active torque to the atom torques:
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {

      quat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat(quat, Q);
      MathExtra::matvec(Q, f_act, f_rot);

      MathExtra::cross3(f_rot, v[i], selfTorque);

      torque[i][0] += selfTorque[0] * magnitude;
      torque[i][1] += selfTorque[1] * magnitude;
      torque[i][2] += selfTorque[2] * magnitude;
    }
  }
}
