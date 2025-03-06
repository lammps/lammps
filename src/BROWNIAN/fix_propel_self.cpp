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
   Contributed by Stefan Paquay @ Brandeis University

   Thanks to Liesbeth Janssen @ Eindhoven University for useful discussions!

   Current maintainer: Sam Cameron @ University of Bristol
----------------------------------------------------------------------- */

#include "fix_propel_self.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { DIPOLE, VELOCITY, QUAT };

static constexpr double TOL = 1e-14;

/* ---------------------------------------------------------------------- */

FixPropelSelf::FixPropelSelf(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), avec(nullptr)
{

  virial_global_flag = virial_peratom_flag = 1;

  if (narg != 5 && narg != 9) error->all(FLERR, "Illegal fix propel/self command");

  if (strcmp(arg[3], "velocity") == 0) {
    mode = VELOCITY;
    thermo_virial = 0;
  } else if (strcmp(arg[3], "dipole") == 0) {
    mode = DIPOLE;
    thermo_virial = 1;
  } else if (strcmp(arg[3], "quat") == 0) {
    mode = QUAT;
    thermo_virial = 1;
  } else {
    error->all(FLERR, "Illegal fix propel/self command");
  }

  magnitude = utils::numeric(FLERR, arg[4], false, lmp);

  // check for keyword

  if (narg == 9) {
    if (mode != QUAT) { error->all(FLERR, "Illegal fix propel/self command"); }
    if (strcmp(arg[5], "qvector") == 0) {
      sx = utils::numeric(FLERR, arg[6], false, lmp);
      sy = utils::numeric(FLERR, arg[7], false, lmp);
      sz = utils::numeric(FLERR, arg[8], false, lmp);
      double snorm = sqrt(sx * sx + sy * sy + sz * sz);
      sx = sx / snorm;
      sy = sy / snorm;
      sz = sz / snorm;
    } else {
      error->all(FLERR, "Illegal fix propel/self command");
    }
  } else {
    sx = 1.0;
    sy = 0.0;
    sz = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

int FixPropelSelf::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::init()
{
  if (mode == DIPOLE && !atom->mu_flag)
    error->all(FLERR, "Fix propel/self requires atom attribute mu with option dipole");

  if (mode == QUAT) {
    avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
    if (!avec) error->all(FLERR, "Fix propel/self requires atom style ellipsoid with option quat");

    // check that all particles are finite-size ellipsoids
    // no point particles allowed, spherical is OK

    int *ellipsoid = atom->ellipsoid;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (ellipsoid[i] < 0)
          error->one(FLERR, "Fix propel/self requires extended particles with option quat");
  }
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::post_force(int vflag)
{
  if (mode == DIPOLE)
    post_force_dipole(vflag);
  else if (mode == VELOCITY)
    post_force_velocity(vflag);
  else if (mode == QUAT)
    post_force_quaternion(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::post_force_dipole(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **mu = atom->mu;
  double fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      fx = magnitude * mu[i][0];
      fy = magnitude * mu[i][1];
      fz = magnitude * mu[i][2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::post_force_velocity(int vflag)
{
  double **f = atom->f;
  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double nv2, fnorm, fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {

      nv2 = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
      fnorm = 0.0;

      if (nv2 > TOL) {

        // Without this check you can run into numerical
        // issues because fnorm will blow up.

        fnorm = magnitude / sqrt(nv2);
      }
      fx = fnorm * v[i][0];
      fy = fnorm * v[i][1];
      fz = fnorm * v[i][2];

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::post_force_quaternion(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *ellipsoid = atom->ellipsoid;

  // ellipsoidal properties
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double f_act[3] = {sx, sy, sz};
  double f_rot[3];
  double *quat;
  double Q[3][3];
  double fx, fy, fz;

  // energy and virial setup
  double vi[6];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  // if domain has PBC, need to unwrap for virial
  double unwrap[3];
  imageint *image = atom->image;

  // Add the active force to the atom force:
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {

      quat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat(quat, Q);
      MathExtra::matvec(Q, f_act, f_rot);

      fx = magnitude * f_rot[0];
      fy = magnitude * f_rot[1];
      fz = magnitude * f_rot[2];

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      if (evflag) {
        domain->unmap(x[i], image[i], unwrap);
        vi[0] = fx * unwrap[0];
        vi[1] = fy * unwrap[1];
        vi[2] = fz * unwrap[2];
        vi[3] = fx * unwrap[1];
        vi[4] = fx * unwrap[2];
        vi[5] = fy * unwrap[2];
        v_tally(i, vi);
      }
    }
  }
}
