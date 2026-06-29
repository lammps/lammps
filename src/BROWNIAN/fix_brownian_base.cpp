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

/* ----------------------------------------------------------------------
   Originally modified from CG-DNA/fix_nve_dotc_langevin.cpp.

   Contributing authors: Sam Cameron (University of Bristol),
                         Arthur Straube (Zuse Institute Berlin)
------------------------------------------------------------------------- */

#include "fix_brownian_base.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
FixBrownianBase::FixBrownianBase(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), gamma_t_inv(nullptr), gamma_r_inv(nullptr), gamma_t_invsqrt(nullptr),
    gamma_r_invsqrt(nullptr), dipole_body(nullptr), rng(nullptr)
{
  time_integrate = 1;

  noise_flag = 1;
  gaussian_noise_flag = 0;
  gamma_t_flag = gamma_r_flag = 0;
  gamma_t_eigen_flag = gamma_r_eigen_flag = 0;
  dipole_flag = 0;
  rot_temp_flag = 0;
  planar_rot_flag = 0;
  rot_style = ROT_GEOMETRIC;
  if (utils::strmatch(style, "^brownian/sphere")) rot_style = ROT_PROJECTION;
  g2 = 0.0;

  std::string mystyle = fmt::format("fix {}", style);
  if (narg < 5) utils::missing_cmd_args(FLERR, mystyle, error);

  temp = utils::numeric(FLERR, arg[3], false, lmp);
  if (temp <= 0) error->all(FLERR, 3, "Fix {} temp must be > 0.0", style);

  seed = utils::inumeric(FLERR, arg[4], false, lmp);
  if (seed <= 0) error->all(FLERR, 4, "Fix {} seed must be > 0", style);

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "rng") == 0) {
      if (narg < iarg + 2) utils::missing_cmd_args(FLERR, mystyle + " rng", error);
      if (strcmp(arg[iarg + 1], "uniform") == 0) {
        noise_flag = 1;
      } else if (strcmp(arg[iarg + 1], "gaussian") == 0) {
        noise_flag = 1;
        gaussian_noise_flag = 1;
      } else if (strcmp(arg[iarg + 1], "none") == 0) {
        noise_flag = 0;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix {} rng keyword {}", style, arg[iarg + 1]);
      }
      iarg = iarg + 2;
    } else if (strcmp(arg[iarg], "dipole") == 0) {
      if (narg < iarg + 4) utils::missing_cmd_args(FLERR, mystyle + " dipole", error);

      dipole_flag = 1;
      delete[] dipole_body;
      dipole_body = new double[3];

      dipole_body[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      dipole_body[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      dipole_body[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg = iarg + 4;

    } else if (strcmp(arg[iarg], "gamma_t_eigen") == 0) {
      if (narg < iarg + 4) utils::missing_cmd_args(FLERR, mystyle + " gamma_t_eigen", error);

      double gamma_t_tmp[3];
      gamma_t_tmp[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      gamma_t_tmp[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (domain->dimension == 2) {
        if (strcmp(arg[iarg + 3], "inf") != 0)
          error->all(FLERR, iarg + 3, "Fix {} gamma_t_eigen third value must be inf for 2D system",
                     style);
        gamma_t_tmp[2] = 1.0;
      } else {
        gamma_t_tmp[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      }
      if ((gamma_t_tmp[0] <= 0.0) || (gamma_t_tmp[1] <= 0.0) || (gamma_t_tmp[2] <= 0.0))
        error->all(FLERR, iarg, "Fix {} gamma_t_eigen values must be > 0", style);

      gamma_t_eigen_flag = 1;
      delete[] gamma_t_inv;
      delete[] gamma_t_invsqrt;
      gamma_t_inv = new double[3];
      gamma_t_invsqrt = new double[3];
      gamma_t_inv[0] = 1.0 / gamma_t_tmp[0];
      gamma_t_inv[1] = 1.0 / gamma_t_tmp[1];
      gamma_t_inv[2] = 1.0 / gamma_t_tmp[2];
      if (domain->dimension == 2) gamma_t_inv[2] = 0.0;

      gamma_t_invsqrt[0] = sqrt(gamma_t_inv[0]);
      gamma_t_invsqrt[1] = sqrt(gamma_t_inv[1]);
      gamma_t_invsqrt[2] = sqrt(gamma_t_inv[2]);
      iarg = iarg + 4;

    } else if (strcmp(arg[iarg], "gamma_r_eigen") == 0) {
      if (narg < iarg + 4) utils::missing_cmd_args(FLERR, mystyle + " gamma_r_eigen", error);

      double gamma_r_tmp[3];
      if (domain->dimension == 2) {
        if (strcmp(arg[iarg + 1], "inf") != 0)
          error->all(FLERR, iarg + 1, "Fix {} gamma_r_eigen first value must be inf for 2D system",
                     style);
        gamma_r_tmp[0] = 1.0;

        if (strcmp(arg[iarg + 2], "inf") != 0)
          error->all(FLERR, iarg + 2, "Fix {} gamma_r_eigen second value must be inf for 2D system",
                     style);
        gamma_r_tmp[1] = 1.0;
      } else {
        gamma_r_tmp[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        gamma_r_tmp[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      }
      gamma_r_tmp[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);

      if ((gamma_r_tmp[0] <= 0.0) || (gamma_r_tmp[1] <= 0) || (gamma_r_tmp[2] <= 0))
        error->all(FLERR, iarg, "Fix {} gamma_r_eigen values must be > 0", style);

      gamma_r_eigen_flag = 1;
      delete[] gamma_r_inv;
      delete[] gamma_r_invsqrt;
      gamma_r_inv = new double[3];
      gamma_r_invsqrt = new double[3];

      gamma_r_inv[0] = 1.0 / gamma_r_tmp[0];
      gamma_r_inv[1] = 1.0 / gamma_r_tmp[1];
      gamma_r_inv[2] = 1.0 / gamma_r_tmp[2];
      if (domain->dimension == 2) {
        gamma_r_inv[0] = 0.0;
        gamma_r_inv[1] = 0.0;
      }
      gamma_r_invsqrt[0] = sqrt(gamma_r_inv[0]);
      gamma_r_invsqrt[1] = sqrt(gamma_r_inv[1]);
      gamma_r_invsqrt[2] = sqrt(gamma_r_inv[2]);
      iarg = iarg + 4;

    } else if (strcmp(arg[iarg], "gamma_t") == 0) {
      if (narg < iarg + 2) utils::missing_cmd_args(FLERR, mystyle + " gamma_t", error);

      gamma_t_flag = 1;
      gamma_t = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (gamma_t <= 0.0) error->all(FLERR, iarg + 1, "Fix {} gamma_t value must be > 0", style);
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "gamma_r") == 0) {
      if (narg < iarg + 2) utils::missing_cmd_args(FLERR, mystyle + " gamma_r", error);

      gamma_r_flag = 1;
      gamma_r = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (gamma_r <= 0) error->all(FLERR, iarg + 1, "Fix {} gamma_r value must be > 0", style);
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "rotation_temp") == 0) {
      if (narg < iarg + 2) utils::missing_cmd_args(FLERR, mystyle + " rotation_temp", error);

      rot_temp_flag = 1;
      rot_temp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (rot_temp <= 0)
        error->all(FLERR, iarg + 1, "Fix {} rotation_temp value must be > 0", style);
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "planar_rotation") == 0) {

      planar_rot_flag = 1;
      if (domain->dimension == 2)
        error->all(FLERR, iarg, "The planar_rotation keyword is not allowed for 2D simulations");
      ++iarg;

    } else if (strcmp(arg[iarg], "rotation_style") == 0) {

      if (!utils::strmatch(style, "^brownian/sphere"))
        error->all(FLERR, "Keyword rotation_style is only supported for fix brownian/sphere");
      if (narg < iarg + 2) utils::missing_cmd_args(FLERR, "fix brownian rotation_style", error);

      if (strcmp(arg[iarg + 1], "projection") == 0)
        rot_style = ROT_PROJECTION;
      else if (strcmp(arg[iarg + 1], "geometric") == 0)
        rot_style = ROT_GEOMETRIC;
      else
        error->all(FLERR, iarg + 1, "Unknown fix {} rotation_style {}", style);
      iarg = iarg + 2;

    } else {
      error->all(FLERR, iarg, "Unknown fix {} keyword {}", style, arg[iarg]);
    }
  }
  if (!rot_temp_flag) rot_temp = temp;

  // initialize Marsaglia RNG with processor-unique seed
  rng = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

int FixBrownianBase::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBrownianBase::~FixBrownianBase()
{

  if (gamma_t_eigen_flag) {
    delete[] gamma_t_inv;
    delete[] gamma_t_invsqrt;
  }
  if (gamma_r_eigen_flag) {
    delete[] gamma_r_inv;
    delete[] gamma_r_invsqrt;
  }

  if (dipole_flag) delete[] dipole_body;
  delete rng;
}

/* ---------------------------------------------------------------------- */

void FixBrownianBase::init()
{
  dt = update->dt;
  sqrtdt = sqrt(dt);
  g1 = force->ftm2v;
  if (noise_flag == 0) {
    g2 = 0.0;
  } else if (gaussian_noise_flag == 1) {
    g2 = sqrt(2 * force->boltz / dt / force->mvv2e);
  } else {
    g2 = sqrt(24 * force->boltz / dt / force->mvv2e);
  }
}

void FixBrownianBase::reset_dt()
{
  double sqrtdt_old = sqrtdt;
  dt = update->dt;
  sqrtdt = sqrt(dt);
  g2 *= sqrtdt_old / sqrtdt;
}
