/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Originally modified from CG-DNA/fix_nve_dotc_langevin.cpp.

   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include "fix_brownian_dipole.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBrownianDipole::FixBrownianDipole(LAMMPS *lmp, int narg, char **arg) :
    FixBrownianBase(lmp, narg, arg)
{
  if (gamma_t_eigen_flag || gamma_r_eigen_flag) {
    error->all(FLERR, "Illegal fix brownian command.");
  }

  if (!gamma_t_flag || !gamma_r_flag) { error->all(FLERR, "Illegal fix brownian command."); }
  if (!atom->mu_flag) error->all(FLERR, "Fix brownian/fluct_dipole requires atom attribute mu");
  if (!atom->sphere_flag) error->all(FLERR, "Fix brownian/fluct_dipole requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

void FixBrownianDipole::init()
{
  FixBrownianBase::init();

  g3 = g1 / gamma_r;
  g4 = g2 * sqrt(rot_temp / gamma_r);
  g1 /= gamma_t;
  g2 *= sqrt(temp / gamma_t);
}

/* ---------------------------------------------------------------------- */

void FixBrownianDipole::initial_integrate(int /*vflag */)
{
  if (domain->dimension == 2) {
    if (!noise_flag) {
      initial_integrate_templated<0, 0, 1, 0>();
    } else if (gaussian_noise_flag) {
      initial_integrate_templated<0, 1, 1, 0>();
    } else {
      initial_integrate_templated<1, 0, 1, 0>();
    }
  } else if (planar_rot_flag) {
    if (!noise_flag) {
      initial_integrate_templated<0, 0, 0, 1>();
    } else if (gaussian_noise_flag) {
      initial_integrate_templated<0, 1, 0, 1>();
    } else {
      initial_integrate_templated<1, 0, 0, 1>();
    }
  } else {
    if (!noise_flag) {
      initial_integrate_templated<0, 0, 0, 0>();
    } else if (gaussian_noise_flag) {
      initial_integrate_templated<0, 1, 0, 0>();
    } else {
      initial_integrate_templated<1, 0, 0, 0>();
    }
  }
  }

/* ---------------------------------------------------------------------- */

template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D, int Tp_2Drot>
void FixBrownianDipole::initial_integrate_templated()
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double wx, wy, wz;
  double **torque = atom->torque;
  double **mu = atom->mu;
  double mulensq;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dx, dy, dz;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (Tp_2D) {
        dz = 0;
	wz = 0;
        if (Tp_UNIFORM) {
          dx = dt * (g1 * f[i][0] + g2 * (rng->uniform() - 0.5));
          dy = dt * (g1 * f[i][1] + g2 * (rng->uniform() - 0.5));
	  wx = (rng->uniform() - 0.5) * g4;
	  wy = (rng->uniform() - 0.5) * g4;

        } else if (Tp_GAUSS) {
          dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
          dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
          wx = rng->gaussian() * g4;
          wy = rng->gaussian() * g4;
        } else {
          dx = dt * g1 * f[i][0];
          dy = dt * g1 * f[i][1];
          wx = 0;
	  wy = 0;
        }
      } else if (Tp_2Drot) {
        wz = 0;
        if (Tp_UNIFORM) {
          dx = dt * (g1 * f[i][0] + g2 * (rng->uniform() - 0.5));
          dy = dt * (g1 * f[i][1] + g2 * (rng->uniform() - 0.5));
          dz = dt * (g1 * f[i][2] + g2 * (rng->uniform() - 0.5));
          wx = (rng->uniform() - 0.5) * g4;
          wy = (rng->uniform() - 0.5) * g4;
        } else if (Tp_GAUSS) {
          dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
          dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
          dz = dt * (g1 * f[i][2] + g2 * rng->gaussian());
          wx = rng->gaussian() * g4;
          wy = rng->gaussian() * g4;
        } else {
          dx = dt * g1 * f[i][0];
          dy = dt * g1 * f[i][1];
          dz = dt * g1 * f[i][2];
          wx = 0;
          wy = 0;
        }
      } else {
        if (Tp_UNIFORM) {
          dx = dt * (g1 * f[i][0] + g2 * (rng->uniform() - 0.5));
          dy = dt * (g1 * f[i][1] + g2 * (rng->uniform() - 0.5));
          dz = dt * (g1 * f[i][2] + g2 * (rng->uniform() - 0.5));
          wx = (rng->uniform() - 0.5) * g4;
          wy = (rng->uniform() - 0.5) * g4;
          wz = (rng->uniform() - 0.5) * g4;
        } else if (Tp_GAUSS) {
          dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
          dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
          dz = dt * (g1 * f[i][2] + g2 * rng->gaussian());
          wx = rng->gaussian() * g4;
          wy = rng->gaussian() * g4;
          wz = rng->gaussian() * g4;
        } else {
          dx = dt * g1 * f[i][0];
          dy = dt * g1 * f[i][1];
          dz = dt * g1 * f[i][2];
          wx = wy = wz = 0;
        }
      }

      x[i][0] += dx;
      v[i][0] = dx / dt;

      x[i][1] += dy;
      v[i][1] = dy / dt;

      x[i][2] += dz;
      v[i][2] = dz / dt;

      wx += g3 * torque[i][0];
      wy += g3 * torque[i][1];
      wz += g3 * torque[i][2];

      mulensq = mu[i][0] * mu[i][0] + mu[i][1] * mu[i][1] + mu[i][2] * mu[i][2];

      mu[i][0] += (mu[i][0]*lambda*(1-mulensq) + wx)*dt;
      mu[i][1] += (mu[i][1]*lambda*(1-mulensq) + wy)*dt;
      mu[i][2] += (mu[i][2]*lambda*(1-mulensq) + wz)*dt;

      
    }
  }
}
