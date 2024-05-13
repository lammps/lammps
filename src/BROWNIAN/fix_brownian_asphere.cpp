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

   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include "fix_brownian_asphere.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBrownianAsphere::FixBrownianAsphere(LAMMPS *lmp, int narg, char **arg) :
    FixBrownianBase(lmp, narg, arg), avec(nullptr)
{
  if (!gamma_t_eigen_flag || !gamma_r_eigen_flag) {
    error->all(FLERR, "Illegal fix brownian command.");
  }

  if (gamma_t_flag || gamma_r_flag) error->all(FLERR, "Illegal fix brownian command.");

  if (dipole_flag && !atom->mu_flag)
    error->all(FLERR, "Fix brownian/asphere dipole requires atom attribute mu");

  if (!atom->ellipsoid_flag)
    error->all(FLERR, "Fix brownian/asphere requires atom style ellipsoid");

  if (planar_rot_flag && (comm->me == 0)) {
    error->warning(FLERR, "Ignoring first two entries of gamma_r_eigen since rotation is planar.");
  }
}

/* ---------------------------------------------------------------------- */

void FixBrownianAsphere::init()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec) error->all(FLERR, "Compute brownian/asphere requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0) error->one(FLERR, "Fix brownian/asphere requires extended particles");

  if (dipole_flag) {

    double f_rot[3];
    double *quat;
    AtomVecEllipsoid::Bonus *bonus = avec->bonus;

    double Q[3][3];
    double **mu = atom->mu;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        quat = bonus[ellipsoid[i]].quat;
        MathExtra::quat_to_mat(quat, Q);
        MathExtra::matvec(Q, dipole_body, f_rot);

        mu[i][0] = f_rot[0];
        mu[i][1] = f_rot[1];
        mu[i][2] = f_rot[2];
      }
    }
  }

  FixBrownianBase::init();

  g4 = g2 * sqrt(rot_temp);
  g2 *= sqrt(temp);
}

/* ---------------------------------------------------------------------- */

void FixBrownianAsphere::initial_integrate(int /*vflag */)
{
  if (domain->dimension == 2) {
    if (dipole_flag) {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 1, 1, 0>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 1, 1, 0>();
      } else {
        initial_integrate_templated<1, 0, 1, 1, 0>();
      }
    } else {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 0, 1, 0>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 0, 1, 0>();
      } else {
        initial_integrate_templated<1, 0, 0, 1, 0>();
      }
    }
  } else if (planar_rot_flag) {
    if (dipole_flag) {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 1, 0, 1>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 1, 0, 1>();
      } else {
        initial_integrate_templated<1, 0, 1, 0, 1>();
      }
    } else {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 0, 0, 1>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 0, 0, 1>();
      } else {
        initial_integrate_templated<1, 0, 0, 0, 1>();
      }
    }
  } else {
    if (dipole_flag) {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 1, 0, 0>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 1, 0, 0>();
      } else {
        initial_integrate_templated<1, 0, 1, 0, 0>();
      }
    } else {
      if (!noise_flag) {
        initial_integrate_templated<0, 0, 0, 0, 0>();
      } else if (gaussian_noise_flag) {
        initial_integrate_templated<0, 1, 0, 0, 0>();
      } else {
        initial_integrate_templated<1, 0, 0, 0, 0>();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_DIPOLE, int Tp_2D, int Tp_2Drot>
void FixBrownianAsphere::initial_integrate_templated()
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  double **mu = atom->mu;
  double **torque = atom->torque;
  double qw[4];
  double *quat;
  int *ellipsoid = atom->ellipsoid;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // project dipole along x axis of quat
  double f_rot[3];
  double rotationmatrix_transpose[3][3];
  double tmp[3];
  double dv[3];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      // update orientation first

      quat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat_trans(quat, rotationmatrix_transpose);

      // tmp holds angular velocity in body frame
      MathExtra::matvec(rotationmatrix_transpose, torque[i], tmp);

      if (Tp_2D) {
        tmp[0] = tmp[1] = 0.0;
        if (Tp_UNIFORM) {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * (rng->uniform() - 0.5) * g4;
        } else if (Tp_GAUSS) {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * rng->gaussian() * g4;
        } else {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2];
        }
      } else if (Tp_2Drot) {
        tmp[0] = tmp[1] = 0.0;
        if (Tp_UNIFORM) {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * (rng->uniform() - 0.5) * g4;
        } else if (Tp_GAUSS) {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * rng->gaussian() * g4;
        } else {
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2];
        }
      } else {
        if (Tp_UNIFORM) {
          tmp[0] = g1 * tmp[0] * gamma_r_inv[0] + gamma_r_invsqrt[0] * (rng->uniform() - 0.5) * g4;
          tmp[1] = g1 * tmp[1] * gamma_r_inv[1] + gamma_r_invsqrt[1] * (rng->uniform() - 0.5) * g4;
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * (rng->uniform() - 0.5) * g4;
        } else if (Tp_GAUSS) {
          tmp[0] = g1 * tmp[0] * gamma_r_inv[0] + gamma_r_invsqrt[0] * rng->gaussian() * g4;
          tmp[1] = g1 * tmp[1] * gamma_r_inv[1] + gamma_r_invsqrt[1] * rng->gaussian() * g4;
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2] + gamma_r_invsqrt[2] * rng->gaussian() * g4;
        } else {
          tmp[0] = g1 * tmp[0] * gamma_r_inv[0];
          tmp[1] = g1 * tmp[1] * gamma_r_inv[1];
          tmp[2] = g1 * tmp[2] * gamma_r_inv[2];
        }
      }

      // convert body frame angular velocity to quaternion
      MathExtra::quatvec(quat, tmp, qw);
      quat[0] = quat[0] + 0.5 * dt * qw[0];
      quat[1] = quat[1] + 0.5 * dt * qw[1];
      quat[2] = quat[2] + 0.5 * dt * qw[2];
      quat[3] = quat[3] + 0.5 * dt * qw[3];

      // normalisation introduces the  stochastic drift term
      // to recover the Boltzmann distribution for the case of conservative torques
      MathExtra::qnormalize(quat);

      // next, update centre of mass positions and velocities

      // tmp now holds force in body frame
      MathExtra::matvec(rotationmatrix_transpose, f[i], tmp);
      // and then converts to gamma_t^{-1} * F (velocity) in body frame

      if (Tp_2D) {
        tmp[2] = 0.0;
        if (Tp_UNIFORM) {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0] + gamma_t_invsqrt[0] * (rng->uniform() - 0.5) * g2;
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1] + gamma_t_invsqrt[1] * (rng->uniform() - 0.5) * g2;
        } else if (Tp_GAUSS) {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0] + gamma_t_invsqrt[0] * rng->gaussian() * g2;
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1] + gamma_t_invsqrt[1] * rng->gaussian() * g2;
        } else {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0];
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1];
        }
      } else {
        if (Tp_UNIFORM) {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0] + gamma_t_invsqrt[0] * (rng->uniform() - 0.5) * g2;
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1] + gamma_t_invsqrt[1] * (rng->uniform() - 0.5) * g2;
          tmp[2] = g1 * tmp[2] * gamma_t_inv[2] + gamma_t_invsqrt[2] * (rng->uniform() - 0.5) * g2;
        } else if (Tp_GAUSS) {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0] + gamma_t_invsqrt[0] * rng->gaussian() * g2;
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1] + gamma_t_invsqrt[1] * rng->gaussian() * g2;
          tmp[2] = g1 * tmp[2] * gamma_t_inv[2] + gamma_t_invsqrt[2] * rng->gaussian() * g2;
        } else {
          tmp[0] = g1 * tmp[0] * gamma_t_inv[0];
          tmp[1] = g1 * tmp[1] * gamma_t_inv[1];
          tmp[2] = g1 * tmp[2] * gamma_t_inv[2];
        }
      }

      // finally, convert this back to lab-frame velocity and store in dv
      MathExtra::transpose_matvec(rotationmatrix_transpose, tmp, dv);

      v[i][0] = dv[0];
      v[i][1] = dv[1];
      v[i][2] = dv[2];

      x[i][0] += dv[0] * dt;
      x[i][1] += dv[1] * dt;
      x[i][2] += dv[2] * dt;

      if (Tp_DIPOLE) {
        MathExtra::quat_to_mat_trans(quat, rotationmatrix_transpose);
        MathExtra::transpose_matvec(rotationmatrix_transpose, dipole_body, f_rot);
        mu[i][0] = f_rot[0];
        mu[i][1] = f_rot[1];
        mu[i][2] = f_rot[2];
      }
    }
  }
}
