// clang-format off
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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "fix_nh_sphere.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

FixNHSphere::FixNHSphere(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!atom->omega_flag)
    error->all(FLERR,"Fix {} requires atom attribute omega", style);
  if (!atom->radius_flag)
    error->all(FLERR,"Fix {} requires atom attribute radius", style);

  // inertia = moment of inertia prefactor for sphere or disc

  inertia = 0.4;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"disc") == 0) {
      inertia = 0.5;
      if (domain->dimension != 2)
        error->all(FLERR, "Fix {} disc option requires 2d simulation", style);
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

void FixNHSphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nvt/npt/nph/sphere require extended particles");

  FixNH::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of rotational velocities
-----------------------------------------------------------------------*/

void FixNHSphere::nve_v()
{
  // standard nve_v velocity update

  FixNH::nve_v();

  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;
  double dtirotate;

  // update omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate*torque[i][0];
      omega[i][1] += dtirotate*torque[i][1];
      omega[i][2] += dtirotate*torque[i][2];
    }
}

/* ----------------------------------------------------------------------
   perform full-step update of position with dipole orientation, if requested
-----------------------------------------------------------------------*/

void FixNHSphere::nve_x()
{
  // standard nve_x position update

  FixNH::nve_x();

  // update mu for dipoles

  if (dipole_flag) {
    double **mu = atom->mu;
    double **omega = atom->omega;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (dlm_flag == 0) {
      // d_mu/dt = omega cross mu
      // renormalize mu to dipole length
      double msq,scale,g[3];

      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          if (mu[i][3] > 0.0) {
            g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
            g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
            g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
            msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
            scale = mu[i][3]/sqrt(msq);
            mu[i][0] = g[0]*scale;
            mu[i][1] = g[1]*scale;
            mu[i][2] = g[2]*scale;
          }
    } else {
      // Integrate orientation following Dullweber-Leimkuhler-Maclachlan scheme
      double w[3], w_temp[3], a[3];
      double Q[3][3], Q_temp[3][3], R[3][3];
      double scale,s2,inv_len_mu;

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && mu[i][3] > 0.0) {

          // Construct Q from dipole:
          // Q is the rotation matrix from space frame to body frame
          // i.e. v_b = Q.v_s

          // Define mu to lie along the z axis in the body frame
          // We take the unit dipole to avoid getting a scaling matrix
          inv_len_mu = 1.0/mu[i][3];
          a[0] = mu[i][0]*inv_len_mu;
          a[1] = mu[i][1]*inv_len_mu;
          a[2] = mu[i][2]*inv_len_mu;

          // v = a x [0 0 1] - cross product of mu in space and body frames
          // s = |v|
          // c = a.[0 0 1] = a[2]
          // vx = [ 0    -v[2]  v[1]
          //        v[2]  0    -v[0]
          //       -v[1]  v[0]  0    ]
          // then
          // Q = I + vx + vx^2 * (1-c)/s^2

          s2 = a[0]*a[0] + a[1]*a[1];
          if (s2 != 0.0) { // i.e. the vectors are not parallel
            scale = (1.0 - a[2])/s2;

            Q[0][0] = 1.0 - scale*a[0]*a[0]; Q[0][1] = -scale*a[0]*a[1];      Q[0][2] = -a[0];
            Q[1][0] = -scale*a[0]*a[1];      Q[1][1] = 1.0 - scale*a[1]*a[1]; Q[1][2] = -a[1];
            Q[2][0] = a[0];                  Q[2][1] = a[1];                  Q[2][2] = 1.0 - scale*(a[0]*a[0] + a[1]*a[1]);
          } else { // if parallel then we just have I or -I
            Q[0][0] = 1.0/a[2];  Q[0][1] = 0.0;       Q[0][2] = 0.0;
            Q[1][0] = 0.0;       Q[1][1] = 1.0/a[2];  Q[1][2] = 0.0;
            Q[2][0] = 0.0;       Q[2][1] = 0.0;       Q[2][2] = 1.0/a[2];
          }

          // Local copy of this particle's angular velocity (in space frame)
          w[0] = omega[i][0]; w[1] = omega[i][1]; w[2] = omega[i][2];

          // Transform omega into body frame: w_temp= Q.w
          matvec(Q,w,w_temp);

          // Construct rotation R1
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);

          // Apply R1 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R1 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Construct rotation R2
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);

          // Apply R2 to w: w_temp = R.w
          matvec(R,w,w_temp);

          // Apply R2 to Q: Q = R^T.Q_temp
          transpose_times3(R,Q_temp,Q);

          // Construct rotation R3
          BuildRzMatrix(R, 2.0*dtf/force->ftm2v*w_temp[2]);

          // Apply R3 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R3 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Construct rotation R4
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);

          // Apply R4 to w: w_temp = R.w
          matvec(R,w,w_temp);

          // Apply R4 to Q: Q = R^T.Q_temp
          transpose_times3(R,Q_temp,Q);

          // Construct rotation R5
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);

          // Apply R5 to w: w = R.w_temp
          matvec(R,w_temp,w);

          // Apply R5 to Q: Q_temp = R^T.Q
          transpose_times3(R,Q,Q_temp);

          // Transform w back into space frame w_temp = Q^T.w
          transpose_matvec(Q_temp,w,w_temp);
          omega[i][0] = w_temp[0]; omega[i][1] = w_temp[1]; omega[i][2] = w_temp[2];

          // Set dipole according to updated Q: mu = Q^T.[0 0 1] * |mu|
          mu[i][0] = Q_temp[2][0] * mu[i][3];
          mu[i][1] = Q_temp[2][1] * mu[i][3];
          mu[i][2] = Q_temp[2][2] * mu[i][3];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of rotatonal velocities
-----------------------------------------------------------------------*/

void FixNHSphere::nh_v_temp()
{
  // standard nh_v_temp scaling

  FixNH::nh_v_temp();

  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      omega[i][0] *= factor_eta;
      omega[i][1] *= factor_eta;
      omega[i][2] *= factor_eta;
    }
  }
}
