// clang-format off
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

#include "omp_compat.h"
#include "fix_nve_sphere_omp.h"

#include "atom.h"
#include "force.h"
#include "math_extra.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};
enum{NODLM,DLM};

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::initial_integrate(int /* vflag */)
{
  double * const * const x = atom->x;
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  double * const * const omega = atom->omega;
  const double * const * const torque = atom->torque;
  const double * const radius = atom->radius;
  const double * const rmass = atom->rmass;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  // set timestep here since dt may have changed or come via rRESPA
  const double dtfrotate = dtf / INERTIA;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE
#endif
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double * const * const mu = atom->mu;
    if (dlm == NODLM) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE
#endif
      for (int i = 0; i < nlocal; i++) {
        double g0,g1,g2,msq,scale;
        if (mask[i] & groupbit) {
          if (mu[i][3] > 0.0) {
            g0 = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
            g1 = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
            g2 = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
            msq = g0*g0 + g1*g1 + g2*g2;
            scale = mu[i][3]/sqrt(msq);
            mu[i][0] = g0*scale;
            mu[i][1] = g1*scale;
            mu[i][2] = g2*scale;
          }
        }
      }
    } else {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE
#endif
      // Integrate orientation following Dullweber-Leimkuhler-Maclachlan scheme
      for (int i = 0; i < nlocal; i++) {
        double w[3], w_temp[3], a[3];
        double Q[3][3], Q_temp[3][3], R[3][3];

        if (mask[i] & groupbit && mu[i][3] > 0.0) {

          // Construct Q from dipole:
          // Q is the rotation matrix from space frame to body frame
          // i.e. v_b = Q.v_s

          // Define mu to lie along the z axis in the body frame
          // We take the unit dipole to avoid getting a scaling matrix
          const double inv_len_mu = 1.0/mu[i][3];
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

          const double s2 = a[0]*a[0] + a[1]*a[1];
          if (s2 != 0.0) { // i.e. the vectors are not parallel
            const double scale = (1.0 - a[2])/s2;

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

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::final_integrate()
{
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  double * const * const omega = atom->omega;
  const double * const * const torque = atom->torque;
  const double * const rmass = atom->rmass;
  const double * const radius = atom->radius;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  // set timestep here since dt may have changed or come via rRESPA

  const double dtfrotate = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE
#endif
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}
