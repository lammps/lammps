/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include "fix_nve_sphere.h"
#include "atom.h"
#include "domain.h"
#include "atom_vec.h"
#include "force.h"
#include "error.h"
#include "math_vector.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

enum{NONE,DIPOLE};
enum{NODLM,DLM};

/* ---------------------------------------------------------------------- */

FixNVESphere::FixNVESphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/sphere command");

  time_integrate = 1;

  // process extra keywords
  // inertia = moment of inertia prefactor for sphere or disc

  extra = NONE;
  dlm = NODLM;
  inertia = 0.4;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/sphere command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else if (strcmp(arg[iarg+1],"dipole/dlm") == 0) {
        extra = DIPOLE;
        dlm = DLM;
      } else error->all(FLERR,"Illegal fix nve/sphere command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"disc")==0) {
      inertia = 0.5;
      if (domain->dimension != 2)
        error->all(FLERR,"Fix nve/sphere disc requires 2d simulation");
      iarg++;
    }
    else error->all(FLERR,"Illegal fix nve/sphere command");
  }

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nve/sphere requires atom style sphere");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/sphere update dipole requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::init()
{
  FixNVE::init();

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::initial_integrate(int /*vflag*/)
{
  double dtfm,dtirotate,msq,scale,s2,inv_len_mu;
  double g[3];
  vector w, w_temp, a;
  matrix Q, Q_temp, R;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    if (dlm == NODLM) {

      // d_mu/dt = omega cross mu
      // renormalize mu to dipole length

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

      // integrate orientation following Dullweber-Leimkuhler-Maclachlan scheme

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && mu[i][3] > 0.0) {

          // Construct Q from dipole:
          // Q is the rotation matrix from space frame to body frame
          // i.e. v_b = Q.v_s

          // define mu to lie along the z axis in the body frame
          // take the unit dipole to avoid getting a scaling matrix

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
          if (s2 != 0.0){ // i.e. the vectors are not parallel
            scale = (1.0 - a[2])/s2;

            Q[0][0] = 1.0 - scale*a[0]*a[0];
            Q[0][1] = -scale*a[0]*a[1];
            Q[0][2] = -a[0];
            Q[1][0] = -scale*a[0]*a[1];
            Q[1][1] = 1.0 - scale*a[1]*a[1];
            Q[1][2] = -a[1];
            Q[2][0] = a[0];
            Q[2][1] = a[1];
            Q[2][2] = 1.0 - scale*(a[0]*a[0] + a[1]*a[1]);
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
          omega[i][0] = w_temp[0];
          omega[i][1] = w_temp[1];
          omega[i][2] = w_temp[2];

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

void FixNVESphere::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  double rke = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
      rke += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
              omega[i][2]*omega[i][2])*radius[i]*radius[i]*rmass[i];
    }

}
