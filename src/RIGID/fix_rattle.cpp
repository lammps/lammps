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

/* ----------------------------------------------------------------------
   Contributing author: Peter Wirnsberger (University of Cambridge)
------------------------------------------------------------------------- */

#include "fix_rattle.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

// set RATTLE_DEBUG  1 to check constraints at end of timestep

#define RATTLE_DEBUG 0

// set RATTLE_RAISE_ERROR 1 if you want this fix to raise
//     an error if the constraints cannot be satisfied

#define RATTLE_RAISE_ERROR 0

// You can enable velocity and coordinate checks separately
// by setting RATTLE_TEST_VEL/POS true

#define RATTLE_TEST_VEL false
#define RATTLE_TEST_POS false

enum{V,VP,XSHAKE};

/* ---------------------------------------------------------------------- */

FixRattle::FixRattle(LAMMPS *lmp, int narg, char **arg) :
  FixShake(lmp, narg, arg)
{
  rattle = 1;

  // allocate memory for unconstrained velocity update

  vp = nullptr;
  FixRattle::grow_arrays(atom->nmax);

  // default communication mode
  // necessary for compatibility with SHAKE
  // see pack_forward and unpack_forward

  comm_mode = XSHAKE;
  vflag_post_force = 0;

  verr_max = 0;
  derr_max = 0;
}

/* ---------------------------------------------------------------------- */

FixRattle::~FixRattle()
{
  memory->destroy(vp);

  if (RATTLE_DEBUG) {

    // communicate maximum distance error

    double global_derr_max, global_verr_max;
    int npid;

    MPI_Reduce(&derr_max, &global_derr_max, 1 , MPI_DOUBLE, MPI_MAX, 0, world);
    MPI_Reduce(&verr_max, &global_verr_max, 1 , MPI_DOUBLE, MPI_MAX, 0, world);

    MPI_Comm_rank (world, &npid); // Find out process rank

    if (npid == 0 && screen) {
      fprintf(screen, "RATTLE: Maximum overall relative position error ( (r_ij-d_ij)/d_ij ): %.10g\n", global_derr_max);
      fprintf(screen, "RATTLE: Maximum overall absolute velocity error (r_ij * v_ij): %.10g\n", global_verr_max);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixRattle::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  if (RATTLE_DEBUG) mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   initialize RATTLE and check that this is the last final_integrate fix
------------------------------------------------------------------------- */

void FixRattle::init() {

  // initialize SHAKE first

  FixShake::init();

  // show a warning if any final-integrate fix comes after this one

  int after = 0;
  int flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(id,modify->fix[i]->id) == 0) after = 1;
    else if ((modify->fmask[i] & FINAL_INTEGRATE) && after) flag = 1;
  }

  if (flag && comm->me == 0)
    error->warning(FLERR,
                   "Fix rattle should come after all other integration fixes ");
}

/* ----------------------------------------------------------------------
   This method carries out an unconstrained velocity update first and
   then applies the velocity corrections directly (v and vp are modified).
------------------------------------------------------------------------- */

void FixRattle::post_force(int vflag)
{
  // remember vflag for the coordinate correction in this->final_integrate

  vflag_post_force = vflag;

  // unconstrained velocity update by half a timestep
  // similar to FixShake::unconstrained_update()

  update_v_half_nocons();

  // communicate the unconstrained velocities

  if (nprocs > 1) {
    comm_mode = VP;
    comm->forward_comm_fix(this);
  }

  // correct the velocity for each molecule accordingly

  int m;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    if      (shake_flag[m] == 2)        vrattle2(m);
    else if (shake_flag[m] == 3)        vrattle3(m);
    else if (shake_flag[m] == 4)        vrattle4(m);
    else                                vrattle3angle(m);
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // remember vflag for the coordinate correction in this->final_integrate

  vflag_post_force = vflag;

  // unconstrained velocity update by half a timestep
  // similar to FixShake::unconstrained_update()

  update_v_half_nocons_respa(ilevel);

  // communicate the unconstrained velocities

  if (nprocs > 1) {
    comm_mode = VP;
    comm->forward_comm_fix(this);
  }

  // correct the velocity for each molecule accordingly

  int m;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    if      (shake_flag[m] == 2)        vrattle2(m);
    else if (shake_flag[m] == 3)        vrattle3(m);
    else if (shake_flag[m] == 4)        vrattle4(m);
    else                                vrattle3angle(m);
  }
}

/* ----------------------------------------------------------------------
   let SHAKE calculate the constraining forces for the coordinates
------------------------------------------------------------------------- */

void FixRattle::final_integrate()
{
  comm_mode = XSHAKE;
  FixShake::post_force(vflag_post_force);
}

/* ---------------------------------------------------------------------- */

void FixRattle::final_integrate_respa(int ilevel, int iloop)
{
  comm_mode = XSHAKE;
  FixShake::post_force_respa(vflag_post_force, ilevel, iloop);
}

/* ----------------------------------------------------------------------
   correct velocities of molecule m with 2 constraints bonds and 1 angle
------------------------------------------------------------------------- */

void FixRattle::vrattle3angle(int m)
{
  tagint i0,i1,i2;
  double c[3], l[3], a[3][3], r01[3], imass[3],
         r02[3], r12[3], vp01[3], vp02[3], vp12[3];

  // local atom IDs and constraint distances

  i0 = atom->map(shake_atom[m][0]);
  i1 = atom->map(shake_atom[m][1]);
  i2 = atom->map(shake_atom[m][2]);

  // r01,r02,r12 = distance vec between atoms

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);
  MathExtra::sub3(x[i2],x[i1],r12);

  // take into account periodicity

  domain->minimum_image(r01);
  domain->minimum_image(r02);
  domain->minimum_image(r12);

  // v01,v02,v12 = velocity differences

  MathExtra::sub3(vp[i1],vp[i0],vp01);
  MathExtra::sub3(vp[i2],vp[i0],vp02);
  MathExtra::sub3(vp[i2],vp[i1],vp12);

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    imass[0] = 1.0/rmass[i0];
    imass[1] = 1.0/rmass[i1];
    imass[2] = 1.0/rmass[i2];
  } else {
    imass[0] = 1.0/mass[type[i0]];
    imass[1] = 1.0/mass[type[i1]];
    imass[2] = 1.0/mass[type[i2]];
  }

  // setup matrix

  a[0][0]   =   (imass[1] + imass[0])   * MathExtra::dot3(r01,r01);
  a[0][1]   =   (imass[0]           )   * MathExtra::dot3(r01,r02);
  a[0][2]   =   (-imass[1]          )   * MathExtra::dot3(r01,r12);
  a[1][0]   =   a[0][1];
  a[1][1]   =   (imass[0] + imass[2])   * MathExtra::dot3(r02,r02);
  a[1][2]   =   (imass[2]           )   * MathExtra::dot3(r02,r12);
  a[2][0]   =   a[0][2];
  a[2][1]   =   a[1][2];
  a[2][2]   =   (imass[2] + imass[1])   * MathExtra::dot3(r12,r12);

  // sestup RHS

  c[0]  = -MathExtra::dot3(vp01,r01);
  c[1]  = -MathExtra::dot3(vp02,r02);
  c[2]  = -MathExtra::dot3(vp12,r12);

  // calculate the inverse matrix exactly

  solve3x3exactly(a,c,l);

  // add corrections to the velocities if processor owns atom

  if (i0 < nlocal) {
    for (int k=0; k<3; k++)
      v[i0][k]  -=  imass[0]*  (  l[0] * r01[k] + l[1] * r02[k] );
  }
  if (i1 < nlocal) {
    for (int k=0; k<3; k++)
      v[i1][k]  -=  imass[1] * ( -l[0] * r01[k] + l[2] * r12[k] );
  }
  if (i2 < nlocal) {
    for (int k=0; k<3; k++)
      v[i2][k] -=   imass[2] * ( -l[1] * r02[k] - l[2] * r12[k] );
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::vrattle2(int m)
{
  tagint    i0, i1;
  double    imass[2], r01[3], vp01[3];

  // local atom IDs and constraint distances

  i0 = atom->map(shake_atom[m][0]);
  i1 = atom->map(shake_atom[m][1]);

  // r01 = distance vec between atoms, with PBC

  MathExtra::sub3(x[i1],x[i0],r01);
  domain->minimum_image(r01);

  // v01 = distance vectors for velocities

  MathExtra::sub3(vp[i1],vp[i0],vp01);

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    imass[0] = 1.0/rmass[i0];
    imass[1] = 1.0/rmass[i1];
  } else {
    imass[0] = 1.0/mass[type[i0]];
    imass[1] = 1.0/mass[type[i1]];
  }

  // Lagrange multiplier: exact solution

  double l01 = - MathExtra::dot3(r01,vp01) /
    (MathExtra::dot3(r01,r01) * (imass[0] + imass[1]));

  // add corrections to the velocities if the process owns this atom

  if (i0 < nlocal) {
    for (int k=0; k<3; k++)
      v[i0][k] -= imass[0] * l01 * r01[k];
  }
  if (i1 < nlocal) {
    for (int k=0; k<3; k++)
      v[i1][k] -= imass[1] * (-l01) * r01[k];
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::vrattle3(int m)
{
  tagint    i0,i1,i2;
  double    imass[3], r01[3], r02[3], vp01[3], vp02[3],
            a[2][2],c[2],l[2];

  // local atom IDs and constraint distances

  i0 = atom->map(shake_atom[m][0]);
  i1 = atom->map(shake_atom[m][1]);
  i2 = atom->map(shake_atom[m][2]);

  // r01,r02 = distance vec between atoms, with PBC

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);

  domain->minimum_image(r01);
  domain->minimum_image(r02);

  // vp01,vp02 =  distance vectors between velocities

  MathExtra::sub3(vp[i1],vp[i0],vp01);
  MathExtra::sub3(vp[i2],vp[i0],vp02);

  if (rmass) {
    imass[0] = 1.0/rmass[i0];
    imass[1] = 1.0/rmass[i1];
    imass[2] = 1.0/rmass[i2];
  } else {
    imass[0] = 1.0/mass[type[i0]];
    imass[1] = 1.0/mass[type[i1]];
    imass[2] = 1.0/mass[type[i2]];
  }

  // setup matrix

  a[0][0]   =   (imass[1] + imass[0])   * MathExtra::dot3(r01,r01);
  a[0][1]   =   (imass[0]           )   * MathExtra::dot3(r01,r02);
  a[1][0]   =   a[0][1];
  a[1][1]   =   (imass[0] + imass[2])   * MathExtra::dot3(r02,r02);

  // setup RHS

  c[0]  = - MathExtra::dot3(vp01,r01);
  c[1]  = - MathExtra::dot3(vp02,r02);

  // calculate the inverse 2x2 matrix exactly

  solve2x2exactly(a,c,l);

  // add corrections to the velocities if the process owns this atom

  if (i0 < nlocal) {
    for (int k=0; k<3; k++)
      v[i0][k] -= imass[0] * (  l[0] * r01[k] + l[1] * r02[k] );
  }
  if (i1 < nlocal)
    for (int k=0; k<3; k++) {
      v[i1][k] -= imass[1] * ( -l[0] * r01[k] );
  }
  if (i2 < nlocal) {
    for (int k=0; k<3; k++)
      v[i2][k] -= imass[2] * ( -l[1] * r02[k] );
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::vrattle4(int m)
{
  tagint    i0,i1,i2,i3;
  double    imass[4], c[3], l[3], a[3][3],
            r01[3], r02[3], r03[3], vp01[3], vp02[3], vp03[3];

  // local atom IDs and constraint distances

  i0 = atom->map(shake_atom[m][0]);
  i1 = atom->map(shake_atom[m][1]);
  i2 = atom->map(shake_atom[m][2]);
  i3 = atom->map(shake_atom[m][3]);

  // r01,r02,r12 = distance vec between atoms, with PBC

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);
  MathExtra::sub3(x[i3],x[i0],r03);

  domain->minimum_image(r01);
  domain->minimum_image(r02);
  domain->minimum_image(r03);

  // vp01,vp02,vp03 = distance vectors between velocities

  MathExtra::sub3(vp[i1],vp[i0],vp01);
  MathExtra::sub3(vp[i2],vp[i0],vp02);
  MathExtra::sub3(vp[i3],vp[i0],vp03);

  // matrix coeffs and rhs for lamda equations

  if (rmass) {
    imass[0] = 1.0/rmass[i0];
    imass[1] = 1.0/rmass[i1];
    imass[2] = 1.0/rmass[i2];
    imass[3] = 1.0/rmass[i3];
  } else {
    imass[0] = 1.0/mass[type[i0]];
    imass[1] = 1.0/mass[type[i1]];
    imass[2] = 1.0/mass[type[i2]];
    imass[3] = 1.0/mass[type[i3]];
  }

  // setup matrix

  a[0][0]   =   (imass[0] + imass[1])   * MathExtra::dot3(r01,r01);
  a[0][1]   =   (imass[0]           )   * MathExtra::dot3(r01,r02);
  a[0][2]   =   (imass[0]           )   * MathExtra::dot3(r01,r03);
  a[1][0]   =   a[0][1];
  a[1][1]   =   (imass[0] + imass[2])   * MathExtra::dot3(r02,r02);
  a[1][2]   =   (imass[0]           )   * MathExtra::dot3(r02,r03);
  a[2][0]   =   a[0][2];
  a[2][1]   =   a[1][2];
  a[2][2]   =   (imass[0] + imass[3])   * MathExtra::dot3(r03,r03);

  // setup RHS

  c[0]  = - MathExtra::dot3(vp01,r01);
  c[1]  = - MathExtra::dot3(vp02,r02);
  c[2]  = - MathExtra::dot3(vp03,r03);

  // calculate the inverse 3x3 matrix exactly

  solve3x3exactly(a,c,l);

  // add corrections to the velocities if the process owns this atom

  if (i0 < nlocal) {
    for (int k=0; k<3; k++)
      v[i0][k] -= imass[0] * (  l[0] * r01[k] + l[1] * r02[k] + l[2] * r03[k]);
  }
  if (i1 < nlocal) {
    for (int k=0; k<3; k++)
      v[i1][k] -= imass[1] * (-l[0] * r01[k]);
  }
  if (i2 < nlocal) {
    for (int k=0; k<3; k++)
      v[i2][k] -= imass[2] * ( -l[1] * r02[k]);
  }
  if (i3 < nlocal) {
    for (int k=0; k<3; k++)
      v[i3][k] -= imass[3] * ( -l[2] * r03[k]);
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::solve2x2exactly(const double a[][2],
                                const double c[], double l[])
{
  double determ, determinv;

  // calculate the determinant of the matrix

  determ = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  // check if matrix is actually invertible

  if (determ == 0.0) error->one(FLERR,"Rattle determinant = 0.0");
  determinv = 1.0/determ;

  // Calculate the solution:  (l01, l02)^T = A^(-1) * c

  l[0] = determinv * ( a[1][1] * c[0]  - a[0][1] * c[1]);
  l[1] = determinv * (-a[1][0] * c[0]  + a[0][0] * c[1]);
}

/* ---------------------------------------------------------------------- */

void FixRattle::solve3x3exactly(const double a[][3],
                                const double c[], double l[])
{
  double ai[3][3];
  double determ, determinv;

  // calculate the determinant of the matrix

  determ = a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] +
    a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] -
    a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0];

  // check if matrix is actually invertible

  if (determ == 0.0) error->one(FLERR,"Rattle determinant = 0.0");

  // calculate the inverse 3x3 matrix: A^(-1) = (ai_jk)

  determinv = 1.0/determ;
  ai[0][0] =  determinv * (a[1][1]*a[2][2] - a[1][2]*a[2][1]);
  ai[0][1] = -determinv * (a[0][1]*a[2][2] - a[0][2]*a[2][1]);
  ai[0][2] =  determinv * (a[0][1]*a[1][2] - a[0][2]*a[1][1]);
  ai[1][0] = -determinv * (a[1][0]*a[2][2] - a[1][2]*a[2][0]);
  ai[1][1] =  determinv * (a[0][0]*a[2][2] - a[0][2]*a[2][0]);
  ai[1][2] = -determinv * (a[0][0]*a[1][2] - a[0][2]*a[1][0]);
  ai[2][0] =  determinv * (a[1][0]*a[2][1] - a[1][1]*a[2][0]);
  ai[2][1] = -determinv * (a[0][0]*a[2][1] - a[0][1]*a[2][0]);
  ai[2][2] =  determinv * (a[0][0]*a[1][1] - a[0][1]*a[1][0]);

  // calculate the solution:  (l01, l02, l12)^T = A^(-1) * c

  for (int i=0; i<3; i++) {
    l[i] = 0;
    for (int j=0; j<3; j++)
      l[i] += ai[i][j] * c[j];
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::reset_dt()
{
  FixShake::reset_dt();
}

/* ----------------------------------------------------------------------
   carry out an unconstrained velocity update (vp is modified)
------------------------------------------------------------------------- */

void FixRattle::update_v_half_nocons()
{
  const double dtfv = 0.5 * update->dt * force->ftm2v;
  double dtfvinvm;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (shake_flag[i]) {
        dtfvinvm = dtfv / rmass[i];
        for (int k=0; k<3; k++)
          vp[i][k] = v[i][k] + dtfvinvm * f[i][k];
      }
      else
        vp[i][0] = vp[i][1] = vp[i][2] = 0;
    }
  }
  else {
    for (int i = 0; i < nlocal; i++) {
      dtfvinvm = dtfv/mass[type[i]];
      if (shake_flag[i]) {
        for (int k=0; k<3; k++)
          vp[i][k] = v[i][k] + dtfvinvm * f[i][k];
      }
      else
        vp[i][0] = vp[i][1] = vp[i][2] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRattle::update_v_half_nocons_respa(int /*ilevel*/)
{
  // carry out unconstrained velocity update

  update_v_half_nocons();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixRattle::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = FixShake::memory_usage();
  bytes += (double)nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixRattle::grow_arrays(int nmax)
{
  FixShake::grow_arrays(nmax);
  memory->destroy(vp);
  memory->create(vp,nmax,3,"rattle:vp");
}

/* ---------------------------------------------------------------------- */

int FixRattle::pack_forward_comm(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  m = 0;

  switch (comm_mode) {
    case XSHAKE:
      m = FixShake::pack_forward_comm(n, list, buf, pbc_flag, pbc);
      break;
    case VP:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = vp[j][0];
        buf[m++] = vp[j][1];
        buf[m++] = vp[j][2];
      }
      break;

    case V:
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
      break;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRattle::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;

  switch (comm_mode) {
    case XSHAKE:
      FixShake::unpack_forward_comm(n, first,buf);
      break;

    case VP:
      for (i = first; i < last; i++) {
        vp[i][0] = buf[m++];
        vp[i][1] = buf[m++];
        vp[i][2] = buf[m++];
      }
      break;

    case V:
      for (i = first; i < last; i++) {
        v[i][0] = buf[m++];
        v[i][1] = buf[m++];
        v[i][2] = buf[m++];
      }
      break;
  }
}


/* ----------------------------------------------------------------------
  Let shake calculate new constraining forces for the coordinates;
  As opposed to the regular shake call, this method is usually called from
  end_of_step fixes after the second velocity integration has happened.
------------------------------------------------------------------------- */

void FixRattle::shake_end_of_step(int vflag) {

  if (nprocs > 1) {
    comm_mode = V;
    comm->forward_comm_fix(this);
  }

  comm_mode = XSHAKE;
  FixShake::shake_end_of_step(vflag);
}


/* ----------------------------------------------------------------------
  Let shake calculate new constraining forces and correct the
  coordinates. Nothing to do for rattle here.
------------------------------------------------------------------------- */

void FixRattle::correct_coordinates(int vflag) {
  comm_mode = XSHAKE;
  FixShake::correct_coordinates(vflag);
}


/* ----------------------------------------------------------------------
   Remove the velocity component along any bond.
------------------------------------------------------------------------- */

void FixRattle::correct_velocities() {

  // Copy current velocities instead of unconstrained_update, because the correction
  // should happen instantaneously and not after the next half step.

  for (int i = 0; i < atom->nlocal; i++) {
    if (shake_flag[i]) {
      for (int k=0; k<3; k++)
        vp[i][k] = v[i][k];
    }
    else
      vp[i][0] = vp[i][1] = vp[i][2] = 0;
  }

  // communicate the unconstrained velocities

  if (nprocs > 1) {
    comm_mode = VP;
    comm->forward_comm_fix(this);
  }

  // correct the velocity for each molecule accordingly

  int m;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    if      (shake_flag[m] == 2)        vrattle2(m);
    else if (shake_flag[m] == 3)        vrattle3(m);
    else if (shake_flag[m] == 4)        vrattle4(m);
    else                                vrattle3angle(m);
  }
}


/* ----------------------------------------------------------------------
   DEBUGGING methods
   The functions below allow you to check whether the
     coordinate and velocity constraints are satisfied at the
     end of the timestep
   only enabled if RATTLE_DEBUG is set to 1 at top of file
   checkX tests if shakeX and vrattleX worked as expected
------------------------------------------------------------------------- */

void FixRattle::end_of_step()
{
  if (nprocs > 1) {
       comm_mode = V;
       comm->forward_comm_fix(this);
  }
  if (!check_constraints(v, RATTLE_TEST_POS, RATTLE_TEST_VEL) && RATTLE_RAISE_ERROR) {
    error->one(FLERR, "Rattle failed ");
  }
}


/* ---------------------------------------------------------------------- */

bool FixRattle::check_constraints(double **v, bool checkr, bool checkv)
{
  int m;
  bool ret = true;
  int i=0;
  while (i < nlist && ret) {
    m = list[i];
    if      (shake_flag[m] == 2)     ret =   check2(v,m,checkr,checkv);
    else if (shake_flag[m] == 3)     ret =   check3(v,m,checkr,checkv);
    else if (shake_flag[m] == 4)     ret =   check4(v,m,checkr,checkv);
    else                             ret =   check3angle(v,m,checkr,checkv);
    i++;
    if (!RATTLE_RAISE_ERROR)         ret = true;
  }
  return ret;
}

/* ---------------------------------------------------------------------- */

bool FixRattle::check2(double **v, int m, bool checkr, bool checkv)
{
  bool      stat;
  double    r01[3],v01[3];
  const double tol = tolerance;
  double bond1 = bond_distance[shake_type[m][0]];

  tagint i0 = atom->map(shake_atom[m][0]);
  tagint i1 = atom->map(shake_atom[m][1]);

  MathExtra::sub3(x[i1],x[i0],r01);
  domain->minimum_image(r01);
  MathExtra::sub3(v[i1],v[i0],v01);

  stat = !(checkr && (fabs(sqrt(MathExtra::dot3(r01,r01)) - bond1) > tol));
  if (!stat)
     error->one(FLERR,"Coordinate constraints are not satisfied "
                "up to desired tolerance ");

  stat = !(checkv && (fabs(MathExtra::dot3(r01,v01)) > tol));
  if (!stat)
     error->one(FLERR,"Velocity constraints are not satisfied "
                "up to desired tolerance ");
  return stat;
}

/* ---------------------------------------------------------------------- */

bool FixRattle::check3(double **v, int m, bool checkr, bool checkv)
{
  bool      stat;
  tagint    i0,i1,i2;
  double    r01[3], r02[3], v01[3], v02[3];
  const double tol = tolerance;
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];

  i0 = atom->map(shake_atom[m][0]);
  i1 = atom->map(shake_atom[m][1]);
  i2 = atom->map(shake_atom[m][2]);

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);

  domain->minimum_image(r01);
  domain->minimum_image(r02);

  MathExtra::sub3(v[i1],v[i0],v01);
  MathExtra::sub3(v[i2],v[i0],v02);

  stat = !(checkr && (fabs(sqrt(MathExtra::dot3(r01,r01)) - bond1) > tol ||
                      fabs(sqrt(MathExtra::dot3(r02,r02))-bond2) > tol));
  if (!stat)
     error->one(FLERR,"Coordinate constraints are not satisfied "
                "up to desired tolerance ");

  stat = !(checkv && (fabs(MathExtra::dot3(r01,v01)) > tol ||
                      fabs(MathExtra::dot3(r02,v02)) > tol));
  if (!stat)
     error->one(FLERR,"Velocity constraints are not satisfied "
                "up to desired tolerance ");
  return stat;
}

/* ---------------------------------------------------------------------- */

bool FixRattle::check4(double **v, int m, bool checkr, bool checkv)
{
  bool stat = true;
  const double tol = tolerance;
  double r01[3], r02[3], r03[3], v01[3], v02[3], v03[3];

  int i0 = atom->map(shake_atom[m][0]);
  int i1 = atom->map(shake_atom[m][1]);
  int i2 = atom->map(shake_atom[m][2]);
  int i3 = atom->map(shake_atom[m][3]);
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];
  double bond3 = bond_distance[shake_type[m][2]];

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);
  MathExtra::sub3(x[i3],x[i0],r03);

  domain->minimum_image(r01);
  domain->minimum_image(r02);
  domain->minimum_image(r03);

  MathExtra::sub3(v[i1],v[i0],v01);
  MathExtra::sub3(v[i2],v[i0],v02);
  MathExtra::sub3(v[i3],v[i0],v03);

  stat = !(checkr && (fabs(sqrt(MathExtra::dot3(r01,r01)) - bond1) > tol ||
                      fabs(sqrt(MathExtra::dot3(r02,r02))-bond2) > tol ||
                      fabs(sqrt(MathExtra::dot3(r03,r03))-bond3) > tol));
  if (!stat)
     error->one(FLERR,"Coordinate constraints are not satisfied "
                "up to desired tolerance ");

  stat = !(checkv && (fabs(MathExtra::dot3(r01,v01)) > tol ||
                      fabs(MathExtra::dot3(r02,v02)) > tol ||
                      fabs(MathExtra::dot3(r03,v03)) > tol));
  if (!stat)
     error->one(FLERR,"Velocity constraints are not satisfied "
                "up to desired tolerance ");
  return stat;
}

/* ---------------------------------------------------------------------- */

bool FixRattle::check3angle(double **v, int m, bool checkr, bool checkv)
{
  bool stat = true;
  const double tol = tolerance;
  double r01[3], r02[3], r12[3], v01[3], v02[3], v12[3];

  int i0 = atom->map(shake_atom[m][0]);
  int i1 = atom->map(shake_atom[m][1]);
  int i2 = atom->map(shake_atom[m][2]);
  double bond1 = bond_distance[shake_type[m][0]];
  double bond2 = bond_distance[shake_type[m][1]];
  double bond12 = angle_distance[shake_type[m][2]];

  MathExtra::sub3(x[i1],x[i0],r01);
  MathExtra::sub3(x[i2],x[i0],r02);
  MathExtra::sub3(x[i2],x[i1],r12);

  domain->minimum_image(r01);
  domain->minimum_image(r02);
  domain->minimum_image(r12);

  MathExtra::sub3(v[i1],v[i0],v01);
  MathExtra::sub3(v[i2],v[i0],v02);
  MathExtra::sub3(v[i2],v[i1],v12);



  double db1 = fabs(sqrt(MathExtra::dot3(r01,r01)) - bond1);
  double db2 = fabs(sqrt(MathExtra::dot3(r02,r02))-bond2);
  double db12 = fabs(sqrt(MathExtra::dot3(r12,r12))-bond12);


  stat = !(checkr && (db1 > tol ||
                      db2 > tol ||
                      db12 > tol));

  if (derr_max < db1/bond1)    derr_max = db1/bond1;
  if (derr_max < db2/bond2)    derr_max = db2/bond2;
  if (derr_max < db12/bond12)  derr_max = db12/bond12;


  if (!stat && RATTLE_RAISE_ERROR)
     error->one(FLERR,"Coordinate constraints are not satisfied "
                "up to desired tolerance ");

  double dv1 = fabs(MathExtra::dot3(r01,v01));
  double dv2 = fabs(MathExtra::dot3(r02,v02));
  double dv12 = fabs(MathExtra::dot3(r12,v12));

  if (verr_max < dv1)    verr_max = dv1;
  if (verr_max < dv2)    verr_max = dv2;
  if (verr_max < dv12)   verr_max = dv12;


  stat = !(checkv && (dv1 > tol ||
                      dv2 > tol ||
                      dv12> tol));


  if (!stat && RATTLE_RAISE_ERROR)
     error->one(FLERR,"Velocity constraints are not satisfied "
                "up to desired tolerance!");


  return stat;
}
