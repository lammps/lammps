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

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_nh_omp.h"
#include <cmath>
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{ISO,ANISO,TRICLINIC};

#define TILTMAX 1.5

typedef struct { double x,y,z; } dbl3_t;

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or dilate group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixNHOMP::remap()
{
  double oldlo,oldhi,expfac;

  double * const * _noalias const x = atom->x;
  const int * _noalias const mask = atom->mask;
  const int nlocal = atom->nlocal;
  double * _noalias const h = domain->h;

  // omega is not used, except for book-keeping

  for (int i = 0; i < 6; i++) omega[i] += dto*omega_dot[i];

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap) domain->x2lamda(nlocal);
  else {
    int i;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->x2lamda(x[i],x[i]);
  }

  if (nrigid)
    for (int i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  // this operation corresponds to applying the
  // translate and scale operations
  // corresponding to the solution of the following ODE:
  //
  // h_dot = omega_dot * h
  //
  // where h_dot, omega_dot and h are all upper-triangular
  // 3x3 tensors. In Voigt notation, the elements of the
  // RHS product tensor are:
  // h_dot = [0*0, 1*1, 2*2, 1*3+3*2, 0*4+5*3+4*2, 0*5+5*1]
  //
  // Ordering of operations preserves time symmetry.

  const double dto2 = dto/2.0;
  const double dto4 = dto/4.0;
  const double dto8 = dto/8.0;

  // off-diagonal components, first half

  if (pstyle == TRICLINIC) {

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

    if (p_flag[3]) {
      expfac = exp(dto4*omega_dot[1]);
      h[3] *= expfac;
      h[3] += dto2*(omega_dot[3]*h[2]);
      h[3] *= expfac;
    }

    if (p_flag[5]) {
      expfac = exp(dto4*omega_dot[0]);
      h[5] *= expfac;
      h[5] += dto2*(omega_dot[5]*h[1]);
      h[5] *= expfac;
    }

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }
  }

  // scale diagonal components
  // scale tilt factors with cell, if set

  if (p_flag[0]) {
    oldlo = domain->boxlo[0];
    oldhi = domain->boxhi[0];
    expfac = exp(dto*omega_dot[0]);
    domain->boxlo[0] = (oldlo-fixedpoint[0])*expfac + fixedpoint[0];
    domain->boxhi[0] = (oldhi-fixedpoint[0])*expfac + fixedpoint[0];
  }

  if (p_flag[1]) {
    oldlo = domain->boxlo[1];
    oldhi = domain->boxhi[1];
    expfac = exp(dto*omega_dot[1]);
    domain->boxlo[1] = (oldlo-fixedpoint[1])*expfac + fixedpoint[1];
    domain->boxhi[1] = (oldhi-fixedpoint[1])*expfac + fixedpoint[1];
    if (scalexy) h[5] *= expfac;
  }

  if (p_flag[2]) {
    oldlo = domain->boxlo[2];
    oldhi = domain->boxhi[2];
    expfac = exp(dto*omega_dot[2]);
    domain->boxlo[2] = (oldlo-fixedpoint[2])*expfac + fixedpoint[2];
    domain->boxhi[2] = (oldhi-fixedpoint[2])*expfac + fixedpoint[2];
    if (scalexz) h[4] *= expfac;
    if (scaleyz) h[3] *= expfac;
  }

  // off-diagonal components, second half

  if (pstyle == TRICLINIC) {

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

    if (p_flag[3]) {
      expfac = exp(dto4*omega_dot[1]);
      h[3] *= expfac;
      h[3] += dto2*(omega_dot[3]*h[2]);
      h[3] *= expfac;
    }

    if (p_flag[5]) {
      expfac = exp(dto4*omega_dot[0]);
      h[5] *= expfac;
      h[5] += dto2*(omega_dot[5]*h[1]);
      h[5] *= expfac;
    }

    if (p_flag[4]) {
      expfac = exp(dto8*omega_dot[0]);
      h[4] *= expfac;
      h[4] += dto4*(omega_dot[5]*h[3]+omega_dot[4]*h[2]);
      h[4] *= expfac;
    }

  }

  domain->yz = h[3];
  domain->xz = h[4];
  domain->xy = h[5];

  // tilt factor to cell length ratio can not exceed TILTMAX in one step

  if (domain->yz < -TILTMAX*domain->yprd ||
      domain->yz > TILTMAX*domain->yprd ||
      domain->xz < -TILTMAX*domain->xprd ||
      domain->xz > TILTMAX*domain->xprd ||
      domain->xy < -TILTMAX*domain->xprd ||
      domain->xy > TILTMAX*domain->xprd)
    error->all(FLERR,"Fix npt/nph has tilted box too far in one step - "
               "periodic cell is too far from equilibrium state");

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap) domain->lamda2x(nlocal);
  else {
    int i;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (int i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}


/* ----------------------------------------------------------------------
   perform half-step barostat scaling of velocities
-----------------------------------------------------------------------*/

void FixNHOMP::nh_v_press()
{
  const double factor0 = exp(-dt4*(omega_dot[0]+mtk_term2));
  const double factor1 = exp(-dt4*(omega_dot[1]+mtk_term2));
  const double factor2 = exp(-dt4*(omega_dot[2]+mtk_term2));
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (which == NOBIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i].x *= factor0;
        v[i].y *= factor1;
        v[i].z *= factor2;
        if (pstyle == TRICLINIC) {
          v[i].x += -dthalf*(v[i].y*omega_dot[5] + v[i].z*omega_dot[4]);
          v[i].y += -dthalf*v[i].z*omega_dot[3];
        }
        v[i].x *= factor0;
        v[i].y *= factor1;
        v[i].z *= factor2;
      }
    }
  } else if (which == BIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      double buf[3];
      if (mask[i] & groupbit) {
        temperature->remove_bias_thr(i,&v[i].x,buf);
        v[i].x *= factor0;
        v[i].y *= factor1;
        v[i].z *= factor2;
        if (pstyle == TRICLINIC) {
          v[i].x += -dthalf*(v[i].y*omega_dot[5] + v[i].z*omega_dot[4]);
          v[i].y += -dthalf*v[i].z*omega_dot[3];
        }
        v[i].x *= factor0;
        v[i].y *= factor1;
        v[i].z *= factor2;
        temperature->restore_bias_thr(i,&v[i].x,buf);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities
-----------------------------------------------------------------------*/

void FixNHOMP::nve_v()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (atom->rmass) {
    const double * _noalias const rmass = atom->rmass;
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        v[i].x += dtfm*f[i].x;
        v[i].y += dtfm*f[i].y;
        v[i].z += dtfm*f[i].z;
      }
    }
  } else {
    const double *_noalias const mass = atom->mass;
    const int * _noalias const type = atom->type;
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        v[i].x += dtfm*f[i].x;
        v[i].y += dtfm*f[i].y;
        v[i].z += dtfm*f[i].z;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions
-----------------------------------------------------------------------*/

void FixNHOMP::nve_x()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  const dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  // x update by full step only for atoms in group

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i].x += dtv * v[i].x;
      x[i].y += dtv * v[i].y;
      x[i].z += dtv * v[i].z;
    }
  }
}
/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of velocities
-----------------------------------------------------------------------*/

void FixNHOMP::nh_v_temp()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (which == NOBIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
      }
    }
  } else if (which == BIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      double buf[3];
      if (mask[i] & groupbit) {
        temperature->remove_bias_thr(i,&v[i].x,buf);
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
        temperature->restore_bias_thr(i,&v[i].x,buf);
      }
    }
  }
}

