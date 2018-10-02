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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "fix_nh_intel.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>
#include <cmath>
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;

#define TILTMAX 1.5

enum{NOBIAS,BIAS};
enum{ISO,ANISO,TRICLINIC};

typedef struct { double x,y,z; } dbl3_t;

/* ----------------------------------------------------------------------
   NVT,NPH,NPT integrators for improved Nose-Hoover equations of motion
 ---------------------------------------------------------------------- */

FixNHIntel::FixNHIntel(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  _dtfm = 0;
  _nlocal3 = 0;
  _nlocal_max = 0;
}

/* ---------------------------------------------------------------------- */

FixNHIntel::~FixNHIntel()
{
}

/* ---------------------------------------------------------------------- */

void FixNHIntel::setup(int vflag)
{
  FixNH::setup(vflag);
  reset_dt();
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or dilate group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixNHIntel::remap()
{
  double oldlo,oldhi;
  double expfac;

  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *h = domain->h;

  // omega is not used, except for book-keeping

  for (int i = 0; i < 6; i++) omega[i] += dto*omega_dot[i];

  // convert pertinent atoms and rigid bodies to lamda coords
  const double hi0 = domain->h_inv[0];
  const double hi1 = domain->h_inv[1];
  const double hi2 = domain->h_inv[2];
  const double hi3 = domain->h_inv[3];
  const double hi4 = domain->h_inv[4];
  const double hi5 = domain->h_inv[5];
  const double b0 = domain->boxlo[0];
  const double b1 = domain->boxlo[1];
  const double b2 = domain->boxlo[2];

  if (allremap) {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      const double d0 = x[i].x - b0;
      const double d1 = x[i].y - b1;
      const double d2 = x[i].z - b2;
      x[i].x = hi0*d0 + hi5*d1 + hi4*d2;
      x[i].y = hi1*d1 + hi3*d2;
      x[i].z = hi2*d2;
    }
  } else {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & dilate_group_bit) {
        const double d0 = x[i].x - b0;
        const double d1 = x[i].y - b1;
        const double d2 = x[i].z - b2;
        x[i].x = hi0*d0 + hi5*d1 + hi4*d2;
        x[i].y = hi1*d1 + hi3*d2;
        x[i].z = hi2*d2;
      }
    }
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

  double dto2 = dto/2.0;
  double dto4 = dto/4.0;
  double dto8 = dto/8.0;

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
  const double h0 = domain->h[0];
  const double h1 = domain->h[1];
  const double h2 = domain->h[2];
  const double h3 = domain->h[3];
  const double h4 = domain->h[4];
  const double h5 = domain->h[5];
  const double nb0 = domain->boxlo[0];
  const double nb1 = domain->boxlo[1];
  const double nb2 = domain->boxlo[2];

  if (allremap) {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      x[i].x = h0*x[i].x + h5*x[i].y + h4*x[i].z + nb0;
      x[i].y = h1*x[i].y + h3*x[i].z + nb1;
      x[i].z = h2*x[i].z + nb2;
    }
  } else {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & dilate_group_bit) {
        x[i].x = h0*x[i].x + h5*x[i].y + h4*x[i].z + nb0;
        x[i].y = h1*x[i].y + h3*x[i].z + nb1;
        x[i].z = h2*x[i].z + nb2;
      }
    }
  }

  if (nrigid)
    for (int i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ---------------------------------------------------------------------- */

void FixNHIntel::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  dto = dthalf;

  // If using respa, then remap is performed in innermost level

  if (strstr(update->integrate_style,"respa"))
    dto = 0.5*step_respa[0];

  if (pstat_flag)
    pdrag_factor = 1.0 - (update->dt * p_freq_max * drag / nc_pchain);

  if (tstat_flag)
    tdrag_factor = 1.0 - (update->dt * t_freq * drag / nc_tchain);

  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (nlocal > _nlocal_max) {
    if (_nlocal_max) memory->destroy(_dtfm);
    _nlocal_max = static_cast<int>(1.20 * nlocal);
    memory->create(_dtfm, _nlocal_max * 3, "fix_nve_intel:dtfm");
  }

  _nlocal3 = nlocal * 3;

  if (igroup == 0) {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      int n = 0;
      for (int i = 0; i < nlocal; i++) {
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
      }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = 0;
      for (int i = 0; i < nlocal; i++) {
        _dtfm[n++] = dtf / mass[type[i]];
        _dtfm[n++] = dtf / mass[type[i]];
        _dtfm[n++] = dtf / mass[type[i]];
      }
    }
  } else {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      int n = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          _dtfm[n++] = dtf / rmass[i];
          _dtfm[n++] = dtf / rmass[i];
          _dtfm[n++] = dtf / rmass[i];
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          _dtfm[n++] = dtf / mass[type[i]];
          _dtfm[n++] = dtf / mass[type[i]];
          _dtfm[n++] = dtf / mass[type[i]];
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step barostat scaling of velocities
-----------------------------------------------------------------------*/

void FixNHIntel::nh_v_press()
{
  if (pstyle == TRICLINIC || which == BIAS) {
    FixNH::nh_v_press();
    return;
  }

  dbl3_t * _noalias const v = (dbl3_t *)atom->v[0];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double f0 = exp(-dt4*(omega_dot[0]+mtk_term2));
  double f1 = exp(-dt4*(omega_dot[1]+mtk_term2));
  double f2 = exp(-dt4*(omega_dot[2]+mtk_term2));
  f0 *= f0;
  f1 *= f1;
  f2 *= f2;

  if (igroup == 0) {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      v[i].x *= f0;
      v[i].y *= f1;
      v[i].z *= f2;
    }
  } else {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i].x *= f0;
        v[i].y *= f1;
        v[i].z *= f2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities
-----------------------------------------------------------------------*/

void FixNHIntel::nve_v()
{
  if (neighbor->ago == 0) reset_dt();
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  #if defined(LMP_SIMD_COMPILER)
  #pragma vector aligned
  #pragma simd
  #endif
  for (int i = 0; i < _nlocal3; i++)
    v[i] += _dtfm[i] * f[i];
}

/* ----------------------------------------------------------------------
   perform full-step update of positions
-----------------------------------------------------------------------*/

void FixNHIntel::nve_x()
{
  double * _noalias const x = atom->x[0];
  double * _noalias const v = atom->v[0];

  // x update by full step only for atoms in group

  if (igroup == 0) {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++)
      x[i] += dtv * v[i];
  } else {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      if (_dtfm[i] != 0.0)
        x[i] += dtv * v[i];
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of velocities
-----------------------------------------------------------------------*/

void FixNHIntel::nh_v_temp()
{
  if (which == BIAS) {
    FixNH::nh_v_temp();
    return;
  }

  double * _noalias const v = atom->v[0];

  if (igroup == 0) {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++)
        v[i] *= factor_eta;
  } else {
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      if (_dtfm[i] != 0.0)
        v[i] *= factor_eta;
    }
  }
}

double FixNHIntel::memory_usage()
{
  return FixNH::memory_usage() + _nlocal_max * 3 * sizeof(double);
}
