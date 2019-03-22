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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_precession_spin.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixPrecessionSpin::FixPrecessionSpin(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal precession/spin command");

  // magnetic interactions coded for cartesian coordinates

  hbar = force->hplanck/MY_2PI;

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  magstr = NULL;
  magfieldstyle = CONSTANT;

  H_field = 0.0;
  nhx = nhy = nhz = 0.0;
  hx = hy = hz = 0.0;
  Ka = 0.0;
  nax = nay = naz = 0.0;
  Kax = Kay = Kaz = 0.0;

  zeeman_flag = aniso_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"zeeman") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      zeeman_flag = 1;
      H_field = force->numeric(FLERR,arg[iarg+1]);
      nhx = force->numeric(FLERR,arg[iarg+2]);
      nhy = force->numeric(FLERR,arg[iarg+3]);
      nhz = force->numeric(FLERR,arg[iarg+4]);
      iarg += 5;
    } else if (strcmp(arg[iarg],"anisotropy") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      aniso_flag = 1;
      Ka = force->numeric(FLERR,arg[iarg+1]);
      nax = force->numeric(FLERR,arg[iarg+2]);
      nay = force->numeric(FLERR,arg[iarg+3]);
      naz = force->numeric(FLERR,arg[iarg+4]);
      iarg += 5;
    } else error->all(FLERR,"Illegal precession/spin command");
  }

  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  emag = 0.0;
}

/* ---------------------------------------------------------------------- */

FixPrecessionSpin::~FixPrecessionSpin()
{
  delete [] magstr;
}

/* ---------------------------------------------------------------------- */

int FixPrecessionSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}


/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::init()
{
  const double hbar = force->hplanck/MY_2PI;    // eV/(rad.THz)
  const double mub = 5.78901e-5;                // in eV/T
  const double gyro = mub/hbar;                 // in rad.THz/T

  H_field *= gyro;                              // in rad.THz
  Ka /= hbar;                                   // in rad.THz

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  if (magstr) {
  magvar = input->variable->find(magstr);
  if (magvar < 0)
        error->all(FLERR,"Illegal precession/spin command");
  if (!input->variable->equalstyle(magvar))
        error->all(FLERR,"Illegal precession/spin command");
  }

  varflag = CONSTANT;
  if (magfieldstyle != CONSTANT) varflag = EQUAL;

  // set magnetic field components

  if (varflag == CONSTANT) set_magneticprecession();

}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::post_force(int vflag)
{

  // update mag field with time (potential improvement)

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    modify->addstep_compute(update->ntimestep + 1);
    set_magneticprecession();                   // update mag. field if time-dep.
  }

  double **sp = atom->sp;
  double **fm = atom->fm;
  double spi[3], fmi[3];
  const int nlocal = atom->nlocal;

  eflag = 0;
  emag = 0.0;

  for (int i = 0; i < nlocal; i++) {
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];
    fmi[0] = fmi[1] = fmi[2] = 0.0;

    if (zeeman_flag) {          // compute Zeeman interaction
      compute_zeeman(i,fmi);
      emag -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
    }

    if (aniso_flag) {           // compute magnetic anisotropy
      compute_anisotropy(spi,fmi);
      emag -= 0.5*(spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
    }

    fm[i][0] += fmi[0];
    fm[i][1] += fmi[1];
    fm[i][2] += fmi[2];
  }
  emag *= hbar;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::compute_single_precession(int i, double spi[3], double fmi[3])
{
  if (zeeman_flag) {
    compute_zeeman(i,fmi);
  }
  if (aniso_flag) {
    compute_anisotropy(spi,fmi);
  }
}


/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::compute_zeeman(int i, double fmi[3])
{
  double **sp = atom->sp;
  fmi[0] += sp[i][3]*hx;
  fmi[1] += sp[i][3]*hy;
  fmi[2] += sp[i][3]*hz;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::compute_anisotropy(double spi[3], double fmi[3])
{
  double scalar = nax*spi[0] + nay*spi[1] + naz*spi[2];
  fmi[0] += scalar*Kax;
  fmi[1] += scalar*Kay;
  fmi[2] += scalar*Kaz;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::set_magneticprecession()
{
  if (zeeman_flag) {
          hx = H_field*nhx;
          hy = H_field*nhy;
          hz = H_field*nhz;
  }
  if (aniso_flag) {
          Kax = 2.0*Ka*nax;
          Kay = 2.0*Ka*nay;
          Kaz = 2.0*Ka*naz;
  }
}

/* ----------------------------------------------------------------------
   potential energy in magnetic field
------------------------------------------------------------------------- */

double FixPrecessionSpin::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&emag,&emag_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return emag_all;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::min_post_force(int vflag)
{
  post_force(vflag);
}
