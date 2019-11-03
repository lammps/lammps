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

#include "fix_precession_spin.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
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
  k1c = k2c = 0.0;
  nc1x = nc1y = nc1z = 0.0;
  nc2x = nc2y = nc2z = 0.0;
  nc3x = nc3y = nc3z = 0.0;

  zeeman_flag = aniso_flag = cubic_flag = 0;

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
    } else if (strcmp(arg[iarg],"cubic") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      cubic_flag = 1;
      k1c = force->numeric(FLERR,arg[iarg+1]);
      k2c = force->numeric(FLERR,arg[iarg+2]);
      nc1x = force->numeric(FLERR,arg[iarg+3]);
      nc1y = force->numeric(FLERR,arg[iarg+4]);
      nc1z = force->numeric(FLERR,arg[iarg+5]);
      nc2x = force->numeric(FLERR,arg[iarg+6]);
      nc2y = force->numeric(FLERR,arg[iarg+7]);
      nc2z = force->numeric(FLERR,arg[iarg+8]);
      nc3x = force->numeric(FLERR,arg[iarg+9]);
      nc3y = force->numeric(FLERR,arg[iarg+10]);
      nc3z = force->numeric(FLERR,arg[iarg+11]);
      iarg += 12;
    } else error->all(FLERR,"Illegal precession/spin command");
  }

  // normalize vectors

  double inorm;
  if (zeeman_flag) {
    inorm = 1.0/sqrt(nhx*nhx + nhy*nhy + nhz*nhz);
    nhx *= inorm;
    nhy *= inorm;
    nhz *= inorm;
  }

  if (aniso_flag) {
    inorm = 1.0/sqrt(nax*nax + nay*nay + naz*naz);
    nax *= inorm;
    nay *= inorm;
    naz *= inorm;
  }
  
  if (cubic_flag) {
    inorm = 1.0/sqrt(nc1x*nc1x + nc1y*nc1y + nc1z*nc1z);
    nc1x *= inorm;
    nc1y *= inorm;
    nc1z *= inorm;
    inorm = 1.0/sqrt(nc2x*nc2x + nc2y*nc2y + nc2z*nc2z);
    nc2x *= inorm;
    nc2y *= inorm;
    nc2z *= inorm;
    inorm = 1.0/sqrt(nc3x*nc3x + nc3y*nc3y + nc3z*nc3z);
    nc3x *= inorm;
    nc3y *= inorm;
    nc3z *= inorm;
  }

  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  eprec = 0.0;
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

  // convert field quantities to rad.THz

  H_field *= gyro;
  Kah = Ka/hbar;
  k1ch = k1c/hbar;
  k2ch = k2c/hbar;

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

void FixPrecessionSpin::post_force(int /* vflag */)
{

  // update mag field with time (potential improvement)

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    modify->addstep_compute(update->ntimestep + 1);
    set_magneticprecession();           // update mag. field if time-dep.
  }

  int *mask = atom->mask;
  double **fm = atom->fm;
  double **sp = atom->sp;
  const int nlocal = atom->nlocal;
  double spi[3], fmi[3], epreci;
  
  eflag = 0;
  eprec = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      epreci = 0.0;
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
      fmi[0] = fmi[1] = fmi[2] = 0.0;

      if (zeeman_flag) {          // compute Zeeman interaction
        compute_zeeman(i,fmi);
        epreci -= hbar*(spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
      }

      if (aniso_flag) {           // compute magnetic anisotropy
        compute_anisotropy(spi,fmi);
        epreci -= compute_anisotropy_energy(spi);
      }

      if (cubic_flag) {         // compute cubic anisotropy
        compute_cubic(spi,fmi);
        epreci -= compute_cubic_energy(spi);
      }

      eprec += epreci;
      fm[i][0] += fmi[0];
      fm[i][1] += fmi[1];
      fm[i][2] += fmi[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::compute_single_precession(int i, double spi[3], double fmi[3])
{
  int *mask = atom->mask;
  if (mask[i] & groupbit) {
    if (zeeman_flag) compute_zeeman(i,fmi);
    if (aniso_flag) compute_anisotropy(spi,fmi);
    if (cubic_flag) compute_cubic(spi,fmi);
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

double FixPrecessionSpin::compute_anisotropy_energy(double spi[3])
{
  double energy = 0.0;
  double scalar = nax*spi[0] + nay*spi[1] + naz*spi[2];
  energy = Ka*scalar*scalar;
  return energy; 
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
    Kax = 2.0*Kah*nax;
    Kay = 2.0*Kah*nay;
    Kaz = 2.0*Kah*naz;
  }
}

/* ----------------------------------------------------------------------
   compute cubic aniso energy of spin i
------------------------------------------------------------------------- */

double FixPrecessionSpin::compute_cubic_energy(double spi[3])
{
  double energy = 0.0;
  double skx,sky,skz;

  skx = spi[0]*nc1x+spi[1]*nc1y+spi[2]*nc1z;
  sky = spi[0]*nc2x+spi[1]*nc2y+spi[2]*nc2z;
  skz = spi[0]*nc3x+spi[1]*nc3y+spi[2]*nc3z;

  energy = k1c*(skx*skx*sky*sky + sky*sky*skz*skz + skx*skx*skz*skz);
  energy += k2c*skx*skx*sky*sky*skz*skz;

  return energy;
}

/* ----------------------------------------------------------------------
   compute cubic anisotropy interaction for spin i
------------------------------------------------------------------------- */

void FixPrecessionSpin::compute_cubic(double spi[3], double fmi[3])
{
  double skx,sky,skz,skx2,sky2,skz2;
  double four1,four2,four3,fourx,foury,fourz;
  double six1,six2,six3,sixx,sixy,sixz;

  skx = spi[0]*nc1x+spi[1]*nc1y+spi[2]*nc1z;
  sky = spi[0]*nc2x+spi[1]*nc2y+spi[2]*nc2z;
  skz = spi[0]*nc3x+spi[1]*nc3y+spi[2]*nc3z;

  skx2 = skx*skx;
  sky2 = sky*sky;
  skz2 = skz*skz;

  four1 = 2.0*skx*(sky2+skz2);
  four2 = 2.0*sky*(skx2+skz2);
  four3 = 2.0*skz*(skx2+sky2);

  fourx = k1ch*(nc1x*four1 + nc2x*four2 + nc3x*four3);
  foury = k1ch*(nc1y*four1 + nc2y*four2 + nc3y*four3);
  fourz = k1ch*(nc1z*four1 + nc2z*four2 + nc3z*four3);

  six1 = 2.0*skx*sky2*skz2;
  six2 = 2.0*sky*skx2*skz2;
  six3 = 2.0*skz*skx2*sky2;
  
  sixx = k2ch*(nc1x*six1 + nc2x*six2 + nc3x*six3);
  sixy = k2ch*(nc1y*six1 + nc2y*six2 + nc3y*six3);
  sixz = k2ch*(nc1z*six1 + nc2z*six2 + nc3z*six3);
  
  fmi[0] += fourx + sixx;
  fmi[1] += foury + sixy;
  fmi[2] += fourz + sixz;
}

/* ----------------------------------------------------------------------
   potential energy in magnetic field
------------------------------------------------------------------------- */

double FixPrecessionSpin::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&eprec,&eprec_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return eprec_all;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::min_post_force(int vflag)
{
  post_force(vflag);
}
