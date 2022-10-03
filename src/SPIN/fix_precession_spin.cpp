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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_precession_spin.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "error.h"
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

FixPrecessionSpin::FixPrecessionSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), emag(nullptr)
{
  if (narg < 7) error->all(FLERR,"Illegal precession/spin command");

  // magnetic interactions coded for cartesian coordinates

  hbar = force->hplanck/MY_2PI;

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  magstr = nullptr;
  magfieldstyle = CONSTANT;

  H_field = 0.0;
  nhx = nhy = nhz = 0.0;
  hx = hy = hz = 0.0;
  stt_field = 0.0;
  nsttx = nstty = nsttz = 0.0;
  sttx = stty = sttz = 0.0;
  Ka = 0.0;
  nax = nay = naz = 0.0;
  Kax = Kay = Kaz = 0.0;
  k1c = k2c = 0.0;
  nc1x = nc1y = nc1z = 0.0;
  nc2x = nc2y = nc2z = 0.0;
  nc3x = nc3y = nc3z = 0.0;
  K6 = 0.0;
  n6x = n6y = n6z = 0.0;
  m6x = m6y = m6z = 0.0;

  zeeman_flag = stt_flag = aniso_flag = cubic_flag = hexaniso_flag =  0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"zeeman") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      zeeman_flag = 1;
      H_field = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      nhx = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      nhy = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      nhz = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;
    } else if (strcmp(arg[iarg],"stt") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      stt_flag = 1;
      stt_field = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      nsttx = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      nstty = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      nsttz = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;
    } else if (strcmp(arg[iarg],"anisotropy") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      aniso_flag = 1;
      Ka = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      nax = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      nay = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      naz = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;
    } else if (strcmp(arg[iarg],"cubic") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      cubic_flag = 1;
      k1c = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      k2c = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      nc1x = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      nc1y = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      nc1z = utils::numeric(FLERR,arg[iarg+5],false,lmp);
      nc2x = utils::numeric(FLERR,arg[iarg+6],false,lmp);
      nc2y = utils::numeric(FLERR,arg[iarg+7],false,lmp);
      nc2z = utils::numeric(FLERR,arg[iarg+8],false,lmp);
      nc3x = utils::numeric(FLERR,arg[iarg+9],false,lmp);
      nc3y = utils::numeric(FLERR,arg[iarg+10],false,lmp);
      nc3z = utils::numeric(FLERR,arg[iarg+11],false,lmp);
      iarg += 12;
    } else if (strcmp(arg[iarg],"hexaniso") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal fix precession/spin command");
      hexaniso_flag = 1;
      K6 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      n6x = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      n6y = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      n6z = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      m6x = utils::numeric(FLERR,arg[iarg+5],false,lmp);
      m6y = utils::numeric(FLERR,arg[iarg+6],false,lmp);
      m6z = utils::numeric(FLERR,arg[iarg+7],false,lmp);
      iarg += 8;
    } else error->all(FLERR,"Illegal precession/spin command");
  }

  // normalize vectors

  double norm2,inorm;
  if (zeeman_flag) {
    norm2 = nhx*nhx + nhy*nhy + nhz*nhz;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nhx *= inorm;
    nhy *= inorm;
    nhz *= inorm;
  }

  if (stt_flag) {
    norm2 = nsttx*nsttx + nstty*nstty + nsttz*nsttz;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nsttx *= inorm;
    nstty *= inorm;
    nsttz *= inorm;
  }

  if (aniso_flag) {
    norm2 = nax*nax + nay*nay + naz*naz;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nax *= inorm;
    nay *= inorm;
    naz *= inorm;
  }

  if (cubic_flag) {
    norm2 = nc1x*nc1x + nc1y*nc1y + nc1z*nc1z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nc1x *= inorm;
    nc1y *= inorm;
    nc1z *= inorm;

    norm2 = nc2x*nc2x + nc2y*nc2y + nc2z*nc2z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nc2x *= inorm;
    nc2y *= inorm;
    nc2z *= inorm;

    norm2 = nc3x*nc3x + nc3y*nc3y + nc3z*nc3z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    nc3x *= inorm;
    nc3y *= inorm;
    nc3z *= inorm;
  }

  if (hexaniso_flag) {
    norm2 = n6x*n6x + n6y*n6y + n6z*n6z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    n6x *= inorm;
    n6y *= inorm;
    n6z *= inorm;

    norm2 = m6x*m6x + m6y*m6y + m6z*m6z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    m6x *= inorm;
    m6y *= inorm;
    m6z *= inorm;
    l6x = (n6z*m6y-n6y*m6z);
    l6y = (n6x*m6z-n6z*m6x);
    l6z = (n6y*m6x-n6x*m6y);

    norm2 = l6x*l6x + l6y*l6y + l6z*l6z;
    if (norm2 == 0.0)
      error->all(FLERR,"Illegal precession/spin command");
    inorm = 1.0/sqrt(norm2);
    l6x *= inorm;
    l6y *= inorm;
    l6z *= inorm;
    m6x = (l6z*n6y-l6y*n6z);
    m6y = (l6x*n6z-l6z*n6x);
    m6z = (l6y*n6x-l6x*n6y);
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
  memory->destroy(emag);
}

/* ---------------------------------------------------------------------- */

int FixPrecessionSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::init()
{
  const double hbar = force->hplanck/MY_2PI;    // eV/(rad.THz)
  const double mub = 5.78901e-5;                // in eV/T
  const double gyro = 2.0*mub/hbar;             // in rad.THz/T

  // convert field quantities to rad.THz

  H_field *= gyro;
  Kah = Ka/hbar;
  k1ch = k1c/hbar;
  k2ch = k2c/hbar;
  K6h = K6/hbar;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
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

  // init. size of energy stacking lists

  nlocal_max = atom->nlocal;
  memory->grow(emag,nlocal_max,"pair/spin:emag");
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
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
  double spi[4], fmi[3], epreci;

  // checking size of emag

  if (nlocal_max < nlocal) {                    // grow emag lists if necessary
    nlocal_max = nlocal;
    memory->grow(emag,nlocal_max,"pair/spin:emag");
  }

  eflag = 0;
  eprec = 0.0;
  for (int i = 0; i < nlocal; i++) {
    emag[i] = 0.0;
    if (mask[i] & groupbit) {
      epreci = 0.0;
      spi[0] = sp[i][0];
      spi[1] = sp[i][1];
      spi[2] = sp[i][2];
      spi[3] = sp[i][3];
      fmi[0] = fmi[1] = fmi[2] = 0.0;

      if (zeeman_flag) {          // compute Zeeman interaction
        compute_zeeman(i,fmi);
        epreci -= compute_zeeman_energy(spi);
      }

      if (stt_flag) {             // compute Spin Transfer Torque
        compute_stt(spi,fmi);
        epreci -= compute_stt_energy(spi);
      }

      if (aniso_flag) {           // compute magnetic anisotropy
        compute_anisotropy(spi,fmi);
        epreci -= compute_anisotropy_energy(spi);
      }

      if (cubic_flag) {           // compute cubic anisotropy
        compute_cubic(spi,fmi);
        epreci -= compute_cubic_energy(spi);
      }

      if (hexaniso_flag) {        // compute hexagonal anisotropy
        compute_hexaniso(spi,fmi);
        epreci -= compute_hexaniso_energy(spi);
      }

      emag[i] += epreci;
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
    if (stt_flag) compute_stt(spi,fmi);
    if (aniso_flag) compute_anisotropy(spi,fmi);
    if (cubic_flag) compute_cubic(spi,fmi);
    if (hexaniso_flag) compute_hexaniso(spi,fmi);
  }
}

/* ----------------------------------------------------------------------
   Zeeman
------------------------------------------------------------------------- */

void FixPrecessionSpin::compute_zeeman(int i, double fmi[3])
{
  double **sp = atom->sp;
  fmi[0] += sp[i][3]*hx;
  fmi[1] += sp[i][3]*hy;
  fmi[2] += sp[i][3]*hz;
}

/* ---------------------------------------------------------------------- */

double FixPrecessionSpin::compute_zeeman_energy(double spi[4])
{
  double energy = 0.0;
  double scalar = nhx*spi[0]+nhy*spi[1]+nhz*spi[2];
  energy = hbar*H_field*spi[3]*scalar;
  return energy;
}

/* ----------------------------------------------------------------------
   STT
------------------------------------------------------------------------- */

void FixPrecessionSpin::compute_stt(double spi[3], double fmi[3])
{
  double sx = spi[0];
  double sy = spi[1];
  double sz = spi[2];
  fmi[0] += 1.0*stt_field*( sy*nsttz-sz*nstty);
  fmi[1] += 1.0*stt_field*(-sx*nsttz+sz*nsttx);
  fmi[2] += 1.0*stt_field*( sx*nstty-sy*nsttx);
}

/* ---------------------------------------------------------------------- */

double FixPrecessionSpin::compute_stt_energy(double * /* spi */)
{
  double energy = 0.0;  // Non-conservative force
  return energy;
}

/* ----------------------------------------------------------------------
   compute uniaxial anisotropy interaction for spin i
------------------------------------------------------------------------- */

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

  fmi[0] += (fourx + sixx);
  fmi[1] += (foury + sixy);
  fmi[2] += (fourz + sixz);
}

/* ---------------------------------------------------------------------- */

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
   compute hexagonal anisotropy interaction for spin i
------------------------------------------------------------------------- */

void FixPrecessionSpin::compute_hexaniso(double spi[3], double fmi[3])
{
  double s_x,s_y;
  double pf, phi, ssint2;

  // changing to the axes' frame

  s_x = l6x*spi[0]+l6y*spi[1]+l6z*spi[2];
  s_y = m6x*spi[0]+m6y*spi[1]+m6z*spi[2];

  // hexagonal anisotropy in the axes' frame

  phi = atan2(s_y,s_x);
  ssint2 = s_x*s_x + s_y*s_y;                 // s^2sin^2(theta)
  pf = 6.0 * K6h * ssint2*ssint2*sqrt(ssint2);   // 6*K_6*s^5*sin^5(theta)
  double fm_x =  pf*cos(5*phi);
  double fm_y = -pf*sin(5*phi);
  double fm_z =  0;

  // back to the lab's frame

  fmi[0] += fm_x*l6x+fm_y*m6x+fm_z*n6x;
  fmi[1] += fm_x*l6y+fm_y*m6y+fm_z*n6y;
  fmi[2] += fm_x*l6z+fm_y*m6z+fm_z*n6z;
}

/* ----------------------------------------------------------------------
   compute hexagonal aniso energy of spin i
------------------------------------------------------------------------- */

double FixPrecessionSpin::compute_hexaniso_energy(double spi[3])
{
  double energy = 0.0;
  double s_x,s_y,s_z, phi,ssint2;

  // changing to the axes' frame

  s_x = l6x*spi[0]+l6y*spi[1]+l6z*spi[2];
  s_y = m6x*spi[0]+m6y*spi[1]+m6z*spi[2];
  s_z = n6x*spi[0]+n6y*spi[1]+n6z*spi[2];

  // hexagonal anisotropy in the axes' frame

  phi = atan2(s_y,s_z);
  ssint2 = s_x*s_x + s_y*s_y;

  energy = K6 * ssint2*ssint2*ssint2*cos(6*phi);

  return 2.0*energy;
}

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::set_magneticprecession()
{
  if (zeeman_flag) {
    hx = H_field*nhx;
    hy = H_field*nhy;
    hz = H_field*nhz;
  }

  if (stt_flag) {
    sttx = stt_field*nsttx;
    stty = stt_field*nstty;
    sttz = stt_field*nsttz;
  }

  if (aniso_flag) {
    Kax = 2.0*Kah*nax;
    Kay = 2.0*Kah*nay;
    Kaz = 2.0*Kah*naz;
  }
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

/* ---------------------------------------------------------------------- */

void FixPrecessionSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}
