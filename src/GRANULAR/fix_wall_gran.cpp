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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{XPLANE,YPLANE,ZPLANE,ZCYLINDER};
enum{NO_HISTORY,HISTORY,HERTZIAN};

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix wall/gran command");

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Fix wall/gran requires atom attributes radius, omega, torque");

  restart_peratom = 1;
  
  int iarg;
  if (strcmp(arg[3],"xplane") == 0) {
    iarg = 8;
    if (narg < iarg) error->all("Illegal fix wall/gran command");
    wallstyle = XPLANE;
    if (strcmp(arg[4],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[4]);
    if (strcmp(arg[5],"NULL") == 0) hi = BIG;
    else hi = atof(arg[5]);
    gamman = atof(arg[6]);
    xmu = atof(arg[7]);
  } else if (strcmp(arg[3],"yplane") == 0) {
    iarg = 8;
    if (narg < iarg) error->all("Illegal fix wall/gran command");
    wallstyle = YPLANE;
    if (strcmp(arg[4],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[4]);
    if (strcmp(arg[5],"NULL") == 0) hi = BIG;
    else hi = atof(arg[5]);
    gamman = atof(arg[6]);
    xmu = atof(arg[7]);
  } else if (strcmp(arg[3],"zplane") == 0) {
    iarg = 8;
    if (narg < iarg) error->all("Illegal fix wall/gran command");
    wallstyle = ZPLANE;
    if (strcmp(arg[4],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[4]);
    if (strcmp(arg[5],"NULL") == 0) hi = BIG;
    else hi = atof(arg[5]);
    gamman = atof(arg[6]);
    xmu = atof(arg[7]);
  } else if (strcmp(arg[3],"zcylinder") == 0) {
    iarg = 7;
    if (narg < iarg) error->all("Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = atof(arg[4]);
    gamman = atof(arg[5]);
    xmu = atof(arg[6]);
  }
  
  // check for trailing keyword/values

  wiggle = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all("Illegal fix wall/gran command");
      amplitude = atof(arg[iarg+2]);
      period = atof(arg[iarg+3]);
      wiggle = 1;
      iarg += 4;
    } else error->all("Illegal fix wall/gran command");
  }

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all("Cannot use wall in periodic dimension");

  if (wallstyle == ZCYLINDER && wiggle)
    if (axis != 2) error->all("Can only wiggle zcylinder wall in z dim");

  // setup oscillations

  if (wiggle) {
    double PI = 4.0 * atan(1.0);
    omega = 2.0*PI / period;
    time_origin = update->ntimestep;
  }

  // perform initial allocation of atom-based arrays
  // register with Atom class

  shear = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize as if particle is not touching wall

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    shear[i][0] = shear[i][1] = shear[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy_2d_double_array(shear);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
  // set local values from Pair values

  if (force->pair == NULL)
    error->all("Fix wall/gran is incompatible with Pair style");
  double *p_xkk = (double *) force->pair->extract("xkk");
  if (!p_xkk)  error->all("Fix wall/gran is incompatible with Pair style");
  xkk = *p_xkk;

  // same initialization as in pair_gran_history::init_style()

  xkkt = xkk * 2.0/7.0;
  dt = update->dt;
  double gammas = 0.5*gamman;
  gamman_dl = gamman/dt;
  gammas_dl = gammas/dt;

  // set pairstyle from granular pair style

  if (force->pair_match("gran/no_history")) pairstyle = NO_HISTORY;
  else if (force->pair_match("gran/history")) pairstyle = HISTORY;
  else if (force->pair_match("gran/hertzian")) pairstyle = HERTZIAN;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force(int vflag)
{
  double vwall[3],dx,dy,dz,del1,del2,delxy,delr,rsq;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    wlo = lo + amplitude - amplitude*cos(arg);
    whi = hi + amplitude - amplitude*cos(arg);
    vwall[axis] = dt * amplitude*omega*sin(arg);
  }

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  //   in cylinder case
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set shear if pair potential stores history

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
	del1 = x[i][0] - wlo;
	del2 = whi - x[i][0];
	if (del1 < del2) dx = del1;
	else dx = -del2;
      } else if (wallstyle == YPLANE) {
	del1 = x[i][1] - wlo;
	del2 = whi - x[i][1];
	if (del1 < del2) dy = del1;
	else dy = -del2;
      } else if (wallstyle == ZPLANE) {
	del1 = x[i][2] - wlo;
	del2 = whi - x[i][2];
	if (del1 < del2) dz = del1;
	else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
	delr = cylradius - delxy;
	if (delr > radius[i]) dz = cylradius;
	else {
	  dx = -delr/delxy * x[i][0];
	  dy = -delr/delxy * x[i][1];
	}
      }

      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq > radius[i]*radius[i]) {
	if (pairstyle != NO_HISTORY) {
	  shear[i][0] = 0.0;
	  shear[i][1] = 0.0;
	  shear[i][2] = 0.0;
	}
      } else {
	if (pairstyle == NO_HISTORY)
	  no_history(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
		     radius[i],rmass[i]);
	else if (pairstyle == HISTORY)
	  history(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
		  radius[i],rmass[i],shear[i]);
	else if (pairstyle == HERTZIAN)
	  hertzian(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
		   radius[i],rmass[i],shear[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::no_history(double rsq, double dx, double dy, double dz,
			     double *vwall, double *v,
			     double *f, double *omega, double *torque,
			     double radius, double mass)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,xmeff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,ccelx,ccely,ccelz,tor1,tor2,tor3;

  r = sqrt(rsq);

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  vr1 *= dt;
  vr2 *= dt;
  vr3 *= dt;

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr / rsq;
  vn2 = dy*vnnr / rsq;
  vn3 = dz*vnnr / rsq;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0];
  wr2 = radius*omega[1];
  wr3 = radius*omega[2];

  wr1 *= dt/r;
  wr2 *= dt/r;
  wr3 *= dt/r;

  // normal damping term
  // this definition of DAMP includes the extra 1/r term

  xmeff = mass;
  damp = xmeff*gamman_dl*vnnr/rsq;
  ccel = xkk*(radius-r)/r - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = xmeff*gammas_dl*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // shear friction forces

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

  // force components

  ccelx = dx*ccel + fs1;
  ccely = dy*ccel + fs2;
  ccelz = dz*ccel + fs3;

  // forces

  f[0] += ccelx;
  f[1] += ccely;
  f[2] += ccelz;

  // torques

  tor1 = dy*fs3 - dz*fs2;
  tor2 = dz*fs1 - dx*fs3;
  tor3 = dx*fs2 - dy*fs1;
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::history(double rsq, double dx, double dy, double dz,
			  double *vwall, double *v,
			  double *f, double *omega, double *torque,
			  double radius, double mass, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,xmeff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,ccelx,ccely,ccelz,tor1,tor2,tor3;
  double shrmag,rsht,rinv;

  r = sqrt(rsq);

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  vr1 *= dt;
  vr2 *= dt;
  vr3 *= dt;

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr / rsq;
  vn2 = dy*vnnr / rsq;
  vn3 = dz*vnnr / rsq;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0];
  wr2 = radius*omega[1];
  wr3 = radius*omega[2];

  wr1 *= dt/r;
  wr2 *= dt/r;
  wr3 *= dt/r;

  // normal damping term
  // this definition of DAMP includes the extra 1/r term

  xmeff = mass;
  damp = xmeff*gamman_dl*vnnr/rsq;
  ccel = xkk*(radius-r)/r - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  shear[0] += vtr1;
  shear[1] += vtr2;
  shear[2] += vtr3;
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements correctly

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht/rsq;
  shear[0] -= rsht*dx;
  shear[1] -= rsht*dy;
  shear[2] -= rsht*dz;

  // tangential forces

  fs1 = - (xkkt*shear[0] + xmeff*gammas_dl*vtr1);
  fs2 = - (xkkt*shear[1] + xmeff*gammas_dl*vtr2);
  fs3 = - (xkkt*shear[2] + xmeff*gammas_dl*vtr3);

  // force normalization

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  // shrmag is magnitude of shearwall
  // rescale frictional displacements and forces if needed

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + xmeff*gammas_dl*vtr1/xkkt) - 
	xmeff*gammas_dl*vtr1/xkkt;
      shear[1] = (fn/fs) * (shear[1] + xmeff*gammas_dl*vtr2/xkkt) -
	xmeff*gammas_dl*vtr2/xkkt;
      shear[2] = (fn/fs) * (shear[2] + xmeff*gammas_dl*vtr3/xkkt) -
	xmeff*gammas_dl*vtr3/xkkt;
      fs1 = fs1 * fn / fs ;
      fs2 = fs2 * fn / fs;
      fs3 = fs3 * fn / fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  ccelx = dx*ccel + fs1;
  ccely = dy*ccel + fs2;
  ccelz = dz*ccel + fs3;

  // forces

  f[0] += ccelx;
  f[1] += ccely;
  f[2] += ccelz;

  // torques

  rinv = 1/r;
  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hertzian(double rsq, double dx, double dy, double dz,
			   double *vwall, double *v,
			   double *f, double *omega, double *torque,
			   double radius, double mass, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,xmeff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,ccelx,ccely,ccelz,tor1,tor2,tor3;
  double shrmag,rsht,rhertz;

  r = sqrt(rsq);

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  vr1 *= dt;
  vr2 *= dt;
  vr3 *= dt;

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr / rsq;
  vn2 = dy*vnnr / rsq;
  vn3 = dz*vnnr / rsq;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0];
  wr2 = radius*omega[1];
  wr3 = radius*omega[2];

  wr1 *= dt/r;
  wr2 *= dt/r;
  wr3 *= dt/r;

  // normal damping term
  // this definition of DAMP includes the extra 1/r term

  xmeff = mass;
  damp = xmeff*gamman_dl*vnnr/rsq;
  ccel = xkk*(radius-r)/r - damp;
  rhertz = sqrt(radius - r);
  ccel = rhertz * ccel;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  shear[0] += vtr1;
  shear[1] += vtr2;
  shear[2] += vtr3;
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements correctly

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht/rsq;
  shear[0] -= rsht*dx;
  shear[1] -= rsht*dy;
  shear[2] -= rsht*dz;

  // tangential forces

  fs1 = -rhertz * (xkkt*shear[0] + xmeff*gammas_dl*vtr1);
  fs2 = -rhertz * (xkkt*shear[1] + xmeff*gammas_dl*vtr2);
  fs3 = -rhertz * (xkkt*shear[2] + xmeff*gammas_dl*vtr3);

  // force normalization

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  // shrmag is magnitude of shearwall
  // rescale frictional displacements and forces if needed

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + xmeff*gammas_dl*vtr1/xkkt) - 
	xmeff*gammas_dl*vtr1/xkkt;
      shear[1] = (fn/fs) * (shear[1] + xmeff*gammas_dl*vtr2/xkkt) -
	xmeff*gammas_dl*vtr2/xkkt;
      shear[2] = (fn/fs) * (shear[2] + xmeff*gammas_dl*vtr3/xkkt) -
	xmeff*gammas_dl*vtr3/xkkt;
      fs1 = fs1 * fn / fs ;
      fs2 = fs2 * fn / fs;
      fs3 = fs3 * fn / fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  ccelx = dx*ccel + fs1;
  ccely = dy*ccel + fs2;
  ccelz = dz*ccel + fs3;

  // forces

  f[0] += ccelx;
  f[1] += ccely;
  f[2] += ccelz;

  // torques

  tor1 = dy*fs3 - dz*fs2;
  tor2 = dz*fs1 - dx*fs3;
  tor3 = dx*fs2 - dy*fs1;
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixWallGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += 3*nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  shear = memory->grow_2d_double_array(shear,nmax,3,"fix_wall_gran:shear");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixWallGran::copy_arrays(int i, int j)
{
  shear[j][0] = shear[i][0];
  shear[j][1] = shear[i][1];
  shear[j][2] = shear[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixWallGran::pack_exchange(int i, double *buf)
{
  buf[0] = shear[i][0];
  buf[1] = shear[i][1];
  buf[2] = shear[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange 
------------------------------------------------------------------------- */

int FixWallGran::unpack_exchange(int nlocal, double *buf)
{
  shear[nlocal][0] = buf[0];
  shear[nlocal][1] = buf[1];
  shear[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file 
------------------------------------------------------------------------- */

int FixWallGran::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 4;
  buf[m++] = shear[i][0];
  buf[m++] = shear[i][1];
  buf[m++] = shear[i][2];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix 
------------------------------------------------------------------------- */

void FixWallGran::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  shear[nlocal][0] = extra[nlocal][m++];
  shear[nlocal][1] = extra[nlocal][m++];
  shear[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data 
------------------------------------------------------------------------- */

int FixWallGran::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data 
------------------------------------------------------------------------- */

int FixWallGran::size_restart(int nlocal)
{
  return 4;
}
