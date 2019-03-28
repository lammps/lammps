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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL),
                         Dan Bolintineanu (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// XYZ PLANE need to be 0,1,2

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,BONDED_HISTORY};
enum{NONE,CONSTANT,EQUAL};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(NULL), shearone(NULL), fix_rigid(NULL), mass_rigid(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix wall/gran command");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix wall/gran requires atom style sphere");

  create_attribute = 1;

  // set interaction style
  // disable bonded/history option for now

  if (strcmp(arg[3],"hooke") == 0) pairstyle = HOOKE;
  else if (strcmp(arg[3],"hooke/history") == 0) pairstyle = HOOKE_HISTORY;
  else if (strcmp(arg[3],"hertz/history") == 0) pairstyle = HERTZ_HISTORY;
  //else if (strcmp(arg[3],"bonded/history") == 0) pairstyle = BONDED_HISTORY;
  else error->all(FLERR,"Invalid fix wall/gran interaction style");

  history = restart_peratom = 1;
  if (pairstyle == HOOKE) history = restart_peratom = 0;

  // wall/particle coefficients

  int iarg;

  if (pairstyle != BONDED_HISTORY) {
    if (narg < 11) error->all(FLERR,"Illegal fix wall/gran command");

    kn = force->numeric(FLERR,arg[4]);
    if (strcmp(arg[5],"NULL") == 0) kt = kn * 2.0/7.0;
    else kt = force->numeric(FLERR,arg[5]);

    gamman = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) gammat = 0.5 * gamman;
    else gammat = force->numeric(FLERR,arg[7]);

    xmu = force->numeric(FLERR,arg[8]);
    int dampflag = force->inumeric(FLERR,arg[9]);
    if (dampflag == 0) gammat = 0.0;

    if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
        xmu < 0.0 || xmu > 10000.0 || dampflag < 0 || dampflag > 1)
      error->all(FLERR,"Illegal fix wall/gran command");

    // convert Kn and Kt from pressure units to force/distance^2 if Hertzian

    if (pairstyle == HERTZ_HISTORY) {
      kn /= force->nktv2p;
      kt /= force->nktv2p;
    }

    iarg = 10;
  }

  else {
    if (narg < 10) error->all(FLERR,"Illegal fix wall/gran command");

    E = force->numeric(FLERR,arg[4]);
    G = force->numeric(FLERR,arg[5]);
    SurfEnergy = force->numeric(FLERR,arg[6]);
    // Note: this doesn't get used, check w/ Jeremy?
    gamman = force->numeric(FLERR,arg[7]);

    xmu = force->numeric(FLERR,arg[8]);
    // pois = E/(2.0*G) - 1.0;
    // kn = 2.0*E/(3.0*(1.0+pois)*(1.0-pois));
    // gammat=0.5*gamman;

    iarg = 9;
  }

  // wallstyle args

  idregion = NULL;

  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
  } else if (strcmp(arg[iarg],"region") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = REGION;
    int n = strlen(arg[iarg+1]) + 1;
    idregion = new char[n];
    strcpy(idregion,arg[iarg+1]);
    iarg += 2;
  }

  // optional args

  wiggle = 0;
  wshear = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      amplitude = force->numeric(FLERR,arg[iarg+2]);
      period = force->numeric(FLERR,arg[iarg+3]);
      wiggle = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"shear") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      vshear = force->numeric(FLERR,arg[iarg+2]);
      wshear = 1;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix wall/gran command");
  }

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all(FLERR,"Cannot use wall in periodic dimension");

  if (wiggle && wshear)
    error->all(FLERR,"Cannot wiggle and shear fix wall/gran");
  if (wiggle && wallstyle == ZCYLINDER && axis != 2)
    error->all(FLERR,"Invalid wiggle direction for fix wall/gran");
  if (wshear && wallstyle == XPLANE && axis == 0)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == YPLANE && axis == 1)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == ZPLANE && axis == 2)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if ((wiggle || wshear) && wallstyle == REGION)
    error->all(FLERR,"Cannot wiggle or shear with fix wall/gran/region");

  // setup oscillations

  if (wiggle) omega = 2.0*MY_PI / period;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  if (pairstyle == BONDED_HISTORY) sheardim = 7;
  else sheardim = 3;

  shearone = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  nmax = 0;
  mass_rigid = NULL;

  // initialize shear history as if particle is not touching region
  // shearone will be NULL for wallstyle = REGION

  if (history && shearone) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < sheardim; j++)
        shearone[i][j] = 0.0;
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
  if (copymode) return;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete local storage

  delete [] idregion;
  memory->destroy(shearone);
  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
  int i;

  dt = update->dt;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = NULL;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force(int /*vflag*/)
{
  int i,j;
  double dx,dy,dz,del1,del2,delxy,delr,rsq,rwall,meff;
  double vwall[3];

  // do not update shear history during setup

  shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"wall/gran:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    }
  }

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
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

  rwall = 0.0;

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
        if (delr > radius[i]) {
          dz = cylradius;
          rwall = 0.0;
        } else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
          // rwall = -2r_c if inside cylinder, 2r_c outside
          rwall = (delxy < cylradius) ? -2*cylradius : 2*cylradius;
          if (wshear && axis != 2) {
            vwall[0] += vshear * x[i][1]/delxy;
            vwall[1] += -vshear * x[i][0]/delxy;
            vwall[2] = 0.0;
          }
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq > radius[i]*radius[i]) {
        if (history)
          for (j = 0; j < sheardim; j++)
            shearone[i][j] = 0.0;

      } else {

        // meff = effective mass of sphere
        // if I is part of rigid body, use body mass

        meff = rmass[i];
        if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

        // invoke sphere/wall interaction

        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,vwall,v[i],f[i],
                omega[i],torque[i],radius[i],meff);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,vwall,v[i],f[i],
                        omega[i],torque[i],radius[i],meff,shearone[i]);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
                        omega[i],torque[i],radius[i],meff,shearone[i]);
        else if (pairstyle == BONDED_HISTORY)
          bonded_history(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
                         omega[i],torque[i],radius[i],meff,shearone[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke(double rsq, double dx, double dy, double dz,
                        double *vwall, double *v,
                        double *f, double *omega, double *torque,
                        double radius, double meff)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = meff*gammat*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // tangential force due to tangential velocity damping

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke_history(double rsq, double dx, double dy, double dz,
                                double *vwall, double *v,
                                double *f, double *omega, double *torque,
                                double radius, double meff, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (shearupdate) {
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  if (shearupdate) {
    shear[0] -= rsht*dx;
    shear[1] -= rsht*dy;
    shear[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = - (kt*shear[0] + meff*gammat*vtr1);
  fs2 = - (kt*shear[1] + meff*gammat*vtr2);
  fs3 = - (kt*shear[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + meff*gammat*vtr1/kt) -
        meff*gammat*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + meff*gammat*vtr2/kt) -
        meff*gammat*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + meff*gammat*vtr3/kt) -
        meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hertz_history(double rsq, double dx, double dy, double dz,
                                double *vwall, double rwall, double *v,
                                double *f, double *omega, double *torque,
                                double radius, double meff, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,polyhertz,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

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

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hertzian contact + normal velocity damping
  // rwall = 0 is flat wall case
  // rwall positive or negative is curved wall
  //   will break (as it should) if rwall is negative and
  //   its absolute value < radius of particle

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;
  if (rwall == 0.0) polyhertz = sqrt((radius-r)*radius);
  else polyhertz = sqrt((radius-r)*radius*rwall/(rwall+radius));
  ccel *= polyhertz;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (shearupdate) {
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  if (shearupdate) {
    shear[0] -= rsht*dx;
    shear[1] -= rsht*dy;
    shear[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt*shear[0] + meff*gammat*vtr1);
  fs2 = -polyhertz * (kt*shear[1] + meff*gammat*vtr2);
  fs3 = -polyhertz * (kt*shear[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + meff*gammat*vtr1/kt) -
        meff*gammat*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + meff*gammat*vtr2/kt) -
        meff*gammat*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + meff*gammat*vtr3/kt) -
        meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}


/* ---------------------------------------------------------------------- */

void FixWallGran::bonded_history(double rsq, double dx, double dy, double dz,
                                 double *vwall, double rwall, double *v,
                                 double *f, double *omega, double *torque,
                                 double radius, double meff, double *shear)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,polyhertz,rinv,rsqinv;

  double pois,E_eff,G_eff,rad_eff;
  double a0,Fcrit,delcrit,delcritinv;
  double overlap,olapsq,olapcubed,sqrtterm,tmp,keyterm,keyterm2,keyterm3;
  double aovera0,foverFc;
  double gammatsuji;

  double ktwist,kroll,twistcrit,rollcrit;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;
  double magtwist,magtortwist;
  double magrollsq,magroll,magrollinv,magtorroll;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

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

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hertzian contact + normal velocity damping
  // material properties: currently assumes identical materials

  pois = E/(2.0*G) - 1.0;
  E_eff=0.5*E/(1.0-pois*pois);
  G_eff=G/(4.0-2.0*pois);

  // rwall = 0 is infinite wall radius of curvature (flat wall)

  if (rwall == 0) rad_eff = radius;
  else rad_eff = radius*rwall/(radius+rwall);

  Fcrit = rad_eff * (3.0 * M_PI * SurfEnergy);
  a0=pow(9.0*M_PI*SurfEnergy*rad_eff*rad_eff/E_eff,1.0/3.0);
  delcrit = 1.0/rad_eff*(0.5 * a0*a0/pow(6.0,1.0/3.0));
  delcritinv = 1.0/delcrit;

  overlap = (radius-r) * delcritinv;
  olapsq = overlap*overlap;
  olapcubed = olapsq*overlap;
  sqrtterm = sqrt(1.0 + olapcubed);
  tmp = 2.0 + olapcubed + 2.0*sqrtterm;
  keyterm = pow(tmp,THIRD);
  keyterm2 = olapsq/keyterm;
  keyterm3 = sqrt(overlap + keyterm2 + keyterm);
  aovera0 = pow(6.0,-TWOTHIRDS) * (keyterm3 +
            sqrt(2.0*overlap - keyterm2 - keyterm + 4.0/keyterm3));
  foverFc = 4.0*((aovera0*aovera0*aovera0) - pow(aovera0,1.5));
  ccel = Fcrit*foverFc*rinv;

  // damp = meff*gamman*vnnr*rsqinv;
  // ccel = kn*(radius-r)*rinv - damp;
  // polyhertz = sqrt((radius-r)*radius);
  // ccel *= polyhertz;

  // use Tsuji et al form

  polyhertz = 1.2728- 4.2783*0.9 + 11.087*0.9*0.9 - 22.348*0.9*0.9*0.9 +
    27.467*0.9*0.9*0.9*0.9 - 18.022*0.9*0.9*0.9*0.9*0.9 +
    4.8218*0.9*0.9*0.9*0.9*0.9*0.9;

  gammatsuji = 0.2*sqrt(meff*kn);
  damp = gammatsuji*vnnr/rsq;
  ccel = ccel - polyhertz * damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (shearupdate) {
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  if (shearupdate) {
    shear[0] -= rsht*dx;
    shear[1] -= rsht*dy;
    shear[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt*shear[0] + meff*gammat*vtr1);
  fs2 = -polyhertz * (kt*shear[1] + meff*gammat*vtr2);
  fs3 = -polyhertz * (kt*shear[2] + meff*gammat*vtr3);

  kt=8.0*G_eff*a0*aovera0;

  // shear damping uses Tsuji et al form also

  fs1 = -kt*shear[0] - polyhertz*gammatsuji*vtr1;
  fs2 = -kt*shear[1] - polyhertz*gammatsuji*vtr2;
  fs3 = -kt*shear[2] - polyhertz*gammatsuji*vtr3;

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r + 2.0*Fcrit);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + polyhertz*gammatsuji*vtr1/kt) -
      polyhertz*gammatsuji*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + polyhertz*gammatsuji*vtr2/kt) -
      polyhertz*gammatsuji*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + polyhertz*gammatsuji*vtr3/kt) -
      polyhertz*gammatsuji*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // calculate twisting and rolling components of torque
  // NOTE: this assumes spheres!

  relrot1 = omega[0];
  relrot2 = omega[1];
  relrot3 = omega[2];

  // rolling velocity
  // NOTE: this assumes mondisperse spheres!

  vrl1 = -rad_eff*rinv * (relrot2*dz - relrot3*dy);
  vrl2 = -rad_eff*rinv * (relrot3*dx - relrot1*dz);
  vrl3 = -rad_eff*rinv * (relrot1*dy - relrot2*dx);
  vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
  if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
  else vrlmaginv = 0.0;

  // bond history effects

  shear[3] += vrl1*dt;
  shear[4] += vrl2*dt;
  shear[5] += vrl3*dt;

  // rotate bonded displacements correctly

  double rlt = shear[3]*dx + shear[4]*dy + shear[5]*dz;
  rlt /= rsq;
  shear[3] -= rlt*dx;
  shear[4] -= rlt*dy;
  shear[5] -= rlt*dz;

  // twisting torque

  magtwist = rinv*(relrot1*dx + relrot2*dy + relrot3*dz);
  shear[6] += magtwist*dt;

  ktwist = 0.5*kt*(a0*aovera0)*(a0*aovera0);
  magtortwist = -ktwist*shear[6] -
    0.5*polyhertz*gammatsuji*(a0*aovera0)*(a0*aovera0)*magtwist;

  twistcrit=TWOTHIRDS*a0*aovera0*Fcrit;
  if (fabs(magtortwist) > twistcrit)
    magtortwist = -twistcrit * magtwist/fabs(magtwist);

  // rolling torque

  magrollsq = shear[3]*shear[3] + shear[4]*shear[4] + shear[5]*shear[5];
  magroll = sqrt(magrollsq);
  if (magroll != 0.0) magrollinv = 1.0/magroll;
  else magrollinv = 0.0;

  kroll = 1.0*4.0*Fcrit*pow(aovera0,1.5);
  magtorroll = -kroll*magroll - 0.1*gammat*vrlmag;

  rollcrit = 0.01;
  if (magroll > rollcrit) magtorroll = -kroll*rollcrit;

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;

  torque[0] += magtortwist * dx*rinv;
  torque[1] += magtortwist * dy*rinv;
  torque[2] += magtortwist * dz*rinv;

  torque[0] += magtorroll * (shear[4]*dz - shear[5]*dy)*rinv*magrollinv;
  torque[1] += magtorroll * (shear[5]*dx - shear[3]*dz)*rinv*magrollinv;
  torque[2] += magtorroll * (shear[3]*dy - shear[4]*dx)*rinv*magrollinv;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  if (history) bytes += nmax*sheardim * sizeof(double);   // shear history
  if (fix_rigid) bytes += nmax * sizeof(int);             // mass_rigid
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  if (history) memory->grow(shearone,nmax,sheardim,"fix_wall_gran:shearone");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::copy_arrays(int i, int j, int /*delflag*/)
{
  if (history)
    for (int m = 0; m < sheardim; m++)
      shearone[j][m] = shearone[i][m];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGran::set_arrays(int i)
{
  if (history)
    for (int m = 0; m < sheardim; m++)
      shearone[i][m] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGran::pack_exchange(int i, double *buf)
{
  if (!history) return 0;

  int n = 0;
  for (int m = 0; m < sheardim; m++)
    buf[n++] = shearone[i][m];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGran::unpack_exchange(int nlocal, double *buf)
{
  if (!history) return 0;

  int n = 0;
  for (int m = 0; m < sheardim; m++)
    shearone[nlocal][m] = buf[n++];
  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGran::pack_restart(int i, double *buf)
{
  if (!history) return 0;

  int n = 0;
  buf[n++] = sheardim + 1;
  for (int m = 0; m < sheardim; m++)
    buf[n++] = shearone[i][m];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGran::unpack_restart(int nlocal, int nth)
{
  if (!history) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int i = 0; i < sheardim; i++)
    shearone[nlocal][i] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGran::maxsize_restart()
{
  if (!history) return 0;
  return 1 + sheardim;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGran::size_restart(int /*nlocal*/)
{
  if (!history) return 0;
  return 1 + sheardim;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::reset_dt()
{
  dt = update->dt;
}
