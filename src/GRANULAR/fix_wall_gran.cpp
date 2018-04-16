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

#include <math.h>
#include <stdlib.h>
#include <string.h>
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
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,JKR_ROLLING,DMT_ROLLING};
enum{NONE,CONSTANT,EQUAL};

enum {TSUJI, BRILLIANTOV};
enum {INDEP, BRILLROLL};

#define BIG 1.0e20
#define EPSILON 1e-10

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
  else if (strcmp(arg[3],"dmt/rolling") == 0) pairstyle = DMT_ROLLING;
  //else if (strcmp(arg[3],"jkr/rolling") == 0) pairstyle = JKR_ROLLING;
  else error->all(FLERR,"Invalid fix wall/gran interaction style");

  history = restart_peratom = 1;
  if (pairstyle == HOOKE) history = restart_peratom = 0;

  // wall/particle coefficients

  int iarg;

  if (pairstyle != JKR_ROLLING && pairstyle != DMT_ROLLING) {
    sheardim = 3;
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
    if (narg < 12) error->all(FLERR,"Illegal fix wall/gran command");

    sheardim = 7;
    Emod = force->numeric(FLERR,arg[4]);
    Gmod = force->numeric(FLERR,arg[5]);
    xmu = force->numeric(FLERR,arg[6]);
    gamman = force->numeric(FLERR,arg[7]);
    Ecoh = force->numeric(FLERR,arg[8]);
    kR = force->numeric(FLERR,arg[9]);
    muR = force->numeric(FLERR,arg[10]);
    etaR = force->numeric(FLERR,arg[11]);

    //Defaults
    normaldamp = TSUJI;
    rollingdamp = INDEP;

    iarg = 12;
    for (int iiarg=iarg; iiarg < narg; ++iiarg){
      if (strcmp(arg[iiarg], "normaldamp") == 0){
        if(iiarg+2 > narg) error->all(FLERR, "Invalid fix/wall/gran region command");
        if (strcmp(arg[iiarg+1],"tsuji") == 0){
          normaldamp = TSUJI;
          alpha = gamman;
        }
        else if (strcmp(arg[iiarg+1],"brilliantov") == 0) normaldamp = BRILLIANTOV;
        else error->all(FLERR, "Invalid normal damping model for fix wall/gran dmt/rolling");
        iarg += 2;
      }
      if (strcmp(arg[iiarg], "rollingdamp") == 0){
        if(iiarg+2 > narg) error->all(FLERR, "Invalid fix/wall/gran region command");
        if (strcmp(arg[iarg+1],"independent") == 0) rollingdamp = INDEP;
        else if (strcmp(arg[iarg+1],"brilliantov") == 0) rollingdamp = BRILLROLL;
        else error->all(FLERR, "Invalid rolling damping model for fix wall/gran dmt/rolling");
        iarg += 2;
      }
    }
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
  peratom_flag = 0;

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
    } else if (strcmp(arg[iarg],"store_contacts") == 0){
      peratom_flag = 1;
      size_peratom_cols = 8; //Could make this a user input option?
      peratom_freq = 1;
      iarg += 1;
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

  if (peratom_flag){
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < size_peratom_cols; m++)
        array_atom[i][m] = 0.0;
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
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

void FixWallGran::post_force(int vflag)
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

        // store contact info
        if (peratom_flag){
          array_atom[i][0] = (double)atom->tag[i];
          array_atom[i][4] = x[i][0] - dx;
          array_atom[i][5] = x[i][1] - dy;
          array_atom[i][6] = x[i][2] - dz;
          array_atom[i][7] = radius[i];
        }

        // invoke sphere/wall interaction
        double *contact;
        if (peratom_flag)
          contact = array_atom[i];
        else
          contact = NULL;

        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,vwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff, contact);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,vwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,shearone[i],
              contact);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,shearone[i],
              contact);
        else if (pairstyle == DMT_ROLLING)
          dmt_rolling(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,shearone[i],
              contact);
        /*else if (pairstyle == JKR_ROLLING)
          jkr_rolling(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,shearone[i],
              contact);*/
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke(double rsq, double dx, double dy, double dz,
    double *vwall, double *v,
    double *f, double *omega, double *torque,
    double radius, double meff, double* contact)
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

  if (peratom_flag){
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

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
    double radius, double meff, double *shear,
    double *contact)
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

  if (peratom_flag){
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

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
    double radius, double meff, double *shear,
    double *contact)
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

  if (peratom_flag){
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

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


void FixWallGran::dmt_rolling(double rsq, double dx, double dy, double dz,
    double *vwall, double rwall, double *v,
    double *f, double *omega, double *torque,
    double radius, double meff, double *shear,
    double *contact)
{
  int i,j,ii,jj,inum,jnum;
  int itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,nx,ny,nz;
  double radi,radj,radsum,r,rinv,rsqinv,R,a;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double kn, kt, k_Q, k_R, eta_N, eta_T, eta_Q, eta_R;
  double Fhz, Fdamp, Fdmt, Fne, Fntot, Fscrit, Frcrit;
  double overlap;
  double mi,mj,damp,ccel,tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double fs,fs1,fs2,fs3,roll1,roll2,roll3,torroll1,torroll2,torroll3;
  double tortwist1, tortwist2, tortwist3;
  double shrmag,rsht;


  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  radsum = radius + rwall;
  if (rwall == 0) R = radius;
  else R = radius*rwall/(radius+rwall);

  nx = delx*rinv;
  ny = dely*rinv;
  nz = delz*rinv;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*nx + vr2*ny + vr3*nz; //v_R . n
  vn1 = nx*vnnr;
  vn2 = ny*vnnr;
  vn3 = nz*vnnr;


  //****************************************
  //Normal force = Hertzian contact + DMT + damping
  //****************************************
  overlap = radsum - r;
  a = sqrt(R*overlap);
  kn = 4.0/3.0*Emod*a;
  Fhz = kn*overlap;

  //Damping (based on Tsuji et al)
  if (normaldamp == BRILLIANTOV) eta_N = a*meff*gamman;
  else if (normaldamp == TSUJI) eta_N=alpha*sqrt(meff*kn);

  Fdamp = -eta_N*vnnr; //F_nd eq 23 and Zhao eq 19

  //DMT
  Fdmt = -4*MY_PI*Ecoh*R;

  Fne = Fhz + Fdmt;
  Fntot = Fne + Fdamp;

  //****************************************
  //Tangential force, including shear history effects
  //****************************************

  // tangential component
  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // relative tangential velocities
  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
      shear[2]*shear[2]);

  // Rotate and update shear displacements.
  // See e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (shearupdate) {
    rsht = shear[0]*nx + shear[1]*ny + shear[2]*nz;
    if (fabs(rsht) < EPSILON) rsht = 0;
    if (rsht > 0){
      scalefac = shrmag/(shrmag - rsht); //if rhst == shrmag, contacting pair has rotated 90 deg. in one step, in which case you deserve a crash!
      shear[0] -= rsht*nx;
      shear[1] -= rsht*ny;
      shear[2] -= rsht*nz;
      //Also rescale to preserve magnitude
      shear[0] *= scalefac;
      shear[1] *= scalefac;
      shear[2] *= scalefac;
    }
    //Update shear history
    shear[0] += vtr1*dt;
    shear[1] += vtr2*dt;
    shear[2] += vtr3*dt;
  }

  // tangential forces = shear + tangential velocity damping
  // following Zhao and Marshall Phys Fluids v20, p043302 (2008)
  kt=8.0*Gmod*a;

  eta_T = eta_N; //Based on discussion in Marshall; eta_T can also be an independent parameter
  fs1 = -kt*shear[0] - eta_T*vtr1; //eq 26
  fs2 = -kt*shear[1] - eta_T*vtr2;
  fs3 = -kt*shear[2] - eta_T*vtr3;

  // rescale frictional displacements and forces if needed
  Fscrit = xmu * fabs(Fne);
  // For JKR, use eq 43 of Marshall. For DMT, use Fne instead

  //Redundant, should be same as above?
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
      shear[2]*shear[2]);
  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  if (fs > Fscrit) {
    if (shrmag != 0.0) {
      //shear[0] = (Fcrit/fs) * (shear[0] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
      //shear[1] = (Fcrit/fs) * (shear[1] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
      //shear[2] = (Fcrit/fs) * (shear[2] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
      shear[0] = -1.0/kt*(Fscrit*fs1/fs + eta_T*vtr1); //Same as above, but simpler (check!)
      shear[1] = -1.0/kt*(Fscrit*fs2/fs + eta_T*vtr2);
      shear[2] = -1.0/kt*(Fscrit*fs3/fs + eta_T*vtr3);
      fs1 *= Fscrit/fs;
      fs2 *= Fscrit/fs;
      fs3 *= Fscrit/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  //****************************************
  // Rolling force, including shear history effects
  //****************************************

  relrot1 = omega[0]; //- omega[j][0]; TODO: figure out how to
  relrot2 = omega[1]; //- omega[j][1];   incorporate wall angular
  relrot3 = omega[2]; //- omega[j][2];   velocity

  // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
  // This is different from the Marshall papers, which use the Bagi/Kuhn formulation
  // for rolling velocity (see Wang et al for why the latter is wrong)
  vrl1 = R*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
  vrl2 = R*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
  vrl3 = R*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;
  vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
  if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
  else vrlmaginv = 0.0;

  // Rolling displacement
  rollmag = sqrt(shear[3]*shear[3] + shear[4]*shear[4] + shear[5]*shear[5]);
  rolldotn = shear[3]*nx + shear[4]*ny + shear[5]*nz;

  if (shearupdate) {
    if (fabs(rolldotn) < EPSILON) rolldotn = 0;
    if (rolldotn > 0){ //Rotate into tangential plane
      scalefac = rollmag/(rollmag - rolldotn);
      shear[3] -= rolldotn*nx;
      shear[4] -= rolldotn*ny;
      shear[5] -= rolldotn*nz;
      //Also rescale to preserve magnitude
      shear[3] *= scalefac;
      shear[4] *= scalefac;
      shear[5] *= scalefac;
    }
    shear[3] += vrl1*dt;
    shear[4] += vrl2*dt;
    shear[5] += vrl3*dt;
  }

  if (rollingdamp == BRILLROLL) etaR = muR*fabs(Fne);
  fr1 = -kR*shear[3] - etaR*vrl1;
  fr2 = -kR*shear[4] - etaR*vrl2;
  fr3 = -kR*shear[5] - etaR*vrl3;

  // rescale frictional displacements and forces if needed
  Frcrit = muR * fabs(Fne);

  fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
  if (fr > Frcrit) {
    if (rollmag != 0.0) {
      shear[3] = -1.0/kR*(Frcrit*fr1/fr + etaR*vrl1);
      shear[4] = -1.0/kR*(Frcrit*fr2/fr + etaR*vrl2);
      shear[5] = -1.0/kR*(Frcrit*fr3/fr + etaR*vrl3);
      fr1 *= Frcrit/fr;
      fr2 *= Frcrit/fr;
      fr3 *= Frcrit/fr;
    } else fr1 = fr2 = fr3 = 0.0;
  }


  //****************************************
  // Twisting torque, including shear history effects
  //****************************************
  magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
  shear[6] += magtwist*dt;
  k_Q = 0.5*kt*a*a;; //eq 32
  eta_Q = 0.5*eta_T*a*a;
  magtortwist = -k_Q*shear[6] - eta_Q*magtwist;//M_t torque (eq 30)

  signtwist = (magtwist > 0) - (magtwist < 0);
  Mtcrit=TWOTHIRDS*a*Fscrit;//critical torque (eq 44)
  if (fabs(magtortwist) > Mtcrit){
    shear[6] = 1.0/k_Q*(Mtcrit*signtwist - eta_Q*magtwist);
    magtortwist = -Mtcrit * signtwist; //eq 34
  }

  // Apply forces & torques

  fx = nx*Fntot + fs1;
  fy = ny*Fntot + fs2;
  fz = nz*Fntot + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = ny*fs3 - nz*fs2;
  tor2 = nz*fs1 - nx*fs3;
  tor3 = nx*fs2 - ny*fs1;

  torque[0] -= radi*tor1;
  torque[1] -= radi*tor2;
  torque[2] -= radi*tor3;

  tortwist1 = magtortwist * nx;
  tortwist2 = magtortwist * ny;
  tortwist3 = magtortwist * nz;

  torque[0] += tortwist1;
  torque[1] += tortwist2;
  torque[2] += tortwist3;

  torroll1 = R*(ny*fr3 - nz*fr2); //n cross fr
  torroll2 = R*(nz*fr1 - nx*fr3);
  torroll3 = R*(nx*fr2 - ny*fr1);

  torque[0] += torroll1;
  torque[1] += torroll2;
  torque[2] += torroll3;

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
  if (peratom_flag) bytes += nmax*size_peratom_cols*sizeof(double); //store contacts
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  if (history) memory->grow(shearone,nmax,sheardim,"fix_wall_gran:shearone");
  if (peratom_flag){
    memory->grow(array_atom,nmax,size_peratom_cols,"fix_wall_gran:array_atom");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::copy_arrays(int i, int j, int delflag)
{
  if (history)
    for (int m = 0; m < sheardim; m++)
      shearone[j][m] = shearone[i][m];
  if (peratom_flag){
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[j][m] = array_atom[i][m];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGran::set_arrays(int i)
{
  if (history)
    for (int m = 0; m < sheardim; m++)
      shearone[i][m] = 0;
  if (peratom_flag){
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[i][m] = 0;
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGran::pack_exchange(int i, double *buf)
{
  int n = 0;
  if (history){
    for (int m = 0; m < sheardim; m++)
      buf[n++] = shearone[i][m];
  }
  if (peratom_flag){
    for (int m = 0; m < size_peratom_cols; m++)
      buf[n++] = array_atom[i][m];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGran::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  if (history){
    for (int m = 0; m < sheardim; m++)
      shearone[nlocal][m] = buf[n++];
  }
  if (peratom_flag){
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[nlocal][m] = buf[n++];
  }
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

int FixWallGran::size_restart(int nlocal)
{
  if (!history) return 0;
  return 1 + sheardim;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::reset_dt()
{
  dt = update->dt;
}

