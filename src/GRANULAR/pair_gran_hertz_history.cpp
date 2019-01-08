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

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "pair_gran_hertz_history.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHertzHistory::PairGranHertzHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp) {}

/* ---------------------------------------------------------------------- */

void PairGranHertzHistory::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht,polyhertz;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi*omega[i][0] + radj*omega[j][0]) * rinv;
        wr2 = (radi*omega[i][1] + radj*omega[j][1]) * rinv;
        wr3 = (radi*omega[i][2] + radj*omega[j][2]) * rinv;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi*mj / (mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        // normal force = Hertzian contact + normal velocity damping

        damp = meff*gamman*vnnr*rsqinv;
        ccel = kn*(radsum-r)*rinv - damp;
        polyhertz = sqrt((radsum-r)*radi*radj / radsum);
        ccel *= polyhertz;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        shear = &allshear[3*jj];
        if (shearupdate) {
          shear[0] += vtr1*dt;
          shear[1] += vtr2*dt;
          shear[2] += vtr3*dt;
        }
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
                      shear[2]*shear[2]);

        // rotate shear displacements

        rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
        rsht *= rsqinv;
        if (shearupdate) {
          shear[0] -= rsht*delx;
          shear[1] -= rsht*dely;
          shear[2] -= rsht*delz;
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
            fs1 *= fn/fs;
            fs2 *= fn/fs;
            fs3 *= fn/fs;
          } else fs1 = fs2 = fs3 = 0.0;
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);
        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= radj*tor1;
          torque[j][1] -= radj*tor2;
          torque[j][2] -= radj*tor3;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHertzHistory::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  kn = force->numeric(FLERR,arg[0]);
  if (strcmp(arg[1],"NULL") == 0) kt = kn * 2.0/7.0;
  else kt = force->numeric(FLERR,arg[1]);

  gamman = force->numeric(FLERR,arg[2]);
  if (strcmp(arg[3],"NULL") == 0) gammat = 0.5 * gamman;
  else gammat = force->numeric(FLERR,arg[3]);

  xmu = force->numeric(FLERR,arg[4]);
  dampflag = force->inumeric(FLERR,arg[5]);
  if (dampflag == 0) gammat = 0.0;

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
      xmu < 0.0 || xmu > 10000.0 || dampflag < 0 || dampflag > 1)
    error->all(FLERR,"Illegal pair_style command");

  // convert Kn and Kt from pressure units to force/distance^2

  kn /= force->nktv2p;
  kt /= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

double PairGranHertzHistory::single(int i, int j, int /*itype*/, int /*jtype*/,
                                    double rsq,
                                    double /*factor_coul*/, double /*factor_lj*/,
                                    double &fforce)
{
  double radi,radj,radsum;
  double r,rinv,rsqinv,delx,dely,delz;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double mi,mj,meff,damp,ccel,polyhertz;
  double vtr1,vtr2,vtr3,vrel,shrmag,rsht;
  double fs1,fs2,fs3,fs,fn;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  if (rsq >= radsum*radsum) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi*omega[i][0] + radj*omega[j][0]) * rinv;
  wr2 = (radi*omega[i][1] + radj*omega[j][1]) * rinv;
  wr3 = (radi*omega[i][2] + radj*omega[j][2]) * rinv;

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: insure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;


  // normal force = Hertzian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radsum-r)*rinv - damp;
  polyhertz = sqrt((radsum-r)*radi*radj / radsum);
  ccel *= polyhertz;

  // relative velocities

  vtr1 = vt1 - (delz*wr2-dely*wr3);
  vtr2 = vt2 - (delx*wr3-delz*wr1);
  vtr3 = vt3 - (dely*wr1-delx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  double *allshear = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  double *shear = &allshear[3*neighprev];
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
                shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
  rsht *= rsqinv;

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt*shear[0] + meff*gammat*vtr1);
  fs2 = -polyhertz * (kt*shear[1] + meff*gammat*vtr2);
  fs3 = -polyhertz * (kt*shear[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      fs1 *= fn/fs;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
      fs *= fn/fs;
    } else fs1 = fs2 = fs3 = fs = 0.0;
  }

  // set force and return no energy

  fforce = ccel;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;

  return 0.0;
}
