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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "pair_gran_hooke.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "fix.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHooke::PairGranHooke(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
  no_virial_fdotr_compute = 0;
  history = 0;
}

/* ---------------------------------------------------------------------- */

void PairGranHooke::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,factor_lj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,ft,fs1,fs2,fs3;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
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
  double *special_lj = force->special_lj;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq < radsum*radsum) {
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

        // normal forces = Hookian contact + normal velocity damping

        damp = meff*gamman*vnnr*rsqinv;
        ccel = kn*(radsum-r)*rinv - damp;
        if (limit_damping && (ccel < 0.0)) ccel = 0.0;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
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

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        fx *= factor_lj;
        fy *= factor_lj;
        fz *= factor_lj;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);
        tor1 *= factor_lj;
        tor2 *= factor_lj;
        tor3 *= factor_lj;
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

/* ---------------------------------------------------------------------- */

double PairGranHooke::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                             double /*factor_coul*/, double /*factor_lj*/,
                             double &fforce)
{
  double radi,radj,radsum,r,rinv,rsqinv;
  double delx,dely,delz;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel;
  double fn,fs,ft;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  // zero out forces if caller requests non-touching pair outside cutoff

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

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radsum-r)*rinv - damp;
  if (limit_damping && (ccel < 0.0)) ccel = 0.0;

  // relative velocities

  vtr1 = vt1 - (delz*wr2-dely*wr3);
  vtr2 = vt2 - (delx*wr3-delz*wr1);
  vtr3 = vt3 - (dely*wr1-delx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = meff*gammat*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // set force and return no energy

  fforce = ccel;


  // set single_extra quantities

  svector[0] = -ft*vtr1;
  svector[1] = -ft*vtr2;
  svector[2] = -ft*vtr3;
  svector[3] = sqrt(svector[0]*svector[0] +
                    svector[1]*svector[1] +
                    svector[2]*svector[2]);
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;

  return 0.0;
}
