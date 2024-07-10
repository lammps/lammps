/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_line_gran_hooke.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLineGranHooke::PairLineGranHooke(LAMMPS *lmp) :
  PairLineGranHookeHistory(lmp)
{
  history = 0;
}

/* ---------------------------------------------------------------------- */

void PairLineGranHooke::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,inum,jnum,jflag,kflag,otherflag;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double factor_couple;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,mline;
  double damp,ccel,tor1,tor2,tor3;
  double fn,fs,ft,fs1,fs2,fs3;
  double dr[3],ds[3],vline[3],contact[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // if just reneighbored:
  // update rigid body masses for owned/ghost atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body
  // forward comm mass_rigid so have it for ghost lines

  if (neighbor->ago == 0) {
    if (fix_rigid) {
      int tmp;
      int *body = (int *) fix_rigid->extract("body",tmp);
      double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
      if (atom->nmax > nmax) {
        memory->destroy(mass_rigid);
        nmax = atom->nmax;
        memory->create(mass_rigid,nmax,"line/gran:mass_rigid");
      }
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++) {
        if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
        else mass_rigid[i] = 0.0;
      }
      comm->forward_comm(this);
    }

    connect2d = fsl->connect2d;
  }

  // pre-calculate current end pts of owned+ghost lines
  // only once per reneighbor if surfs not moving

  if (surfmoveflag || neighbor->ago == 0) calculate_endpts();

  // loop over neighbors of my atoms
  // I is always sphere, J is always line

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *line = atom->line;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq < radsum*radsum) {

        // sanity check that neighbor list is built correctly

        if (line[i] >= 0 || line[j] < 0)
          error->one(FLERR,"Pair line/gran iteraction is invalid");

        // check for overlap of sphere and line segment
        // jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
        // if no overlap, just continue
        // for overlap, also return:
        //   contact = nearest point on line to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = overlap_sphere_line(i,j,contact,dr,rsq);

        if (!jflag) continue;

        // if contact = line end pt:
        // check overlap status of line adjacent to the end pt
        // otherflag = 0/1 for this/other line performs calculation

        if (jflag < 0) {
          otherflag = endpt_neigh_check(i,j,jflag);
          if (otherflag) continue;
        }

        //printf("CONTACT %ld: %d %d: %d %d\n",update->ntimestep,i,j,
        //      tag[i],tag[j]);

        // NOTE: add logic to check for coupled contacts and weight them

        factor_couple = 1.0;

        // ds = vector from line center to contact pt

        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        ds[0] = contact[0] - x[j][0];
        ds[1] = contact[1] - x[j][1];
        ds[2] = contact[2] - x[j][2];

        // vline = velocity of contact pt on line, translation + rotation

        vline[0] = v[j][0] + (omega[j][1]*ds[2] - omega[j][2]*ds[1]);
        vline[1] = v[j][1] + (omega[j][2]*ds[0] - omega[j][0]*ds[2]);
        vline[2] = v[j][2] + (omega[j][0]*ds[1] - omega[j][1]*ds[0]);

        // relative translational velocity

        vr1 = v[i][0] - vline[0];
        vr2 = v[i][1] - vline[1];
        vr3 = v[i][2] - vline[2];

        // normal component

        vnnr = vr1*dr[0] + vr2*dr[1] + vr3*dr[2];
        vn1 = dr[0]*vnnr * rsqinv;
        vn2 = dr[1]*vnnr * rsqinv;
        vn3 = dr[2]*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi*omega[i][0]) * rinv;
        wr2 = (radi*omega[i][1]) * rinv;
        wr3 = (radi*omega[i][2]) * rinv;

        // meff = effective mass of sphere and line
        // if I or J is part of rigid body, use body mass
        // if line is not part of rigid body assume infinite mass

        meff = rmass[i];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) meff = mass_rigid[i];
          if (mass_rigid[j] > 0.0) {
            mline = mass_rigid[j];
            meff = meff*mline / (meff+mline);
          }
        }

        // normal forces = Hookian contact + normal velocity damping

        damp = meff*gamman*vnnr*rsqinv;
        ccel = kn*(radi-r)*rinv - damp;
        ccel *= factor_couple;

        // relative velocities

        vtr1 = vt1 - (dr[2]*wr2-dr[1]*wr3);
        vtr2 = vt2 - (dr[0]*wr3-dr[2]*wr1);
        vtr3 = vt3 - (dr[1]*wr1-dr[0]*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // force normalization

        fn = xmu * fabs(ccel*r);
        fs = meff*gammat*vrel;
        if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
        else ft = 0.0;

        // tangential force due to tangential velocity damping

        fs1 = -ft*vtr1 * factor_couple;
        fs2 = -ft*vtr2 * factor_couple;
        fs3 = -ft*vtr3 * factor_couple;

        // total force on sphere

        fx = dr[0]*ccel + fs1;
        fy = dr[1]*ccel + fs2;
        fz = dr[2]*ccel + fs3;

        // sphere force & torque

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dr[1]*fs3 - dr[2]*fs2);
        tor2 = rinv * (dr[2]*fs1 - dr[0]*fs3);
        tor3 = rinv * (dr[0]*fs2 - dr[1]*fs1);
        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        // line force & torque
        // torque applied at contact pt
        // use total force for torque
        //   since opposite force is not norm/tang to line at its end pt

        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        torque[j][0] -= ds[1]*fz - ds[2]*fy;
        torque[j][1] -= ds[2]*fx - ds[0]*fz;
        torque[j][2] -= ds[0]*fy - ds[1]*fx;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,dr[0],dr[1],dr[2]);
      }
    }
  }

  // NOTE: should there be virial contributions from boundary tris?

  if (vflag_fdotr) virial_fdotr_compute();
}
