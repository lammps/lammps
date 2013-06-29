/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_gran_hooke_history_omp.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"

#include "string.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryOMP::PairGranHookeHistoryOMP(LAMMPS *lmp) :
  PairGranHookeHistory(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  // trigger use of OpenMP version of FixShearHistory
  suffix = new char[4];
  memcpy(suffix,"omp",4);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  computeflag = 1;
  const int shearupdate = (update->setupflag) ? 0 : 1;

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
    for (int i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag)
      if (shearupdate) eval<1,1>(ifrom, ito, thr);
      else eval<1,0>(ifrom, ito, thr);
    else
      if (shearupdate) eval<0,1>(ifrom, ito, thr);
      else eval<0,0>(ifrom, ito, thr);

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int SHEARUPDATE>
void PairGranHookeHistoryOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double myshear[3];
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *allshear,**firstshear;

  const double * const * const x = atom->x;
  const double * const * const v = atom->v;
  const double * const * const omega = atom->omega;
  const double * const radius = atom->radius;
  const double * const rmass = atom->rmass;
  const double * const mass = atom->mass;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const int * const type = atom->type;
  const int * const mask = atom->mask;
  const int nlocal = atom->nlocal;
  double fxtmp,fytmp,fztmp;
  double t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=t1tmp=t2tmp=t3tmp=0.0;

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
        myshear[0] = 0.0;
        myshear[1] = 0.0;
        myshear[2] = 0.0;

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

        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          mi = mass[type[i]];
          mj = mass[type[j]];
        }
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

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        memcpy(myshear,allshear + 3*jj, 3*sizeof(double));

        if (SHEARUPDATE) {
          myshear[0] += vtr1*dt;
          myshear[1] += vtr2*dt;
          myshear[2] += vtr3*dt;
        }
        shrmag = sqrt(myshear[0]*myshear[0] + myshear[1]*myshear[1] +
                      myshear[2]*myshear[2]);

        // rotate shear displacements

        rsht = myshear[0]*delx + myshear[1]*dely + myshear[2]*delz;
        rsht *= rsqinv;
        if (SHEARUPDATE) {
          myshear[0] -= rsht*delx;
          myshear[1] -= rsht*dely;
          myshear[2] -= rsht*delz;
        }

        // tangential forces = shear + tangential velocity damping

        fs1 = - (kt*myshear[0] + meff*gammat*vtr1);
        fs2 = - (kt*myshear[1] + meff*gammat*vtr2);
        fs3 = - (kt*myshear[2] + meff*gammat*vtr3);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        if (fs > fn) {
          if (shrmag != 0.0) {
            const double fnfs = fn/fs;
            const double mgkt = meff*gammat/kt;
            myshear[0] = fnfs * (myshear[0] + mgkt*vtr1) - mgkt*vtr1;
            myshear[1] = fnfs * (myshear[1] + mgkt*vtr2) - mgkt*vtr2;
            myshear[2] = fnfs * (myshear[2] + mgkt*vtr3) - mgkt*vtr3;
            fs1 *= fnfs;
            fs2 *= fnfs;
            fs3 *= fnfs;
          } else fs1 = fs2 = fs3 = 0.0;
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        fxtmp  += fx;
        fytmp  += fy;
        fztmp  += fz;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);
        t1tmp -= radi*tor1;
        t2tmp -= radi*tor2;
        t3tmp -= radi*tor3;

        if (j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= radj*tor1;
          torque[j][1] -= radj*tor2;
          torque[j][2] -= radj*tor3;
        }

        if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,/* newton_pair */ 0,
                                     0.0,0.0,fx,fy,fz,delx,dely,delz,thr);

      }
      memcpy(allshear + 3*jj, myshear, 3*sizeof(double));
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
    torque[i][0] += t1tmp;
    torque[i][1] += t2tmp;
    torque[i][2] += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairGranHookeHistoryOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairGranHookeHistory::memory_usage();

  return bytes;
}
