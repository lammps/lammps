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
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"

#include "string.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryOMP::PairGranHookeHistoryOMP(LAMMPS *lmp) :
  PairGranHookeHistory(lmp), ThrOMP(lmp, PAIR)
{
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
    ev_setup_thr(this);
  } else evflag = vflag_fdotr = 0;

  const int shearupdate = (update->ntimestep > laststep) ? 1 : 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f, **torque;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);
    torque = atom->torque + tid*nall;

    if (evflag)
      if (shearupdate) eval<1,1>(f, torque, ifrom, ito, tid);
      else eval<1,0>(f, torque, ifrom, ito, tid);
    else 
      if (shearupdate) eval<0,1>(f, torque, ifrom, ito, tid);
      else eval<0,0>(f, torque, ifrom, ito, tid);

    // reduce per thread forces and torque into global arrays.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
    data_reduce_thr(&(atom->torque[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);

  laststep = update->ntimestep;
}

template <int EVFLAG, int SHEARUPDATE>
void PairGranHookeHistoryOMP::eval(double **f, double **torque, int iifrom, int iito, int tid)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  double **x = atom->x;
  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
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

	// normal forces = Hookian contact + normal velocity damping

	if (rmass) {
	  meff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	  if (mask[i] & freeze_group_bit) meff = rmass[j];
	  if (mask[j] & freeze_group_bit) meff = rmass[i];
	} else {
	  itype = type[i];
	  jtype = type[j];
	  meff = mass[itype]*mass[jtype] / (mass[itype]+mass[jtype]);
	  if (mask[i] & freeze_group_bit) meff = mass[jtype];
	  if (mask[j] & freeze_group_bit) meff = mass[itype];
	}

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
	shear = &allshear[3*jj];

	if (SHEARUPDATE) {
	  shear[0] += vtr1*dt;
	  shear[1] += vtr2*dt;
	  shear[2] += vtr3*dt;
	}
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
		      shear[2]*shear[2]);

	// rotate shear displacements

	rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
	rsht *= rsqinv;
	if (SHEARUPDATE) {
	  shear[0] -= rsht*delx;
	  shear[1] -= rsht*dely;
	  shear[2] -= rsht*delz;
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
	    const double fnfs = fn/fs;
	    const double mgkt = meff*gammat/kt;
	    shear[0] = fnfs * (shear[0] + mgkt*vtr1) - mgkt*vtr1;
	    shear[1] = fnfs * (shear[1] + mgkt*vtr2) - mgkt*vtr2;
	    shear[2] = fnfs * (shear[2] + mgkt*vtr3) - mgkt*vtr3;
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
				     0.0,0.0,fx,fy,fz,delx,dely,delz,tid);

      }
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
