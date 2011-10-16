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
#include "pair_lubricate_omp.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "neighbor.h"
#include "random_mars.h"
#include "neigh_list.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLubricateOMP::PairLubricateOMP(LAMMPS *lmp) :
  PairLubricate(lmp), ThrOMP(lmp, THR_PAIR)
{
  respa_enable = 0;
  random_thr = NULL;
}

/* ---------------------------------------------------------------------- */

PairLubricateOMP::~PairLubricateOMP()
{
  if (random_thr) {
    for (int i=1; i < comm->nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void PairLubricateOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  if (!random_thr)
    random_thr = new RanMars*[nthreads];
  
  random_thr[0] = random;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (random_thr && tid > 0)
      random_thr[tid] = new RanMars(Pair::lmp, seed + comm->me 
				    + comm->nprocs*tid);

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
	else eval<1,1,0>(ifrom, ito, thr);
      } else {
	if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
	else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    reduce_thr(eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLubricateOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,w1,w2,w3,v_shear1,v_shear2,v_shear3;
  double omega_t_1,omega_t_2,omega_t_3;
  double n_cross_omega_t_1,n_cross_omega_t_2,n_cross_omega_t_3;
  double wr1,wr2,wr3,wnnr,wn1,wn2,wn3;
  double P_dot_wrel_1,P_dot_wrel_2,P_dot_wrel_3;
  double a_squeeze,a_shear,a_pump,a_twist;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * const * const x = atom->x;
  const double * const * const v = atom->v;
  const double * const * const omega = atom->omega;
  const double * const radius = atom->radius;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double vxmu2f = force->vxmu2f;
  RanMars &rng = *random_thr[thr->get_tid()];

  double prethermostat = sqrt(2.0 * force->boltz * t_target / update->dt);
  prethermostat *= sqrt(force->vxmu2f/force->ftm2v/force->mvv2e);

  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  a_squeeze = a_shear = a_pump = a_twist = 0.0;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    itype = type[i];
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
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

	r = sqrt(rsq);

        // relative translational velocity 

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component N.(v1-v2) = nn.(v1-v2)

        vnnr = vr1*delx + vr2*dely + vr3*delz;
	vnnr /= r;
	vn1 = delx*vnnr / r;
        vn2 = dely*vnnr / r;
        vn3 = delz*vnnr / r;

        // tangential component -P.(v1-v2)
	// P = (I - nn) where n is vector between centers
     
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // additive rotational velocity = omega_1 + omega_2

	w1 = omega[i][0] + omega[j][0];
	w2 = omega[i][1] + omega[j][1];
	w3 = omega[i][2] + omega[j][2];

        // relative velocities n X P . (v1-v2) = n X (I-nn) . (v1-v2)

        v_shear1 = (dely*vt3 - delz*vt2) / r;
        v_shear2 = -(delx*vt3 - delz*vt1) / r;
        v_shear3 = (delx*vt2 - dely*vt1) / r;

        // relative rotation rate P.(omega1 + omega2)

	omega_t_1 = w1 - delx*(delx*w1) / rsq;
	omega_t_2 = w2 - dely*(dely*w2) / rsq;
	omega_t_3 = w3 - delz*(delz*w3) / rsq;

        // n X omega_t

        n_cross_omega_t_1 =  (dely*omega_t_3 - delz*omega_t_2) / r;
        n_cross_omega_t_2 =  -(delx*omega_t_3 - delz*omega_t_1) / r;
        n_cross_omega_t_3 =  (delx*omega_t_2 - dely*omega_t_1) / r;

        // N.(w1-w2) and P.(w1-w2)

	wr1 = omega[i][0] - omega[j][0];
	wr2 = omega[i][1] - omega[j][1];
	wr3 = omega[i][2] - omega[j][2];
 
	wnnr = wr1*delx + wr2*dely + wr3*delz;
	wn1 = delx*wnnr / rsq;
	wn2 = dely*wnnr / rsq;
	wn3 = delz*wnnr / rsq;

        P_dot_wrel_1 = wr1 - delx*(delx*wr1)/rsq; 
        P_dot_wrel_2 = wr2 - dely*(dely*wr2)/rsq; 
        P_dot_wrel_3 = wr3 - delz*(delz*wr3)/rsq; 

        // compute components of pair-hydro

        h_sep = r - 2.0*radi;

	if (flag1)
	  a_squeeze = (3.0*MY_PI*mu*2.0*radi/2.0) * (2.0*radi/4.0/h_sep);
	if (flag2) 
	  a_shear = (MY_PI*mu*2.*radi/2.0) *
	    log(2.0*radi/2.0/h_sep)*(2.0*radi+h_sep)*(2.0*radi+h_sep)/4.0;
	if (flag3) 
	  a_pump = (MY_PI*mu*pow(2.0*radi,4)/8.0) *
	    ((3.0/20.0) * log(2.0*radi/2.0/h_sep) + 
	     (63.0/250.0) * (h_sep/2.0/radi) * log(2.0*radi/2.0/h_sep));
	if (flag4)
	  a_twist = (MY_PI*mu*pow(2.0*radi,4)/4.0) *
	    (h_sep/2.0/radi) * log(2.0/(2.0*h_sep));

        if (h_sep >= cut_inner[itype][jtype]) {
          fx = -a_squeeze*vn1 - a_shear*(2.0/r)*(2.0/r)*vt1 + 
	    (2.0/r)*a_shear*n_cross_omega_t_1;
          fy = -a_squeeze*vn2 - a_shear*(2.0/r)*(2.0/r)*vt2 + 
	    (2.0/r)*a_shear*n_cross_omega_t_2;
          fz = -a_squeeze*vn3 - a_shear*(2.0/r)*(2.0/r)*vt3 +
	    (2.0/r)*a_shear*n_cross_omega_t_3;
	  fx *= vxmu2f;
	  fy *= vxmu2f;
	  fz *= vxmu2f;

	  // add in thermostat force

	  tfmag = prethermostat*sqrt(a_squeeze)*(rng.uniform()-0.5);
	  fx -= tfmag * delx/r;
	  fy -= tfmag * dely/r;
	  fz -= tfmag * delz/r;
	  
	  tx = -(2.0/r)*a_shear*v_shear1 - a_shear*omega_t_1 - 
	    a_pump*P_dot_wrel_1 - a_twist*wn1;
	  ty = -(2.0/r)*a_shear*v_shear2 - a_shear*omega_t_2 - 
	    a_pump*P_dot_wrel_2 - a_twist*wn2;
	  tz = -(2.0/r)*a_shear*v_shear3 - a_shear*omega_t_3 - 
	    a_pump*P_dot_wrel_3 - a_twist*wn3;
	  t1tmp += vxmu2f * tx;
	  t2tmp += vxmu2f * ty;
	  t3tmp += vxmu2f * tz;

        } else {
	  a_squeeze = (3.0*MY_PI*mu*2.0*radi/2.0) * 
	    (2.0*radi/4.0/cut_inner[itype][jtype]);
	  fpair = -a_squeeze*vnnr;
	  fpair *= vxmu2f;

	  // add in thermostat force

	  fpair -= prethermostat*sqrt(a_squeeze)*(rng.uniform()-0.5);

	  fx = fpair * delx/r;
	  fy = fpair * dely/r;
	  fz = fpair * delz/r;
	}

    	fxtmp += fx;
	fytmp += fy;
	fztmp += fz;

	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;

	  if (h_sep >= cut_inner[itype][jtype]) {
	    tx = -(2.0/r)*a_shear*v_shear1 - a_shear*omega_t_1 + 
	      a_pump*P_dot_wrel_1 + a_twist*wn1;
	    ty = -(2.0/r)*a_shear*v_shear2 - a_shear*omega_t_2 + 
	      a_pump*P_dot_wrel_2 + a_twist*wn2;
	    tz = -(2.0/r)*a_shear*v_shear3 - a_shear*omega_t_3 + 
	      a_pump*P_dot_wrel_3 + a_twist*wn3;
	    torque[j][0] += vxmu2f * tx;
	    torque[j][1] += vxmu2f * ty;
	    torque[j][2] += vxmu2f * tz;
	  }
	}

	if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
				     0.0,0.0,fx,fy,fz,delx,dely,delz,thr);
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

double PairLubricateOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLubricate::memory_usage();

  return bytes;
}
