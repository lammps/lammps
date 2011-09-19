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
#include "pair_gayberne_omp.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGayBerneOMP::PairGayBerneOMP(LAMMPS *lmp) :
  PairGayBerne(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairGayBerneOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(this);
  } else evflag = vflag_fdotr = 0;

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

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, torque, ifrom, ito, tid);
	else eval<1,1,0>(f, torque, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, torque, ifrom, ito, tid);
	else eval<1,0,0>(f, torque, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, torque, ifrom, ito, tid);
      else eval<0,0,0>(f, torque, ifrom, ito, tid);
    }

    // reduce per thread forces and torques into global arrays.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
    data_reduce_thr(&(atom->torque[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairGayBerneOMP::eval(double **f, double **tor, int iifrom, int iito, int tid)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  double a1[3][3],b1[3][3],g1[3][3],a2[3][3],b2[3][3],g2[3][3],temp[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *iquat,*jquat;

  double **x = atom->x;
  int *ellipsoid = atom->ellipsoid;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itype = type[i];

    if (form[itype][itype] == ELLIPSE_ELLIPSE) {
      iquat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat_trans(iquat,a1);
      MathExtra::diag_times3(well[itype],a1,temp);
      MathExtra::transpose_times3(a1,temp,b1);
      MathExtra::diag_times3(shape2[itype],a1,temp);
      MathExtra::transpose_times3(a1,temp,g1);
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

	switch (form[itype][jtype]) {
	case SPHERE_SPHERE:
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  forcelj *= -r2inv;
	  if (EFLAG) 
	    one_eng = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
	      offset[itype][jtype];
	  fforce[0] = r12[0]*forcelj;
	  fforce[1] = r12[1]*forcelj;
	  fforce[2] = r12[2]*forcelj;
	  ttor[0] = ttor[1] = ttor[2] = 0.0;
	  rtor[0] = rtor[1] = rtor[2] = 0.0;
	  break;

        case SPHERE_ELLIPSE:
	  jquat = bonus[ellipsoid[j]].quat;
	  MathExtra::quat_to_mat_trans(jquat,a2);
	  MathExtra::diag_times3(well[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,b2);
	  MathExtra::diag_times3(shape2[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,g2);
	  one_eng = gayberne_lj(j,i,a2,b2,g2,r12,rsq,fforce,rtor);
	  ttor[0] = ttor[1] = ttor[2] = 0.0;
	  break;

        case ELLIPSE_SPHERE:
	  one_eng = gayberne_lj(i,j,a1,b1,g1,r12,rsq,fforce,ttor);
	  rtor[0] = rtor[1] = rtor[2] = 0.0;
	  break;

	default:
	  jquat = bonus[ellipsoid[j]].quat;
	  MathExtra::quat_to_mat_trans(jquat,a2);
	  MathExtra::diag_times3(well[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,b2);
	  MathExtra::diag_times3(shape2[jtype],a2,temp);
	  MathExtra::transpose_times3(a2,temp,g2);
	  one_eng = gayberne_analytic(i,j,a1,a2,b1,b2,g1,g2,r12,rsq,
				      fforce,ttor,rtor);
	  break;
	}

        fforce[0] *= factor_lj;
	fforce[1] *= factor_lj;
	fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
	ttor[1] *= factor_lj;
	ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
	f[i][1] += fforce[1];
	f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
	tor[i][1] += ttor[1];
	tor[i][2] += ttor[2];

        if (NEWTON_PAIR || j < nlocal) {
          rtor[0] *= factor_lj;
	  rtor[1] *= factor_lj;
	  rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
	  f[j][1] -= fforce[1];
	  f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
	  tor[j][1] += rtor[1];
	  tor[j][2] += rtor[2];
        }

        if (EFLAG) evdwl = factor_lj*one_eng;

	if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
				     evdwl,0.0,fforce[0],fforce[1],fforce[2],
				     -r12[0],-r12[1],-r12[2],tid);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairGayBerneOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairGayBerne::memory_usage();

  return bytes;
}
