/* -------------------------------------------------------------------------
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
   OpenMP based threading support for LAMMPS
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include "thr_omp.h"

#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ThrOMP::ThrOMP(LAMMPS *ptr, int style) : lmp(ptr), fix(NULL), thr_style(style)
{
  // register fix omp with this class
  int ifix = lmp->modify->find_fix("package_omp");
  if (ifix < 0)
    lmp->error->all(FLERR,"The 'package omp' command is required for /omp styles");
  fix = static_cast<FixOMP *>(lmp->modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

ThrOMP::~ThrOMP()
{
  // nothing to do?
}

/* ----------------------------------------------------------------------
   Hook up per thread per atom arrays into the tally infrastructure
   ---------------------------------------------------------------------- */

void ThrOMP::ev_setup_thr(int eflag, int vflag, int nall, double *eatom,
			  double **vatom, ThrData *thr)
{
  const int tid = thr->get_tid();
  
  if (eflag & 2)
    thr->_eatom = eatom + tid*nall;

  if (vflag & 4)
    thr->_vatom = vatom + tid*nall;
  
}

/* ----------------------------------------------------------------------
   Reduce per thread data into the regular structures
   ---------------------------------------------------------------------- */

void ThrOMP::reduce_thr(const int eflag, const int vflag, ThrData *const thr)
{
  const int nlocal = lmp->atom->nlocal;
  const int nghost = lmp->atom->nghost;
  const int nall = nlocal + nghost;
  const int nfirst = lmp->atom->nfirst;
  const int nthreads = lmp->comm->nthreads;
  const int evflag = eflag | vflag;
  
  const int tid = thr->get_tid();
  double **f = lmp->atom->f;
  double **x = lmp->atom->x;

  switch (thr_style) {

  case THR_PAIR: {
    Pair * const pair = lmp->force->pair;
    
    if (pair->vflag_fdotr) {
      sync_threads();

      if (lmp->neighbor->includegroup == 0)
	thr->virial_fdotr_compute(x, nlocal, nghost, -1);
      else
	thr->virial_fdotr_compute(x, nlocal, nghost, nfirst);
    }

    if (evflag) {
      sync_threads();
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
	if (eflag & 1) {
	  pair->eng_vdwl += thr->eng_vdwl;
	  pair->eng_coul += thr->eng_coul;
	}
	if (vflag & 3)
	for (int i=0; i < 6; ++i)
	  pair->virial[i] += thr->virial_pair[i];
      }
    }
  }
    break;

  case THR_BOND: {
    Bond *bond = lmp->force->bond;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      bond->energy += thr->eng_bond;
      for (int i=0; i < 6; ++i)
	bond->virial[i] += thr->virial_bond[i];
    }
  }
    break;

  case THR_ANGLE: {
    Angle *angle = lmp->force->angle;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      angle->energy += thr->eng_angle;
      for (int i=0; i < 6; ++i)
	angle->virial[i] += thr->virial_angle[i];
    }
  }
    break;

  case THR_DIHEDRAL: {
    Dihedral *dihedral = lmp->force->dihedral;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      dihedral->energy += thr->eng_dihed;
      for (int i=0; i < 6; ++i)
	dihedral->virial[i] += thr->virial_dihed[i];
    }
  }
    break;

  case THR_IMPROPER: {
    Improper *improper = lmp->force->improper;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      improper->energy += thr->eng_imprp;
      for (int i=0; i < 6; ++i)
	improper->virial[i] += thr->virial_imprp[i];
    }
  }
    break;

  case THR_KSPACE: {
    KSpace *kspace = lmp->force->kspace;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      kspace->energy += thr->eng_kspce;
      for (int i=0; i < 6; ++i)
	kspace->virial[i] += thr->virial_kspce[i];
    }
  }
    break;

  }
    
  if (thr_style == fix->last_omp_style) {
    sync_threads();
    data_reduce_thr(&(f[0][0]), nall, nthreads, 3, tid);
    if (lmp->atom->torque)
      data_reduce_thr(&(lmp->atom->torque[0][0]), nall, nthreads, 3, tid);
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and eng_coul into per thread global and per-atom accumulators
------------------------------------------------------------------------- */

void ThrOMP::e_tally_thr(Pair * const pair, const int i, const int j, 
			 const int nlocal, const int newton_pair,
			 const double evdwl, const double ecoul, ThrData * const thr)
{
  if (pair->eflag_global) {
    if (newton_pair) {
      thr->eng_vdwl += evdwl;
      thr->eng_coul += ecoul;
    } else {
      const double evdwlhalf = 0.5*evdwl;
      const double ecoulhalf = 0.5*ecoul;
      if (i < nlocal) {
	thr->eng_vdwl += evdwlhalf;
	thr->eng_coul += ecoulhalf;
      }
      if (j < nlocal) {
	thr->eng_vdwl += evdwlhalf;
	thr->eng_coul += ecoulhalf;
      }
    }
  }
  if (pair->eflag_atom) {
    const double epairhalf = 0.5 * (evdwl + ecoul);
    if (newton_pair || i < nlocal) thr->_eatom[i] += epairhalf;
    if (newton_pair || j < nlocal) thr->_eatom[j] += epairhalf;
  }
}

/* ----------------------------------------------------------------------
   tally virial into per thread global and per-atom accumulators
------------------------------------------------------------------------- */
void ThrOMP::v_tally_thr(Pair * const pair, const int i, const int j, 
			 const int nlocal, const int newton_pair,
			 const double * const v, ThrData * const thr)
{
  if (pair->vflag_global) {
    double * const va = thr->virial_pair;
    if (newton_pair) {
      va[0] += v[0];
      va[1] += v[1];
      va[2] += v[2];
      va[3] += v[3];
      va[4] += v[4];
      va[5] += v[5];
    } else {
      if (i < nlocal) {
	va[0] += 0.5*v[0];
	va[1] += 0.5*v[1];
	va[2] += 0.5*v[2];
	va[3] += 0.5*v[3];
	va[4] += 0.5*v[4];
	va[5] += 0.5*v[5];
      }
      if (j < nlocal) {
	va[0] += 0.5*v[0];
	va[1] += 0.5*v[1];
	va[2] += 0.5*v[2];
	va[3] += 0.5*v[3];
	va[4] += 0.5*v[4];
	va[5] += 0.5*v[5];
      }
    }
  }

  if (pair->vflag_atom) {
    if (newton_pair || i < nlocal) {
      double * const va = thr->_vatom[i];
      va[0] += 0.5*v[0];
      va[1] += 0.5*v[1];
      va[2] += 0.5*v[2];
      va[3] += 0.5*v[3];
      va[4] += 0.5*v[4];
      va[5] += 0.5*v[5];
    }
    if (newton_pair || j < nlocal) {
      double * const va = thr->_vatom[j];
      va[0] += 0.5*v[0];
      va[1] += 0.5*v[1];
      va[2] += 0.5*v[2];
      va[3] += 0.5*v[3];
      va[4] += 0.5*v[4];
      va[5] += 0.5*v[5];
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per thread global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Pair * const pair, const int i, const int j, const int nlocal,
			  const int newton_pair, const double evdwl, const double ecoul,
			  const double fpair, const double delx, const double dely,
			  const double delz, ThrData * const thr)
{

  if (pair->eflag_either)
    e_tally_thr(pair, i, j, nlocal, newton_pair, evdwl, ecoul, thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    v_tally_thr(pair, i, j, nlocal, newton_pair, v, thr);
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_xyz_thr(Pair * const pair, const int i, const int j,
			      const int nlocal, const int newton_pair, 
			      const double evdwl, const double ecoul,
			      const double fx, const double fy, const double fz,
			      const double delx, const double dely, const double delz,
			      ThrData * const thr)
{

  if (pair->eflag_either)
    e_tally_thr(pair, i, j, nlocal, newton_pair, evdwl, ecoul, thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    v_tally_thr(pair, i, j, nlocal, newton_pair, v, thr);
  }
}

#if 0
/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally3_thr(Pair *pair, int i, int j, int k, double evdwl, double ecoul,
			   double *fj, double *fk, double *drji, double *drki, int tid)
{
  double epairthird,v[6];

  if (pair->eflag_either) {
    if (pair->eflag_global) {
      eng_vdwl_thr[tid] += evdwl;
      eng_coul_thr[tid] += ecoul;
    }
    if (pair->eflag_atom) {
      epairthird = THIRD * (evdwl + ecoul);
      eatom_thr[tid][i] += epairthird;
      eatom_thr[tid][j] += epairthird;
      eatom_thr[tid][k] += epairthird;
    }
  }

  if (pair->vflag_either) {
    v[0] = drji[0]*fj[0] + drki[0]*fk[0];
    v[1] = drji[1]*fj[1] + drki[1]*fk[1];
    v[2] = drji[2]*fj[2] + drki[2]*fk[2];
    v[3] = drji[0]*fj[1] + drki[0]*fk[1];
    v[4] = drji[0]*fj[2] + drki[0]*fk[2];
    v[5] = drji[1]*fj[2] + drki[1]*fk[2];
      
    if (pair->vflag_global) {
      virial_thr[tid][0] += v[0];
      virial_thr[tid][1] += v[1];
      virial_thr[tid][2] += v[2];
      virial_thr[tid][3] += v[3];
      virial_thr[tid][4] += v[4];
      virial_thr[tid][5] += v[5];
    }

    if (pair->vflag_atom) {
      for (int n=0; n < 6; ++n) {
	vatom_thr[tid][i][n] += THIRD*v[n];
	vatom_thr[tid][j][n] += THIRD*v[n];
	vatom_thr[tid][k][n] += THIRD*v[n];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by AIREBO potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally4_thr(Pair *pair, int i, int j, int k, int m, double evdwl,
			   double *fi, double *fj, double *fk,
			   double *drim, double *drjm, double *drkm,int tid)
{
  double epairfourth,v[6];

  if (pair->eflag_either) {
    if (pair->eflag_global) eng_vdwl_thr[tid] += evdwl;
    if (pair->eflag_atom) {
      epairfourth = 0.25 * evdwl;
      eatom_thr[tid][i] += epairfourth;
      eatom_thr[tid][j] += epairfourth;
      eatom_thr[tid][k] += epairfourth;
      eatom_thr[tid][m] += epairfourth;
    }
  }

  if (pair->vflag_atom) {
    v[0] = 0.25 * (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
    v[1] = 0.25 * (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
    v[2] = 0.25 * (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
    v[3] = 0.25 * (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
    v[4] = 0.25 * (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
    v[5] = 0.25 * (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);
    
    vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1]; vatom_thr[tid][i][2] += v[2];
    vatom_thr[tid][i][3] += v[3]; vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
    vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1]; vatom_thr[tid][j][2] += v[2];
    vatom_thr[tid][j][3] += v[3]; vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
    vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1]; vatom_thr[tid][k][2] += v[2];
    vatom_thr[tid][k][3] += v[3]; vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
    vatom_thr[tid][m][0] += v[0]; vatom_thr[tid][m][1] += v[1]; vatom_thr[tid][m][2] += v[2];
    vatom_thr[tid][m][3] += v[3]; vatom_thr[tid][m][4] += v[4]; vatom_thr[tid][m][5] += v[5];
  }
}

/* ----------------------------------------------------------------------
   tally ecoul and virial into each of n atoms in list
   called by TIP4P potential, newton_pair is always on
   changes v values by dividing by n
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally_list_thr(Pair *pair, int n, int *list, double ecoul, double *v, int tid)
{
  int i,j;

  if (pair->eflag_either) {
    if (pair->eflag_global) eng_coul_thr[tid] += ecoul;
    if (pair->eflag_atom) {
      double epairatom = ecoul/n;
      for (i = 0; i < n; i++) eatom_thr[tid][list[i]] += epairatom;
    }
  }

  if (pair->vflag_either) {
    if (pair->vflag_global) {
      virial_thr[tid][0] += v[0];
      virial_thr[tid][1] += v[1];
      virial_thr[tid][2] += v[2];
      virial_thr[tid][3] += v[3];
      virial_thr[tid][4] += v[4];
      virial_thr[tid][5] += v[5];
    }

    if (pair->vflag_atom) {
      v[0] /= n;
      v[1] /= n;
      v[2] /= n;
      v[3] /= n;
      v[4] /= n;
      v[5] /= n;
      for (i = 0; i < n; i++) {
	j = list[i];
	vatom_thr[tid][j][0] += v[0];
	vatom_thr[tid][j][1] += v[1];
	vatom_thr[tid][j][2] += v[2];
	vatom_thr[tid][j][3] += v[3];
	vatom_thr[tid][j][4] += v[4];
	vatom_thr[tid][j][5] += v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
	  = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Dihedral *dihed, int i1, int i2, int i3, int i4,
			  int nlocal, int newton_bond,
			  double edihedral, double *f1, double *f3, double *f4,
			  double vb1x, double vb1y, double vb1z,
			  double vb2x, double vb2y, double vb2z,
			  double vb3x, double vb3y, double vb3z, int tid)
{
  double edihedralquarter,v[6];
  int cnt;

  if (dihed->eflag_either) {
    if (dihed->eflag_global) {
      if (newton_bond) {
	eng_bond_thr[tid] += edihedral;
      } else {
	edihedralquarter = 0.25*edihedral;
	cnt = 0;
	if (i1 < nlocal) ++cnt;
	if (i2 < nlocal) ++cnt;
	if (i3 < nlocal) ++cnt;
	if (i4 < nlocal) ++cnt;
	eng_bond_thr[tid] += static_cast<double>(cnt) * edihedralquarter;
      }
    }
    if (dihed->eflag_atom) {
      edihedralquarter = 0.25*edihedral;
      if (newton_bond || i1 < nlocal) eatom_thr[tid][i1] += edihedralquarter;
      if (newton_bond || i2 < nlocal) eatom_thr[tid][i2] += edihedralquarter;
      if (newton_bond || i3 < nlocal) eatom_thr[tid][i3] += edihedralquarter;
      if (newton_bond || i4 < nlocal) eatom_thr[tid][i4] += edihedralquarter;
    }
  }

  if (dihed->vflag_either) {
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (dihed->vflag_global) {
      if (newton_bond) {
	virial_thr[tid][0] += v[0];
	virial_thr[tid][1] += v[1];
	virial_thr[tid][2] += v[2];
	virial_thr[tid][3] += v[3];
	virial_thr[tid][4] += v[4];
	virial_thr[tid][5] += v[5];
      } else {
	if (i1 < nlocal) {
	  virial_thr[tid][0] += 0.25*v[0];
	  virial_thr[tid][1] += 0.25*v[1];
	  virial_thr[tid][2] += 0.25*v[2];
	  virial_thr[tid][3] += 0.25*v[3];
	  virial_thr[tid][4] += 0.25*v[4];
	  virial_thr[tid][5] += 0.25*v[5];
	}
	if (i2 < nlocal) {
	  virial_thr[tid][0] += 0.25*v[0];
	  virial_thr[tid][1] += 0.25*v[1];
	  virial_thr[tid][2] += 0.25*v[2];
	  virial_thr[tid][3] += 0.25*v[3];
	  virial_thr[tid][4] += 0.25*v[4];
	  virial_thr[tid][5] += 0.25*v[5];
	}
	if (i3 < nlocal) {
	  virial_thr[tid][0] += 0.25*v[0];
	  virial_thr[tid][1] += 0.25*v[1];
	  virial_thr[tid][2] += 0.25*v[2];
	  virial_thr[tid][3] += 0.25*v[3];
	  virial_thr[tid][4] += 0.25*v[4];
	  virial_thr[tid][5] += 0.25*v[5];
	}
	if (i4 < nlocal) {
	  virial_thr[tid][0] += 0.25*v[0];
	  virial_thr[tid][1] += 0.25*v[1];
	  virial_thr[tid][2] += 0.25*v[2];
	  virial_thr[tid][3] += 0.25*v[3];
	  virial_thr[tid][4] += 0.25*v[4];
	  virial_thr[tid][5] += 0.25*v[5];
	}
      }
    }

    if (dihed->vflag_atom) {
      if (newton_bond || i1 < nlocal) {
	vatom_thr[tid][i1][0] += 0.25*v[0];
	vatom_thr[tid][i1][1] += 0.25*v[1];
	vatom_thr[tid][i1][2] += 0.25*v[2];
	vatom_thr[tid][i1][3] += 0.25*v[3];
	vatom_thr[tid][i1][4] += 0.25*v[4];
	vatom_thr[tid][i1][5] += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
	vatom_thr[tid][i2][0] += 0.25*v[0];
	vatom_thr[tid][i2][1] += 0.25*v[1];
	vatom_thr[tid][i2][2] += 0.25*v[2];
	vatom_thr[tid][i2][3] += 0.25*v[3];
	vatom_thr[tid][i2][4] += 0.25*v[4];
	vatom_thr[tid][i2][5] += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
	vatom_thr[tid][i3][0] += 0.25*v[0];
	vatom_thr[tid][i3][1] += 0.25*v[1];
	vatom_thr[tid][i3][2] += 0.25*v[2];
	vatom_thr[tid][i3][3] += 0.25*v[3];
	vatom_thr[tid][i3][4] += 0.25*v[4];
	vatom_thr[tid][i3][5] += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
	vatom_thr[tid][i4][0] += 0.25*v[0];
	vatom_thr[tid][i4][1] += 0.25*v[1];
	vatom_thr[tid][i4][2] += 0.25*v[2];
	vatom_thr[tid][i4][3] += 0.25*v[3];
	vatom_thr[tid][i4][4] += 0.25*v[4];
	vatom_thr[tid][i4][5] += 0.25*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void ThrOMP::v_tally2_thr(int i, int j, double fpair, double *drij, int tid)
{
  double v[6];
  
  v[0] = 0.5 * drij[0]*drij[0]*fpair;
  v[1] = 0.5 * drij[1]*drij[1]*fpair;
  v[2] = 0.5 * drij[2]*drij[2]*fpair;
  v[3] = 0.5 * drij[0]*drij[1]*fpair;
  v[4] = 0.5 * drij[0]*drij[2]*fpair;
  v[5] = 0.5 * drij[1]*drij[2]*fpair;

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1]; vatom_thr[tid][i][2] += v[2];
  vatom_thr[tid][i][3] += v[3]; vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1]; vatom_thr[tid][j][2] += v[2];
  vatom_thr[tid][j][3] += v[3]; vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void ThrOMP::v_tally3_thr(int i, int j, int k, double *fi, double *fj,
			  double *drik, double *drjk, int tid)
{
  double v[6];
  
  v[0] = THIRD * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = THIRD * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = THIRD * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = THIRD * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = THIRD * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = THIRD * (drik[1]*fi[2] + drjk[1]*fj[2]);

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1]; vatom_thr[tid][i][2] += v[2];
  vatom_thr[tid][i][3] += v[3]; vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1]; vatom_thr[tid][j][2] += v[2];
  vatom_thr[tid][j][3] += v[3]; vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
  vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1]; vatom_thr[tid][k][2] += v[2];
  vatom_thr[tid][k][3] += v[3]; vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
------------------------------------------------------------------------- */

void ThrOMP::v_tally4_thr(int i, int j, int k, int m,
			  double *fi, double *fj, double *fk,
			  double *drim, double *drjm, double *drkm, int tid)
{
  double v[6];

  v[0] = 0.25 * (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
  v[1] = 0.25 * (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
  v[2] = 0.25 * (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
  v[3] = 0.25 * (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
  v[4] = 0.25 * (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
  v[5] = 0.25 * (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1]; vatom_thr[tid][i][2] += v[2];
  vatom_thr[tid][i][3] += v[3]; vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1]; vatom_thr[tid][j][2] += v[2];
  vatom_thr[tid][j][3] += v[3]; vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
  vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1]; vatom_thr[tid][k][2] += v[2];
  vatom_thr[tid][k][3] += v[3]; vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
  vatom_thr[tid][m][0] += v[0]; vatom_thr[tid][m][1] += v[1]; vatom_thr[tid][m][2] += v[2];
  vatom_thr[tid][m][3] += v[3]; vatom_thr[tid][m][4] += v[4]; vatom_thr[tid][m][5] += v[5];
}

#endif

/* ---------------------------------------------------------------------- */

double ThrOMP::memory_usage_thr() 
{
  const int nthreads=lmp->comm->nthreads;
  double bytes=0.0;
  
  return bytes;
}
