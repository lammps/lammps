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
   OPT version: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_lj_cut_coul_long_tip4p_opt.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutCoulLongTIP4POpt::PairLJCutCoulLongTIP4POpt(LAMMPS *lmp) : 
  PairLJCutCoulLongTIP4P(lmp)
{
  single_enable = 0;
  respa_enable = 0;

  // TIP4P cannot compute virial as F dot r
  // due to find_M() finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

void PairLJCutCoulLongTIP4POpt::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  if (!ncoultablebits) {
    if (evflag) {
      if (eflag) {
	if (vflag) {
	  if (force->newton_pair) return eval<1,1,1,1,0>();
	  else return eval<1,1,1,0,0>();
	} else {
	  if (force->newton_pair) return eval<1,1,0,1,0>();
	  else return eval<1,1,0,0,0>();
	} 
      } else {
	if (vflag) {
	  if (force->newton_pair) return eval<1,0,1,1,0>();
	  else return eval<1,0,1,0,0>();
	} else {
	  if (force->newton_pair) return eval<1,0,0,1,0>();
	  else return eval<1,0,0,0,0>();
	}
      }
    } else {
      if (force->newton_pair) return eval<0,0,0,1,0>();
      else return eval<0,0,0,0,0>();
    }
  } else {
    if (evflag) {
      if (eflag) {
	if (vflag) {
	  if (force->newton_pair) return eval<1,1,1,1,1>();
	  else return eval<1,1,1,0,1>();
	} else {
	  if (force->newton_pair) return eval<1,1,0,1,1>();
	  else return eval<1,1,0,0,1>();
	}
      } else {
	if (vflag) {
	  if (force->newton_pair) return eval<1,0,1,1,1>();
	  else return eval<1,0,1,0,1>();
	} else {
	  if (force->newton_pair) return eval<1,0,0,1,1>();
	  else return eval<1,0,0,0,1>();
	}
      }
    } else {
      if (force->newton_pair) return eval<0,0,0,1,1>();
      else return eval<0,0,0,0,1>();
    }
  }
}

/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG, const int VFLAG,
           const int NEWTON_PAIR, const int CTABLE >
void PairLJCutCoulLongTIP4POpt::eval()
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double delxOM, delyOM, delzOM;
  double r,r2inv,r6inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,ddotf;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],f1[3],v[6],xH1[3],xH2[3];
  double *x1,*x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    if (itype == typeO) {
      find_M(i,iH1,iH2,xiM);
      x1 = xiM;
    } else x1 = x[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cutsq[itype][jtype]) {

	r2inv = 1.0/rsq;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  forcelj *= factor_lj * r2inv;

	  f[i][0] += delx*forcelj;
	  f[i][1] += dely*forcelj;
	  f[i][2] += delz*forcelj;
	  f[j][0] -= delx*forcelj;
	  f[j][1] -= dely*forcelj;
	  f[j][2] -= delz*forcelj;

	  if (EFLAG) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;

	  if (EVFLAG) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,forcelj,delx,dely,delz);
	}

	// adjust rsq and delxyz for off-site O charge(s)

	if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
	    find_M(j,jH1,jH2,xjM);
	    x2 = xjM;
	  } else x2 = x[j];
	  delx = x1[0] - x2[0];
	  dely = x1[1] - x2[1];
	  delz = x1[2] - x2[2];
	  rsq = delx*delx + dely*dely + delz*delz;
	}

	// test current rsq against cutoff and compute Coulombic force

	if (rsq < cut_coulsq) {
	  r2inv = 1 / rsq;
	  if (!CTABLE || rsq <= tabinnersq) {
	    r = sqrt(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    if (factor_coul < 1.0) {
	      forcecoul -= (1.0-factor_coul)*prefactor; 
	    }
	  } else {
	    union_int_float_t rsq_lookup;
	    rsq_lookup.f = rsq;
	    itable = rsq_lookup.i & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }

	  cforce = forcecoul * r2inv;

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  n = 0;

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

	    if (VFLAG) {
	      v[0] = x[i][0] * delx * cforce;
	      v[1] = x[i][1] * dely * cforce;
	      v[2] = x[i][2] * delz * cforce;
	      v[3] = x[i][0] * dely * cforce;
	      v[4] = x[i][0] * delz * cforce;
	      v[5] = x[i][1] * delz * cforce;
	      vlist[n++] = i;
	    }

	  } else {

            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            delxOM = x[i][0] - x1[0];
            delyOM = x[i][1] - x1[1];
            delzOM = x[i][2] - x1[2];

            ddotf = (delxOM * fd[0] + delyOM * fd[1] + delzOM * fd[2]) /
            (qdist*qdist);

            f1[0] = ddotf * delxOM;
            f1[1] = ddotf * delyOM;
            f1[2] = ddotf * delzOM;

            fO[0] = fd[0] - alpha * (fd[0] - f1[0]);
            fO[1] = fd[1] - alpha * (fd[1] - f1[1]);
            fO[2] = fd[2] - alpha * (fd[2] - f1[2]);

            fH[0] = 0.5 * alpha * (fd[0] - f1[0]);
            fH[1] = 0.5 * alpha * (fd[1] - f1[1]);
            fH[2] = 0.5 * alpha * (fd[2] - f1[2]);

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

	    if (VFLAG) {
	      domain->closest_image(x[i],x[iH1],xH1);
	      domain->closest_image(x[i],x[iH2],xH2);

	      v[0] = x[i][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] = x[i][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] = x[i][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] = x[i][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] = x[i][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] = x[i][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];

	      vlist[n++] = i;
	      vlist[n++] = iH1;
	      vlist[n++] = iH2;
	    }
	  }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	    if (VFLAG) {
	      v[0] -= x[j][0] * delx * cforce;
	      v[1] -= x[j][1] * dely * cforce;
	      v[2] -= x[j][2] * delz * cforce;
	      v[3] -= x[j][0] * dely * cforce;
	      v[4] -= x[j][0] * delz * cforce;
	      v[5] -= x[j][1] * delz * cforce;
	      vlist[n++] = j;
	    }

	  } else {

	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

	    delxOM = x[j][0] - x2[0];
	    delyOM = x[j][1] - x2[1];
	    delzOM = x[j][2] - x2[2];

	    ddotf = (delxOM * fd[0] + delyOM * fd[1] + delzOM * fd[2]) /
	      (qdist*qdist);

	    f1[0] = ddotf * delxOM;
	    f1[1] = ddotf * delyOM;
	    f1[2] = ddotf * delzOM;

	    fO[0] = fd[0] - alpha * (fd[0] - f1[0]);
	    fO[1] = fd[1] - alpha * (fd[1] - f1[1]);
	    fO[2] = fd[2] - alpha * (fd[2] - f1[2]);

	    fH[0] = 0.5 * alpha * (fd[0] - f1[0]);
	    fH[1] = 0.5 * alpha * (fd[1] - f1[1]);
	    fH[2] = 0.5 * alpha * (fd[2] - f1[2]);

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];

	    if (VFLAG) {
	      domain->closest_image(x[j],x[jH1],xH1);
	      domain->closest_image(x[j],x[jH2],xH2);

	      v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];

	      vlist[n++] = j;
	      vlist[n++] = jH1;
	      vlist[n++] = jH2;
	    }
	  }

	  if (EFLAG) {
	    if (!CTABLE || rsq <= tabinnersq)
	      ecoul = prefactor*erfc;
	    else {
	      table = etable[itable] + fraction*detable[itable];
	      ecoul = qtmp*q[j] * table;
	    }
	    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	  } else ecoul = 0.0;

	  if (EVFLAG) ev_tally_list(n,vlist,ecoul,v);
	}
      }
    }
  }
}

