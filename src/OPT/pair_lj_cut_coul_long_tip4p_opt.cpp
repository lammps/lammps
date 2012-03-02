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
#include "error.h"
#include "memory.h"
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

  // for caching m-shift corrected positions
  maxmpos = 0;
  h1idx = h2idx = NULL;
  mpos = NULL;
}

PairLJCutCoulLongTIP4POpt::~PairLJCutCoulLongTIP4POpt()
{
  memory->destroy(h1idx);
  memory->destroy(h2idx);
  memory->destroy(mpos);
}

void PairLJCutCoulLongTIP4POpt::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  // reallocate per-atom arrays, if necessary
  if (nall > maxmpos) {
    maxmpos = nall;
    memory->grow(mpos,maxmpos,3,"pair:mpos");
    memory->grow(h1idx,maxmpos,"pair:h1idx");
    memory->grow(h2idx,maxmpos,"pair:h2idx");
  }

  // cache corrected M positions in mpos[]
  const double * const * const x = atom->x;
  const int * const type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,h1idx[i],h2idx[i],mpos[i]);
    } else {
      mpos[i][0] = x[i][0];
      mpos[i][1] = x[i][1];
      mpos[i][2] = x[i][2];
    }
  }
  for (int i = nlocal; i < nall; i++) {
    if (type[i] == typeO) {
      find_M_permissive(i,h1idx[i],h2idx[i],mpos[i]);
    } else {
      mpos[i][0] = x[i][0];
      mpos[i][1] = x[i][1];
      mpos[i][2] = x[i][2];
    }
  }

  if (!ncoultablebits) {
    if (evflag) {
      if (eflag) {
	if (vflag) return eval<1,1,1,1>();
	else return eval<1,1,1,0>();
      } else {
	if (vflag) return eval<1,1,0,1>();
	else return eval<1,1,0,0>();
      }
    } else return eval<1,0,0,0>();
  } else {
    if (evflag) {
      if (eflag) {
	if (vflag) return eval<0,1,1,1>();
	else return eval<0,1,1,0>();
      } else {
	if (vflag) return eval<0,1,0,1>();
	else return eval<0,1,0,0>();
      }
    } else return eval<0,0,0,0>();
  }
}

/* ---------------------------------------------------------------------- */

template < const int CTABLE, const int EVFLAG, 
	   const int EFLAG, const int VFLAG>
void PairLJCutCoulLongTIP4POpt::eval()
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double delxOM,delyOM,delzOM;
  double r,rsq,r2inv,r6inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,ddotf;
  double v[6],xH1[3],xH2[3];
  double fdx,fdy,fdz,f1x,f1y,f1z,fOx,fOy,fOz,fHx,fHy,fHz;
  double *x1,*x2;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);

  double fxtmp,fytmp,fztmp;

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

    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;
    x1 = mpos[i];
    iH1 = h1idx[i];
    iH2 = h2idx[i];

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

      // LJ interaction based on true rsq

      if (rsq < cut_ljsq[itype][jtype]) {
	r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	forcelj *= factor_lj * r2inv;

	fxtmp += delx*forcelj;
	fytmp += dely*forcelj;
	fztmp += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;

	if (EFLAG) {
	  evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	    offset[itype][jtype];
	  evdwl *= factor_lj;
	} else evdwl = 0.0;

	if (EVFLAG) ev_tally(i,j,nlocal,/* newton_pair = */ 1,
			     evdwl,0.0,forcelj,delx,dely,delz);
      }

      // adjust rsq and delxyz for off-site O charge(s),
      // but only if they are within reach

      if (rsq < cut_coulsqplus) {

	if (itype == typeO || jtype == typeO) {
	  x2 = mpos[j];
	  jH1 = h1idx[j];
	  jH2 = h2idx[j];
	  if (jtype == typeO && ( jH1 < 0 || jH2 < 0))
	    error->one(FLERR,"TIP4P hydrogen is missing");
	  delx = x1[0] - x2[0];
	  dely = x1[1] - x2[1];
	  delz = x1[2] - x2[2];
	  rsq = delx*delx + dely*dely + delz*delz;
	}
      
	// Coulombic interaction based on modified rsq

	if (rsq < cut_coulsq) {
	  r2inv = 1 / rsq;
	  if (CTABLE || rsq <= tabinnersq) {
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
	    fxtmp += delx * cforce;
	    fytmp += dely * cforce;
	    fztmp += delz * cforce;

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

	    fdx = delx*cforce;
	    fdy = dely*cforce;
	    fdz = delz*cforce;

	    delxOM = x[i][0] - x1[0];
	    delyOM = x[i][1] - x1[1];
	    delzOM = x[i][2] - x1[2];

	    ddotf = (delxOM * fdx + delyOM * fdy + delzOM * fdz) /
	      (qdist*qdist);

	    f1x = alpha * (fdx - ddotf * delxOM);
	    f1y = alpha * (fdy - ddotf * delyOM);
	    f1z = alpha * (fdz - ddotf * delzOM);

	    fOx = fdx - f1x;
	    fOy = fdy - f1y;
	    fOz = fdz - f1z;

	    fHx = 0.5 * f1x;
	    fHy = 0.5 * f1y;
	    fHz = 0.5 * f1z;

	    fxtmp += fOx;
	    fytmp += fOy;
	    fztmp += fOz;

	    f[iH1][0] += fHx;
	    f[iH1][1] += fHy;
	    f[iH1][2] += fHz;

	    f[iH2][0] += fHx;
	    f[iH2][1] += fHy;
	    f[iH2][2] += fHz;

	    if (VFLAG) {
	      domain->closest_image(x[i],x[iH1],xH1);
	      domain->closest_image(x[i],x[iH2],xH2);

	      v[0] = x[i][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
	      v[1] = x[i][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
	      v[2] = x[i][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
	      v[3] = x[i][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
	      v[4] = x[i][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
	      v[5] = x[i][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;

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

	    fdx = -delx*cforce;
	    fdy = -dely*cforce;
	    fdz = -delz*cforce;

	    delxOM = x[j][0] - x2[0];
	    delyOM = x[j][1] - x2[1];
	    delzOM = x[j][2] - x2[2];

	    ddotf = (delxOM * fdx + delyOM * fdy + delzOM * fdz) /
	      (qdist*qdist);

	    f1x = alpha * (fdx - ddotf * delxOM);
	    f1y = alpha * (fdy - ddotf * delyOM);
	    f1z = alpha * (fdz - ddotf * delzOM);

	    fOx = fdx - f1x;
	    fOy = fdy - f1y;
	    fOz = fdz - f1z;

	    fHx = 0.5 * f1x;
	    fHy = 0.5 * f1y;
	    fHz = 0.5 * f1z;

	    f[j][0] += fOx;
	    f[j][1] += fOy;
	    f[j][2] += fOz;

	    f[jH1][0] += fHx;
	    f[jH1][1] += fHy;
	    f[jH1][2] += fHz;

	    f[jH2][0] += fHx;
	    f[jH2][1] += fHy;
	    f[jH2][2] += fHz;

	    if (VFLAG) {
	      domain->closest_image(x[j],x[jH1],xH1);
	      domain->closest_image(x[j],x[jH2],xH2);

	      v[0] += x[j][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
	      v[1] += x[j][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
	      v[2] += x[j][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
	      v[3] += x[j][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
	      v[4] += x[j][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
	      v[5] += x[j][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;

	      vlist[n++] = j;
	      vlist[n++] = jH1;
	      vlist[n++] = jH2;
	    }
	  }

	  if (EFLAG) {
	    if (CTABLE || rsq <= tabinnersq)
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
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4POpt::find_M_permissive(int i, int &iH1, int &iH2, double *xM)
{
  // test that O is correctly bonded to 2 succesive H atoms

   iH1 = atom->map(atom->tag[i] + 1);
   iH2 = atom->map(atom->tag[i] + 2);

   if (iH1 == -1 || iH2 == -1)
      return;
   else
      find_M(i,iH1,iH2,xM);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulLongTIP4POpt::memory_usage()
{
  double bytes = PairLJCutCoulLongTIP4P::memory_usage();
  bytes += 2 * maxmpos * sizeof(int);
  bytes += 3 * maxmpos * sizeof(double);
  bytes += maxmpos * sizeof(double *);

  return bytes;
}
