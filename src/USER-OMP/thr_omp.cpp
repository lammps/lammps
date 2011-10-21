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

#include "thr_omp.h"

#include "memory.h"

#include "atom.h"
#include "comm.h"
#include "force.h"

#include "pair.h"
#include "dihedral.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ThrOMP::ThrOMP(LAMMPS *ptr, int style) : thr_style(style), lmp(ptr)
{
  // initialize fixed size per thread storage
  eng_vdwl_thr = eng_coul_thr = eng_bond_thr = NULL;
  virial_thr = NULL;

  lmp->memory->create(eng_vdwl_thr,lmp->comm->nthreads,"thr_omp:eng_vdwl_thr");
  lmp->memory->create(eng_coul_thr,lmp->comm->nthreads,"thr_omp:eng_coul_thr");
  lmp->memory->create(eng_bond_thr,lmp->comm->nthreads,"thr_omp:eng_bond_thr");
  lmp->memory->create(virial_thr,lmp->comm->nthreads,6,"thr_omp:virial_thr");

  // variable size per thread, per atom storage
  // the actually allocation happens via memory->grow() in ev_steup_thr()
  maxeatom_thr = maxvatom_thr = 0;
  evflag_global = evflag_atom = 0;
  eatom_thr = NULL;
  vatom_thr = NULL;
}

/* ---------------------------------------------------------------------- */

ThrOMP::~ThrOMP()
{
  lmp->memory->destroy(eng_vdwl_thr);
  lmp->memory->destroy(eng_coul_thr);
  lmp->memory->destroy(eng_bond_thr);
  lmp->memory->destroy(virial_thr);
  lmp->memory->destroy(eatom_thr);
  lmp->memory->destroy(vatom_thr);
}

/* ---------------------------------------------------------------------- */

void ThrOMP::ev_setup_acc_thr(int ntotal, int eflag_global, int vflag_global,
			     int eflag_atom, int vflag_atom, int nthreads)
{
  int t,i;

  evflag_global = (eflag_global || vflag_global);
  evflag_atom = (eflag_atom || vflag_atom);
  
  for (t = 0; t < nthreads; ++t) {

    if (eflag_global) 
      eng_vdwl_thr[t] = eng_coul_thr[t] = eng_bond_thr[t] = 0.0;

    if (vflag_global) 
      for (i = 0; i < 6; ++i)
	virial_thr[t][i] = 0.0;

    if (eflag_atom)
      for (i = 0; i < ntotal; ++i)
	eatom_thr[t][i] = 0.0;
    
    if (vflag_atom)
      for (i = 0; i < ntotal; ++i) {
        vatom_thr[t][i][0] = 0.0;
        vatom_thr[t][i][1] = 0.0;
        vatom_thr[t][i][2] = 0.0;
        vatom_thr[t][i][3] = 0.0;
        vatom_thr[t][i][4] = 0.0;
        vatom_thr[t][i][5] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void ThrOMP::ev_setup_thr(Dihedral *dihed)
{
  int nthreads = lmp->comm->nthreads;

  // reallocate per-atom arrays if necessary
  if (dihed->eflag_atom && lmp->atom->nmax > maxeatom_thr) {
    maxeatom_thr = lmp->atom->nmax;
    lmp->memory->grow(eatom_thr,nthreads,maxeatom_thr,"thr_omp:eatom_thr");
  }
  if (dihed->vflag_atom && lmp->atom->nmax > maxvatom_thr) {
    maxvatom_thr = lmp->atom->nmax;
    lmp->memory->grow(vatom_thr,nthreads,maxeatom_thr,6,"thr_omp:vatom_thr");
  }

  int ntotal = (lmp->force->newton_bond) ? 
    (lmp->atom->nlocal + lmp->atom->nghost) : lmp->atom->nlocal;

  // set up per thread accumulators
  ev_setup_acc_thr(ntotal, dihed->eflag_global, dihed->vflag_global,
		   dihed->eflag_atom, dihed->vflag_atom, nthreads);
}

/* ---------------------------------------------------------------------- */

void ThrOMP::ev_setup_thr(Pair *pair)
{
  int nthreads = lmp->comm->nthreads;

  // reallocate per-atom arrays if necessary
  if (pair->eflag_atom && lmp->atom->nmax > maxeatom_thr) {
    maxeatom_thr = lmp->atom->nmax;
    lmp->memory->grow(eatom_thr,nthreads,maxeatom_thr,"thr_omp:eatom_thr");
  }
  if (pair->vflag_atom && lmp->atom->nmax > maxvatom_thr) {
    maxvatom_thr = lmp->atom->nmax;
    lmp->memory->grow(vatom_thr,nthreads,maxeatom_thr,6,"thr_omp:vatom_thr");
  }

  int ntotal = (lmp->force->newton) ?
    (lmp->atom->nlocal + lmp->atom->nghost) : lmp->atom->nlocal;

  // set up per thread accumulators
  ev_setup_acc_thr(ntotal, pair->eflag_global, pair->vflag_global,
		   pair->eflag_atom, pair->vflag_atom, nthreads);
}

/* ----------------------------------------------------------------------
   reduce the per thread accumulated E/V data into the canonical accumulators.
------------------------------------------------------------------------- */
void ThrOMP::ev_reduce_thr(Dihedral *dihed)
{
  int nthreads = lmp->comm->nthreads;
  int ntotal = (lmp->force->newton_bond) ?
    (lmp->atom->nlocal + lmp->atom->nghost) : lmp->atom->nlocal;

  for (int n = 0; n < nthreads; ++n) {
    dihed->energy += eng_bond_thr[n];
    if (dihed->vflag_either) {
      dihed->virial[0] += virial_thr[n][0];
      dihed->virial[1] += virial_thr[n][1];
      dihed->virial[2] += virial_thr[n][2];
      dihed->virial[3] += virial_thr[n][3];
      dihed->virial[4] += virial_thr[n][4];
      dihed->virial[5] += virial_thr[n][5];
      if (dihed->vflag_atom) {
        for (int i = 0; i < ntotal; ++i) {
          dihed->vatom[i][0] += vatom_thr[n][i][0];
          dihed->vatom[i][1] += vatom_thr[n][i][1];
          dihed->vatom[i][2] += vatom_thr[n][i][2];
          dihed->vatom[i][3] += vatom_thr[n][i][3];
          dihed->vatom[i][4] += vatom_thr[n][i][4];
          dihed->vatom[i][5] += vatom_thr[n][i][5];
        }
      }
    }
    if (dihed->eflag_atom) {
      for (int i = 0; i < ntotal; ++i) {
        dihed->eatom[i] += eatom_thr[n][i];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reduce the per thread accumulated E/V data into the canonical accumulators.
------------------------------------------------------------------------- */
void ThrOMP::ev_reduce_thr(Pair *pair)
{
  const int nthreads = lmp->comm->nthreads;
  const int ntotal = (lmp->force->newton) ? 
    (lmp->atom->nlocal + lmp->atom->nghost) : lmp->atom->nlocal;

  for (int n = 0; n < nthreads; ++n) {
    pair->eng_vdwl += eng_vdwl_thr[n];
    pair->eng_coul += eng_coul_thr[n];
    if (pair->vflag_either) {
      pair->virial[0] += virial_thr[n][0];
      pair->virial[1] += virial_thr[n][1];
      pair->virial[2] += virial_thr[n][2];
      pair->virial[3] += virial_thr[n][3];
      pair->virial[4] += virial_thr[n][4];
      pair->virial[5] += virial_thr[n][5];
      if (pair->vflag_atom) {
        for (int i = 0; i < ntotal; ++i) {
          pair->vatom[i][0] += vatom_thr[n][i][0];
          pair->vatom[i][1] += vatom_thr[n][i][1];
          pair->vatom[i][2] += vatom_thr[n][i][2];
          pair->vatom[i][3] += vatom_thr[n][i][3];
          pair->vatom[i][4] += vatom_thr[n][i][4];
          pair->vatom[i][5] += vatom_thr[n][i][5];
        }
      }
    }
    if (pair->eflag_atom) {
      for (int i = 0; i < ntotal; ++i) {
        pair->eatom[i] += eatom_thr[n][i];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per thread global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Pair *pair, int i, int j, int nlocal,
			  int newton_pair, double evdwl, double ecoul,
			  double fpair, double delx, double dely,
			  double delz, int tid)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (pair->eflag_either) {
    if (pair->eflag_global) {
      if (newton_pair) {
	eng_vdwl_thr[tid] += evdwl;
	eng_coul_thr[tid] += ecoul;
      } else {
	evdwlhalf = 0.5*evdwl;
	ecoulhalf = 0.5*ecoul;
	if (i < nlocal) {
	  eng_vdwl_thr[tid] += evdwlhalf;
	  eng_coul_thr[tid] += ecoulhalf;
	}
	if (j < nlocal) {
	  eng_vdwl_thr[tid] += evdwlhalf;
	  eng_coul_thr[tid] += ecoulhalf;
	}
      }
    }
    if (pair->eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom_thr[tid][i] += epairhalf;
      if (newton_pair || j < nlocal) eatom_thr[tid][j] += epairhalf;
    }
  }

  if (pair->vflag_either) {
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (pair->vflag_global) {
      if (newton_pair) {
	virial_thr[tid][0] += v[0];
	virial_thr[tid][1] += v[1];
	virial_thr[tid][2] += v[2];
	virial_thr[tid][3] += v[3];
	virial_thr[tid][4] += v[4];
	virial_thr[tid][5] += v[5];
      } else {
	if (i < nlocal) {
	  virial_thr[tid][0] += 0.5*v[0];
	  virial_thr[tid][1] += 0.5*v[1];
	  virial_thr[tid][2] += 0.5*v[2];
	  virial_thr[tid][3] += 0.5*v[3];
	  virial_thr[tid][4] += 0.5*v[4];
	  virial_thr[tid][5] += 0.5*v[5];
	}
	if (j < nlocal) {
	  virial_thr[tid][0] += 0.5*v[0];
	  virial_thr[tid][1] += 0.5*v[1];
	  virial_thr[tid][2] += 0.5*v[2];
	  virial_thr[tid][3] += 0.5*v[3];
	  virial_thr[tid][4] += 0.5*v[4];
	  virial_thr[tid][5] += 0.5*v[5];
	}
      }
    }

    if (pair->vflag_atom) {
      if (newton_pair || i < nlocal) {
	vatom_thr[tid][i][0] += 0.5*v[0];
	vatom_thr[tid][i][1] += 0.5*v[1];
	vatom_thr[tid][i][2] += 0.5*v[2];
	vatom_thr[tid][i][3] += 0.5*v[3];
	vatom_thr[tid][i][4] += 0.5*v[4];
	vatom_thr[tid][i][5] += 0.5*v[5];
      }
      if (newton_pair || j < nlocal) {
	vatom_thr[tid][j][0] += 0.5*v[0];
	vatom_thr[tid][j][1] += 0.5*v[1];
	vatom_thr[tid][j][2] += 0.5*v[2];
	vatom_thr[tid][j][3] += 0.5*v[3];
	vatom_thr[tid][j][4] += 0.5*v[4];
	vatom_thr[tid][j][5] += 0.5*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_xyz_thr(Pair *pair, int i, int j, int nlocal,
			      int newton_pair, double evdwl, double ecoul,
			      double fx, double fy, double fz,
			      double delx, double dely, double delz, int tid)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (pair->eflag_either) {
    if (pair->eflag_global) {
      if (newton_pair) {
	eng_vdwl_thr[tid] += evdwl;
	eng_coul_thr[tid] += ecoul;
      } else {
	evdwlhalf = 0.5*evdwl;
	ecoulhalf = 0.5*ecoul;
	if (i < nlocal) {
	  eng_vdwl_thr[tid] += evdwlhalf;
	  eng_coul_thr[tid] += ecoulhalf;
	}
	if (j < nlocal) {
	  eng_vdwl_thr[tid] += evdwlhalf;
	  eng_coul_thr[tid] += ecoulhalf;
	}
      }
    }
    if (pair->eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom_thr[tid][i] += epairhalf;
      if (newton_pair || j < nlocal) eatom_thr[tid][j] += epairhalf;
    }
  }

  if (pair->vflag_either) {
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    if (pair->vflag_global) {
      if (newton_pair) {
	virial_thr[tid][0] += v[0];
	virial_thr[tid][1] += v[1];
	virial_thr[tid][2] += v[2];
	virial_thr[tid][3] += v[3];
	virial_thr[tid][4] += v[4];
	virial_thr[tid][5] += v[5];
      } else {
	if (i < nlocal) {
	  virial_thr[tid][0] += 0.5*v[0];
	  virial_thr[tid][1] += 0.5*v[1];
	  virial_thr[tid][2] += 0.5*v[2];
	  virial_thr[tid][3] += 0.5*v[3];
	  virial_thr[tid][4] += 0.5*v[4];
	  virial_thr[tid][5] += 0.5*v[5];
	}
	if (j < nlocal) {
	  virial_thr[tid][0] += 0.5*v[0];
	  virial_thr[tid][1] += 0.5*v[1];
	  virial_thr[tid][2] += 0.5*v[2];
	  virial_thr[tid][3] += 0.5*v[3];
	  virial_thr[tid][4] += 0.5*v[4];
	  virial_thr[tid][5] += 0.5*v[5];
	}
      }
    }

    if (pair->vflag_atom) {
      if (newton_pair || i < nlocal) {
	vatom_thr[tid][i][0] += 0.5*v[0];
	vatom_thr[tid][i][1] += 0.5*v[1];
	vatom_thr[tid][i][2] += 0.5*v[2];
	vatom_thr[tid][i][3] += 0.5*v[3];
	vatom_thr[tid][i][4] += 0.5*v[4];
	vatom_thr[tid][i][5] += 0.5*v[5];
      }
      if (newton_pair || j < nlocal) {
	vatom_thr[tid][j][0] += 0.5*v[0];
	vatom_thr[tid][j][1] += 0.5*v[1];
	vatom_thr[tid][j][2] += 0.5*v[2];
	vatom_thr[tid][j][3] += 0.5*v[3];
	vatom_thr[tid][j][4] += 0.5*v[4];
	vatom_thr[tid][j][5] += 0.5*v[5];
      }
    }
  }
}

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

/* ---------------------------------------------------------------------- */

// set loop range thread id, and force array offset for threaded runs.
double **ThrOMP::loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
				int inum, int nall, int nthreads)
{
#if defined(_OPENMP)
  tid = omp_get_thread_num();

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + inum/nthreads;
  ifrom = tid*idelta;
  ito   = ifrom + idelta;
  if (ito > inum)
    ito = inum;

  return f + nall*tid;
#else
  tid = 0;
  ifrom = 0;
  ito = inum;
  return f;
#endif
}

/* ---------------------------------------------------------------------- */

// reduce per thread data into the first part of the data
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on
// multi-dimensional arrays like force stored in this order
// x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with writing to the array .
void ThrOMP::data_reduce_thr(double *dall, int nall, int nthreads,
			     int ndim, int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    const int nvals = ndim*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nvals) ? nvals : (ifrom + idelta);

    for (int m = ifrom; m < ito; ++m) {
      for (int n = 1; n < nthreads; ++n) {
	dall[m] += dall[n*nvals + m];
	dall[n*nvals + m] = 0.0;
      }
    }
  }
#else
  // NOOP in non-threaded execution.
  return;
#endif
}

/* ---------------------------------------------------------------------- */

double ThrOMP::memory_usage_thr() 
{
  const int nthreads=lmp->comm->nthreads;

  double bytes = nthreads * (3 + 7) * sizeof(double);
  bytes += nthreads * maxeatom_thr * sizeof(double);
  bytes += nthreads * maxvatom_thr * 6 * sizeof(double);
  return bytes;
}
