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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair_omp.h"
#include "memory.h"

#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOMP::PairOMP(LAMMPS *lmp) : Pair(lmp)
{
  // for hybrid OpenMP/MPI we need multiple copies
  // of some accumulators to avoid race conditions
  const int nthreads = comm->nthreads;
  is_omp = 1;
  eng_vdwl_thr = (double *)memory->smalloc(nthreads*sizeof(double),
					   "pair:eng_vdwl_thr");
  eng_coul_thr = (double *)memory->smalloc(nthreads*sizeof(double),
					   "pair:eng_coul_thr");
  virial_thr = memory->create_2d_double_array(nthreads,6,"pair:virial_thr");
  
  maxeatom_thr = maxvatom_thr = 0;
  eatom_thr = NULL;
  vatom_thr = NULL;
}

/* ---------------------------------------------------------------------- */

PairOMP::~PairOMP()
{
  memory->sfree(eng_vdwl_thr);
  memory->sfree(eng_coul_thr);
  memory->destroy_2d_double_array(virial_thr);
  memory->destroy_2d_double_array(eatom_thr);
  memory->destroy_3d_double_array(vatom_thr);
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation. additional code for multi-threading
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void PairOMP::ev_setup_thr(int eflag, int vflag)
{
  int i,n,t;
  const int nthreads = comm->nthreads;

  // reallocate per-atom arrays if necessary
  if (eflag_atom && atom->nmax > maxeatom_thr) {
    maxeatom_thr = atom->nmax;
    memory->destroy_2d_double_array(eatom_thr);
    eatom_thr = memory->create_2d_double_array(nthreads,
					       maxeatom_thr,"pair:eatom_thr");
  }
  if (vflag_atom && atom->nmax > maxvatom_thr) {
    maxvatom_thr = atom->nmax;
    memory->destroy_3d_double_array(vatom_thr);
    vatom_thr = memory->create_3d_double_array(nthreads,
					       maxvatom_thr,6,"pair:vatom_thr");
  }
  
  // zero per thread accumulators
  // use force->newton instead of newton_pair
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  const int ntotal = (force->newton) ? 
    (atom->nlocal + atom->nghost) : atom->nlocal;
  for (t = 0; t < nthreads; ++t) {
    if (eflag_global) eng_vdwl_thr[t] = eng_coul_thr[t] = 0.0;
    if (vflag_global) for (i = 0; i < 6; ++i) virial_thr[t][i] = 0.0;
    if (eflag_atom) {
      for (i = 0; i < ntotal; ++i) eatom_thr[t][i] = 0.0;
    }
    if (vflag_atom) {
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
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per thread global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void PairOMP::ev_tally_thr(int i, int j, int nlocal, int newton_pair,
			   double evdwl, double ecoul, double fpair,
			   double delx, double dely, double delz, int tid)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
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
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom_thr[tid][i] += epairhalf;
      if (newton_pair || j < nlocal) eatom_thr[tid][j] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (vflag_global) {
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

    if (vflag_atom) {
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

void PairOMP::ev_tally_xyz_thr(int i, int j, int nlocal, int newton_pair,
			    double evdwl, double ecoul,
			    double fx, double fy, double fz,
			    double delx, double dely, double delz, int tid)
{
  double evdwlhalf,ecoulhalf,epairhalf,v[6];
  
  if (eflag_either) {
    if (eflag_global) {
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
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) eatom_thr[tid][i] += epairhalf;
      if (newton_pair || j < nlocal) eatom_thr[tid][j] += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    if (vflag_global) {
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

    if (vflag_atom) {
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
   called by SW potential, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

void PairOMP::ev_tally3_thr(int i, int j, int k, double evdwl, double ecoul,
			    double *fj, double *fk, double *drji, double *drki, int tid)
{
  double epairthird,v[6];

  if (eflag_either) {
    if (eflag_global) {
      eng_vdwl_thr[tid] += evdwl;
      eng_coul_thr[tid] += ecoul;
    }
    if (eflag_atom) {
      epairthird = THIRD * (evdwl + ecoul);
      eatom_thr[tid][i] += epairthird;
      eatom_thr[tid][j] += epairthird;
      eatom_thr[tid][k] += epairthird;
    }
  }

  if (vflag_atom) {
    v[0] = THIRD * (drji[0]*fj[0] + drki[0]*fk[0]);
    v[1] = THIRD * (drji[1]*fj[1] + drki[1]*fk[1]);
    v[2] = THIRD * (drji[2]*fj[2] + drki[2]*fk[2]);
    v[3] = THIRD * (drji[0]*fj[1] + drki[0]*fk[1]);
    v[4] = THIRD * (drji[0]*fj[2] + drki[0]*fk[2]);
    v[5] = THIRD * (drji[1]*fj[2] + drki[1]*fk[2]);

    vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1];
    vatom_thr[tid][i][2] += v[2]; vatom_thr[tid][i][3] += v[3];
    vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
    vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1];
    vatom_thr[tid][j][2] += v[2]; vatom_thr[tid][j][3] += v[3];
    vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
    vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1];
    vatom_thr[tid][k][2] += v[2]; vatom_thr[tid][k][3] += v[3];
    vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by AIREBO potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void PairOMP::ev_tally4_thr(int i, int j, int k, int m, double evdwl,
			    double *fi, double *fj, double *fk,
			    double *drim, double *drjm, double *drkm, int tid)
{
  double epairfourth,v[6];

  if (eflag_either) {
    if (eflag_global) eng_vdwl_thr[tid] += evdwl;
    if (eflag_atom) {
      epairfourth = 0.25 * evdwl;
      eatom_thr[tid][i] += epairfourth;
      eatom_thr[tid][j] += epairfourth;
      eatom_thr[tid][k] += epairfourth;
      eatom_thr[tid][m] += epairfourth;
    }
  }

  if (vflag_atom) {
    v[0] = 0.25 * (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
    v[1] = 0.25 * (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
    v[2] = 0.25 * (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
    v[3] = 0.25 * (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
    v[4] = 0.25 * (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
    v[5] = 0.25 * (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);
    
    vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1];
    vatom_thr[tid][i][2] += v[2]; vatom_thr[tid][i][3] += v[3];
    vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
    vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1];
    vatom_thr[tid][j][2] += v[2]; vatom_thr[tid][j][3] += v[3];
    vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
    vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1];
    vatom_thr[tid][k][2] += v[2]; vatom_thr[tid][k][3] += v[3];
    vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
    vatom_thr[tid][m][0] += v[0]; vatom_thr[tid][m][1] += v[1];
    vatom_thr[tid][m][2] += v[2]; vatom_thr[tid][m][3] += v[3];
    vatom_thr[tid][m][4] += v[4]; vatom_thr[tid][m][5] += v[5];
  }
}

/* ----------------------------------------------------------------------
   tally ecoul and virial into each of n atoms in list
   called by TIP4P potential, newton_pair is always on
   changes v values by dividing by n
 ------------------------------------------------------------------------- */

void PairOMP::ev_tally_list_thr(int n, int *list, double ecoul, double *v, int tid)
{
  int i,j;

  if (eflag_either) {
    if (eflag_global) eng_coul_thr[tid] += ecoul;
    if (eflag_atom) {
      double epairatom = ecoul/n;
      for (i = 0; i < n; i++) eatom_thr[tid][list[i]] += epairatom;
    }
  }

  if (vflag_either) {
    if (vflag_global) {
      virial_thr[tid][0] += v[0];
      virial_thr[tid][1] += v[1];
      virial_thr[tid][2] += v[2];
      virial_thr[tid][3] += v[3];
      virial_thr[tid][4] += v[4];
      virial_thr[tid][5] += v[5];
    }

    if (vflag_atom) {
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
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void PairOMP::v_tally2_thr(int i, int j, double fpair, double *drij, int tid)
{
  double v[6];

  v[0] = 0.5 * drij[0]*drij[0]*fpair;
  v[1] = 0.5 * drij[1]*drij[1]*fpair;
  v[2] = 0.5 * drij[2]*drij[2]*fpair;
  v[3] = 0.5 * drij[0]*drij[1]*fpair;
  v[4] = 0.5 * drij[0]*drij[2]*fpair;
  v[5] = 0.5 * drij[1]*drij[2]*fpair;

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1];
  vatom_thr[tid][i][2] += v[2]; vatom_thr[tid][i][3] += v[3];
  vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1];
  vatom_thr[tid][j][2] += v[2]; vatom_thr[tid][j][3] += v[3];
  vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void PairOMP::v_tally3_thr(int i, int j, int k, double *fi, double *fj,
                           double *drik, double *drjk, int tid)
{
  double v[6];

  v[0] = THIRD * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = THIRD * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = THIRD * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = THIRD * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = THIRD * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = THIRD * (drik[1]*fi[2] + drjk[1]*fj[2]);

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1];
  vatom_thr[tid][i][2] += v[2]; vatom_thr[tid][i][3] += v[3];
  vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1];
  vatom_thr[tid][j][2] += v[2]; vatom_thr[tid][j][3] += v[3];
  vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
  vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1];
  vatom_thr[tid][k][2] += v[2]; vatom_thr[tid][k][3] += v[3];
  vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
------------------------------------------------------------------------- */

void PairOMP::v_tally4_thr(int i, int j, int k, int m,
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

  vatom_thr[tid][i][0] += v[0]; vatom_thr[tid][i][1] += v[1];
  vatom_thr[tid][i][2] += v[2]; vatom_thr[tid][i][3] += v[3];
  vatom_thr[tid][i][4] += v[4]; vatom_thr[tid][i][5] += v[5];
  vatom_thr[tid][j][0] += v[0]; vatom_thr[tid][j][1] += v[1];
  vatom_thr[tid][j][2] += v[2]; vatom_thr[tid][j][3] += v[3];
  vatom_thr[tid][j][4] += v[4]; vatom_thr[tid][j][5] += v[5];
  vatom_thr[tid][k][0] += v[0]; vatom_thr[tid][k][1] += v[1];
  vatom_thr[tid][k][2] += v[2]; vatom_thr[tid][k][3] += v[3];
  vatom_thr[tid][k][4] += v[4]; vatom_thr[tid][k][5] += v[5];
  vatom_thr[tid][m][0] += v[0]; vatom_thr[tid][m][1] += v[1];
  vatom_thr[tid][m][2] += v[2]; vatom_thr[tid][m][3] += v[3];
  vatom_thr[tid][m][4] += v[4]; vatom_thr[tid][m][5] += v[5];
}

/* ----------------------------------------------------------------------
   reduce the per thread accumulated E/V data into the canonical accumulators.
------------------------------------------------------------------------- */
void PairOMP::ev_reduce_thr()
{
  const int nthreads=comm->nthreads;
  const int ntotal = (force->newton) ? 
    (atom->nlocal + atom->nghost) : atom->nlocal;

  for (int n = 0; n < nthreads; ++n) {
    eng_vdwl += eng_vdwl_thr[n];
    eng_coul += eng_coul_thr[n];
    if (vflag_either) {
      virial[0] += virial_thr[n][0];
      virial[1] += virial_thr[n][1];
      virial[2] += virial_thr[n][2];
      virial[3] += virial_thr[n][3];
      virial[4] += virial_thr[n][4];
      virial[5] += virial_thr[n][5];
      if (vflag_atom) {
	for (int i = 0; i < ntotal; ++i) {
	  vatom[i][0] += vatom_thr[n][i][0];
	  vatom[i][1] += vatom_thr[n][i][1];
	  vatom[i][2] += vatom_thr[n][i][2];
	  vatom[i][3] += vatom_thr[n][i][3];
	  vatom[i][4] += vatom_thr[n][i][4];
	  vatom[i][5] += vatom_thr[n][i][5];
	}
      }
    }
    if (eflag_atom) {
      for (int i = 0; i < ntotal; ++i) {
	eatom[i] += eatom_thr[n][i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void *PairOMP::extract(char *str, int &dim) {
    if (strcmp(str,"omp") == 0) return (void *) &is_omp;
    return NULL;
}

/* ---------------------------------------------------------------------- */

double PairOMP::memory_usage()
{
  const int nthreads=comm->nthreads;
  
  double bytes = Pair::memory_usage();

  bytes += nthreads * (2 + 7) * sizeof(double);
  bytes += nthreads * maxeatom_thr * sizeof(double);
  bytes += nthreads * maxvatom_thr * 6 * sizeof(double);

  return bytes;
}
