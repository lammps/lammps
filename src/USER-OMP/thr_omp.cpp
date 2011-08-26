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

using namespace LAMMPS_NS;

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

void ThrOMP::ev_zero_acc_thr(int ntotal, int eflag_global, int vflag_global,
			     int eflag_atom, int vflag_atom, int nthreads)
{
  int t,i;
  
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

  // zero per thread accumulators
  ev_zero_acc_thr(ntotal, dihed->eflag_global, dihed->vflag_global,
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

  // zero per thread accumulators
  ev_zero_acc_thr(ntotal, pair->eflag_global, pair->vflag_global,
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

/* ---------------------------------------------------------------------- */

// set loop range thread id, and force array offset for threaded runs.
double **ThrOMP::loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
				int inum, int nall, int nthreads)
{
#if defined(_OPENMP)
  if (nthreads > 1) {
    tid = omp_get_thread_num();

    // each thread works on a fixed chunk of atoms.
    const int idelta = 1 + inum/nthreads;
    ifrom = tid*idelta;
    ito   = ifrom + idelta;
    if (ito > inum)
      ito = inum;

    return f + nall*tid;

  } else {
#endif
    tid = 0;
    ifrom = 0;
    ito = inum;
    return f;
#if defined(_OPENMP)
  }
#endif
}

/* ---------------------------------------------------------------------- */

// reduce per thread forces into the first part of the force
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on the
// forces arrays stored in this order x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with computing forces.
void ThrOMP::force_reduce_thr(double *fall, int nall,
			      int nthreads, int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    double *f;
    const int idelta = 1 + nall/nthreads;
    const int ifrom = 3*tid*idelta;
    const int ito   = 3*(((ifrom + idelta) > nall) ? nall : (ifrom + idelta));

    for (int n = 1; n < nthreads; ++n) {
      const int toffs = 3*n*nall;
      f = fall + toffs;
      for (int m = ifrom; m < ito; ++m) {
	fall[m] += f[m];
	f[m] = 0.0;
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
