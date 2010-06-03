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

#include "math.h"
#include "dihedral_omp.h"
#include "atom.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   set dihedral contribution to Vdwl and Coulombic energy to 0.0
   DihedralCharmm will override this
------------------------------------------------------------------------- */

DihedralOMP::DihedralOMP(LAMMPS *lmp) : Dihedral(lmp)
{
  const int nthreads = comm->nthreads;
  energy_thr = (double *)memory->smalloc(nthreads*sizeof(double),
					   "pair:eng_dihed_thr");
  virial_thr = memory->create_2d_double_array(nthreads,6,"pair:virial_thr");
  maxeatom_thr = maxvatom_thr = 0;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

DihedralOMP::~DihedralOMP()
{
  memory->sfree(energy_thr);
  memory->destroy_2d_double_array(virial_thr);
  memory->destroy_2d_double_array(eatom_thr);
  memory->destroy_3d_double_array(vatom_thr);
}

/*------------------------------------------------------------------------- */

void DihedralOMP::ev_setup_thr(int eflag, int vflag)
{
  int i,n,t;
  const int nthreads = comm->nthreads;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memory->sfree(eatom);
    eatom = (double *) memory->smalloc(maxeatom*sizeof(double),"bond:eatom");
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy_2d_double_array(vatom);
    vatom = memory->create_2d_double_array(maxvatom,6,"bond:vatom");
  }

  // zero accumulators

for (t = 0; t < nthreads; ++t) {
    if (eflag_global) energy_thr[t] = 0.0;
    if (vflag_global) for (i = 0; i < 6; i++) virial_thr[t][i] = 0.0;
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton) n += atom->nghost;
      for (i = 0; i < n; i++) eatom_thr[t][i] = 0.0;
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      for (i = 0; i < n; i++) {
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
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
	  = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

void DihedralOMP::ev_tally_thr(int i1, int i2, int i3, int i4,
			int nlocal, int newton_bond,
			double edihedral, double *f1, double *f3, double *f4,
			double vb1x, double vb1y, double vb1z,
			double vb2x, double vb2y, double vb2z,
			double vb3x, double vb3y, double vb3z, int tid)
{
  double edihedralquarter,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy_thr[tid] += edihedral;
      else {
	edihedralquarter = 0.25*edihedral;
	if (i1 < nlocal) energy_thr[tid] += edihedralquarter;
	if (i2 < nlocal) energy_thr[tid] += edihedralquarter;
	if (i3 < nlocal) energy_thr[tid] += edihedralquarter;
	if (i4 < nlocal) energy_thr[tid] += edihedralquarter;
      }
    }
    if (eflag_atom) {
      edihedralquarter = 0.25*edihedral;
      if (newton_bond || i1 < nlocal) eatom_thr[tid][i1] += edihedralquarter;
      if (newton_bond || i2 < nlocal) eatom_thr[tid][i2] += edihedralquarter;
      if (newton_bond || i3 < nlocal) eatom_thr[tid][i3] += edihedralquarter;
      if (newton_bond || i4 < nlocal) eatom_thr[tid][i4] += edihedralquarter;
    }
  }

  if (vflag_either) {
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (vflag_global) {
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

    if (vflag_atom) {
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
   reduce the per thread accumulated E/V data into the canonical accumulators.
------------------------------------------------------------------------- */
void DihedralOMP::ev_reduce_thr()
{
  const int nthreads=comm->nthreads;

  for (int n = 0; n < nthreads; ++n) {
    energy += energy_thr[n];
    if (vflag_either) {
      virial[0] += virial_thr[n][0];
      virial[1] += virial_thr[n][1];
      virial[2] += virial_thr[n][2];
      virial[3] += virial_thr[n][3];
      virial[4] += virial_thr[n][4];
      virial[5] += virial_thr[n][5];
      if (vflag_atom) {
	for (int i = 0; i < atom->nmax; ++i) {
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
      for (int i = 0; i < atom->nmax; ++i) {
	eatom[i] += eatom_thr[n][i];
      }
    }
  }
}
/* ---------------------------------------------------------------------- */

double DihedralOMP::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  return bytes;
}

