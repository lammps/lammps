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
   Contributing author : Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "domain.h"
#include "accelerator_omp.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void DomainOMP::pbc()
{
  int i;
  tagint idim,otherdims;
  double *lo,*hi,*period;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  tagint *image = atom->image;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0],lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i][1] < lo[1]) {
        x[i][1] += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[5];
          v[i][1] += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
        x[i][1] -= period[1];
        x[i][1] = MAX(x[i][1],lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[5];
          v[i][1] -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }

    if (zperiodic) {
      if (x[i][2] < lo[2]) {
        x[i][2] += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[4];
          v[i][1] += h_rate[3];
          v[i][2] += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i][2] >= hi[2]) {
        x[i][2] -= period[2];
        x[i][2] = MAX(x[i][2],lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[4];
          v[i][1] -= h_rate[3];
          v[i][2] -= h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void DomainOMP::lamda2x(int n)
{
  double * const * const x = atom->x;
  const int num = n;
  int i;

#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
  for (i = 0; i < num; i++) {
    x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
    x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
    x[i][2] = h[2]*x[i][2] + boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void DomainOMP::x2lamda(int n)
{
  double * const * const x = atom->x;
  const int num = n;
  int i;

#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
  for (i = 0; i < num; i++) {
    double delta0 = x[i][0] - boxlo[0];
    double delta1 = x[i][1] - boxlo[1];
    double delta2 = x[i][2] - boxlo[2];

    x[i][0] = h_inv[0]*delta0 + h_inv[5]*delta1 + h_inv[4]*delta2;
    x[i][1] = h_inv[1]*delta1 + h_inv[3]*delta2;
    x[i][2] = h_inv[2]*delta2;
  }
}

