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

#include "accelerator_omp.h"
#include "atom.h"

using namespace LAMMPS_NS;

typedef struct { double x,y,z; } dbl3_t;

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
  dbl3_t * _noalias const x = (dbl3_t *)&atom->x[0][0];
  dbl3_t * _noalias const v = (dbl3_t *)&atom->v[0][0];
  const double * _noalias const lo     = (triclinic == 0) ? boxlo : boxlo_lamda;
  const double * _noalias const hi     = (triclinic == 0) ? boxhi : boxhi_lamda;
  const double * _noalias const period = (triclinic == 0) ? prd   : prd_lamda;
  const int * _noalias const mask  = atom->mask;
  imageint    * _noalias const image = atom->image;
  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    imageint idim,otherdims;

    if (xperiodic) {
      if (x[i].x < lo[0]) {
        x[i].x += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i].x += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i].x >= hi[0]) {
        x[i].x -= period[0];
        x[i].x = MAX(x[i].x,lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i].x -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i].y < lo[1]) {
        x[i].y += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i].x += h_rate[5];
          v[i].y += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i].y >= hi[1]) {
        x[i].y -= period[1];
        x[i].y = MAX(x[i].y,lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i].x -= h_rate[5];
          v[i].y -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }

    if (zperiodic) {
      if (x[i].z < lo[2]) {
        x[i].z += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i].x += h_rate[4];
          v[i].y += h_rate[3];
          v[i].z += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i].z >= hi[2]) {
        x[i].z -= period[2];
        x[i].z = MAX(x[i].z,lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i].x -= h_rate[4];
          v[i].y -= h_rate[3];
          v[i].z -= h_rate[2];
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
  dbl3_t * _noalias const x = (dbl3_t *)&atom->x[0][0];
  const int num = n;

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int i = 0; i < num; i++) {
    x[i].x = h[0]*x[i].x + h[5]*x[i].y + h[4]*x[i].z + boxlo[0];
    x[i].y = h[1]*x[i].y + h[3]*x[i].z + boxlo[1];
    x[i].z = h[2]*x[i].z + boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void DomainOMP::x2lamda(int n)
{
  dbl3_t * _noalias const x = (dbl3_t *)&atom->x[0][0];
  const int num = n;

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int i = 0; i < num; i++) {
    double delta0 = x[i].x - boxlo[0];
    double delta1 = x[i].y - boxlo[1];
    double delta2 = x[i].z - boxlo[2];

    x[i].x = h_inv[0]*delta0 + h_inv[5]*delta1 + h_inv[4]*delta2;
    x[i].y = h_inv[1]*delta1 + h_inv[3]*delta2;
    x[i].z = h_inv[2]*delta2;
  }
}

