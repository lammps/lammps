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

#include "ns_half_bin_3d_no_newton.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighStencilHalfBin3dNoNewton::NeighStencilHalfBin3dNoNewton(LAMMPS *lmp) : 
  NeighStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NeighStencilHalfBin3dNoNewton::create()
{
  int i,j,k;

  nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}
