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
   Contributing authors:
   James Larentzos and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "nstencil_half_bin_3d_newton_ssa.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfBin3dNewtonSSA::NStencilHalfBin3dNewtonSSA(LAMMPS *lmp) :
  NStencilSSA(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
   3d creates xyz stencil, 2d creates xy stencil
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
   additionally, includes the bins beyond nstencil that are needed
     to locate all the Active Interaction Region (AIR) ghosts for SSA
------------------------------------------------------------------------- */

void NStencilHalfBin3dNewtonSSA::create()
{
  int i,j,k,pos = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (k > 0 || j > 0 || (j == 0 && i > 0))
          if (bin_distance(i,j,k) < cutneighmaxsq)
            stencil[pos++] = k*mbiny*mbinx + j*mbinx + i;

  nstencil_half = pos; // record where normal half stencil ends

  // include additional bins for AIR ghosts only

  for (k = -sz; k < 0; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[pos++] = k*mbiny*mbinx + j*mbinx + i;

  // For k==0, make sure to skip already included bins

  k = 0;
  for (j = -sy; j <= 0; j++)
    for (i = -sx; i <= sx; i++) {
      if (j == 0 && i > 0) continue;
      if (bin_distance(i,j,k) < cutneighmaxsq)
        stencil[pos++] = k*mbiny*mbinx + j*mbinx + i;
    }

  nstencil = pos; // record where full stencil ends
}
