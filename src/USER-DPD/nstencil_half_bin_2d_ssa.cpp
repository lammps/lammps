// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "nstencil_half_bin_2d_ssa.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfBin2dSSA::NStencilHalfBin2dSSA(LAMMPS *lmp) :
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

void NStencilHalfBin2dSSA::create()
{
  int i,j,pos = 0;
  nstencil_ssa[0] = 0; // redundant info, but saves a conditional

  // Include the centroid at the start.
  // It will be handled as part of Subphase 0.
  stencilxyz[pos][0] = 0;
  stencilxyz[pos][1] = 0;
  stencilxyz[pos][2] = 0;
  stencil[pos++] = 0;

  // Subphase 0: upper right front bins (red)
  for (j = 0; j <= sy; j++)
    for (i = 0; i <= sx; i++)
      if (j > 0 || i > 0) // skip the centroid
        if (bin_distance(i,j,0) < cutneighmaxsq) {
          stencilxyz[pos][0] = i;
          stencilxyz[pos][1] = j;
          stencilxyz[pos][2] = 0;
          stencil[pos++] = j*mbinx + i;
        }

  nstencil_ssa[1] = pos;
  // Subphase 1: upper left front bins (light blue)
  for (j = 1; j <= sy; j++)
    for (i = -sx; i < 0; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq) {
        stencilxyz[pos][0] = i;
        stencilxyz[pos][1] = j;
        stencilxyz[pos][2] = 0;
        stencil[pos++] = j*mbinx + i;
      }

  nstencil_ssa[2] = pos;
  // Subphase 2: lower right front bins (yellow)

  nstencil_ssa[3] = pos;
  // Subphase 3: lower left front bins (blue)

  nstencil_ssa[4] = pos; // record end of half stencil
  // Now include additional bins for AIR ghosts, and impure-to-pure locals
  // Subphase 4: upper right back bins (pink)

  // nstencil_ssa[5] = pos;
  // Subphase 5: upper left back bins (light green)

  // nstencil_ssa[6] = pos;
  // Subphase 6: lower right back bins (white)
  for (j = -sy; j < 0; j++)
    for (i = 0; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq) {
        stencilxyz[pos][0] = i;
        stencilxyz[pos][1] = j;
        stencilxyz[pos][2] = 0;
        stencil[pos++] = j*mbinx + i;
      }

  // nstencil_ssa[7] = pos;
  // Subphase 7: lower left back bins (purple)
  for (j = -sy; j <= 0; j++)
    for (i = -sx; i < 0; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq) {
        stencilxyz[pos][0] = i;
        stencilxyz[pos][1] = j;
        stencilxyz[pos][2] = 0;
        stencil[pos++] = j*mbinx + i;
      }
  // nstencil_ssa[8] = pos;

  nstencil = pos; // record where full stencil ends
}
