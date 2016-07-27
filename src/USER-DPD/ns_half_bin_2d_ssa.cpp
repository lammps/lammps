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

#include "ns_half_bin_2d_ssa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighStencilHalfBin2dSSA::NeighStencilHalfBin2dSSA(LAMMPS *lmp) : 
  NeighStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
   3d creates xyz stencil, 2d creates xy stencil
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
   Additionally, includes the bins beyond nstencil that are needed
   to locate all the Active Interaction Region (AIR) ghosts for SSA
------------------------------------------------------------------------- */

void NeighStencilHalfBin2dSSA::create()
{
  int i,j;

  if (naux_nstencil < 1) {
    naux_nstencil = 1;
    if (aux_nstencil) memory->destroy(aux_nstencil);
    memory->create(aux_nstencil,naux_nstencil,"neighstencil:aux_nstencil");
  }

  nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (j > 0 || (j == 0 && i > 0))
        if (bin_distance(i,j,0) < cutneighmaxsq)
          stencil[nstencil++] = j*mbinx + i;

  aux_nstencil[0] = nstencil; // record where non-ghost half stencil ends

  // include additional bins for AIR ghosts only

  for (j = -sy; j <= 0; j++)
    for (i = -sx; i <= sx; i++) {
      if (j == 0 && i > 0) continue;
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;
    }
}
