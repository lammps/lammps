/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "nstencil_ghost_bin.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
NStencilGhostBin<HALF, DIM_3D, TRI>::NStencilGhostBin(LAMMPS *lmp) : NStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
void NStencilGhostBin<HALF, DIM_3D, TRI>::create()
{
  int i, j, k;

  // For half stencils, only the upper plane is needed
  int sy_min = sy;
  int sz_min = sz;
  if (HALF && (!DIM_3D)) sy_min = 0;
  if (HALF && DIM_3D) sz_min = 0;

  nstencil = 0;

  // Half and ortho stencils include central bin first
  // This preserves the historical order of the neighbor list
  // as the old npair classes used to separately parse the central bin first
  if (HALF && (!TRI)) stencil[nstencil++] = 0;

  for (k = -sz_min; k <= sz; k++) {
    for (j = -sy_min; j <= sy; j++) {
      for (i = -sx; i <= sx; i++) {

        // Now only include "upper right" bins for half and ortho stencils
        if (HALF && (!DIM_3D) && (!TRI))
          if (! (j > 0 || (j == 0 && i > 0))) continue;
        if (HALF && DIM_3D && (!TRI))
          if (! (k > 0 || j > 0 || (j == 0 && i > 0))) continue;

        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencilxyz[nstencil][0] = i;
          stencilxyz[nstencil][1] = j;
          stencilxyz[nstencil][2] = k;
          stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
      }
    }
  }
}

namespace LAMMPS_NS {
template class NStencilGhostBin<0,0,0>;
template class NStencilGhostBin<0,1,0>;
template class NStencilGhostBin<1,0,0>;
template class NStencilGhostBin<1,0,1>;
template class NStencilGhostBin<1,1,0>;
template class NStencilGhostBin<1,1,1>;
}
