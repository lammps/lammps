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

#include "nstencil_bin_intel.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
NStencilBinIntel<HALF, DIM_3D, TRI>::NStencilBinIntel(LAMMPS *lmp) : NStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
void NStencilBinIntel<HALF, DIM_3D, TRI>::create()
{
  int i, j, k;

  // For half stencils, only the upper plane is needed
  int sy_min = sy;
  int sz_min = sz;
  if ((!TRI) && HALF && (!DIM_3D)) sy_min = 0;
  if ((!TRI) && HALF && DIM_3D) sz_min = 0;

  nstencil = 0;

  // For Intel, half and ortho stencils do not include central bin
  // as, historically, this was never included in a stencil.
  // Non-Intel npair classes were updated to account for this change,
  // but the Intel npair classes have not yet been updated
  // if (HALF && (!TRI)) stencil[nstencil++] = 0;

  for (k = -sz_min; k <= sz; k++) {
    for (j = -sy_min; j <= sy; j++) {
      for (i = -sx; i <= sx; i++) {

        // Now only include "upper right" bins for half and ortho stencils
        if (HALF && (!DIM_3D) && (!TRI))
          if (! (j > 0 || (j == 0 && i > 0))) continue;
        if (HALF && DIM_3D && (!TRI))
          if (! (k > 0 || j > 0 || (j == 0 && i > 0))) continue;

        if (bin_distance(i, j, k) < cutneighmaxsq)
          stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
      }
    }
  }
}

namespace LAMMPS_NS {
template class NStencilBinIntel<0,0,0>;
template class NStencilBinIntel<0,1,0>;
template class NStencilBinIntel<1,0,0>;
template class NStencilBinIntel<1,0,1>;
template class NStencilBinIntel<1,1,0>;
template class NStencilBinIntel<1,1,1>;
}
