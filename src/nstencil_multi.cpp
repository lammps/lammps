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

#include "nstencil_multi.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
NStencilMulti<HALF, DIM_3D, TRI>::NStencilMulti(LAMMPS *lmp) : NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
void NStencilMulti<HALF, DIM_3D, TRI>::set_stencil_properties()
{
  int n = ncollections;
  int i, j;

  // FULL
  // Always look up neighbor using full stencil and neighbor's bin
  // Stencil cutoff set by i-j cutoff

  // HALF
  // Cross collections: use full stencil, looking one way through hierarchy
  // smaller -> larger => use full stencil in larger bin
  // larger -> smaller => no nstencil required
  // If cut offs are same, use half stencil
  // If triclinic, need full stencil

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (HALF)
        if (cutcollectionsq[i][i] > cutcollectionsq[j][j]) continue;

      flag_skip_multi[i][j] = false;
      flag_half_multi[i][j] = false;
      flag_same_multi[i][j] = false;
      bin_collection_multi[i][j] = j;

      if (HALF) {
        if (cutcollectionsq[i][i] == cutcollectionsq[j][j]) {
          if (!TRI) flag_half_multi[i][j] = true;
          flag_same_multi[i][j] = true;
          bin_collection_multi[i][j] = i;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
void NStencilMulti<HALF, DIM_3D, TRI>::create()
{
  int icollection, jcollection, bin_collection, i, j, k, ns, half_flag;
  int n = ncollections;
  double cutsq;

  for (icollection = 0; icollection < n; icollection++) {
    for (jcollection = 0; jcollection < n; jcollection++) {
      if (flag_skip_multi[icollection][jcollection]) {
        nstencil_multi[icollection][jcollection] = 0;
        continue;
      }

      ns = 0;

      sx = stencil_sx_multi[icollection][jcollection];
      sy = stencil_sy_multi[icollection][jcollection];
      sz = stencil_sz_multi[icollection][jcollection];

      mbinx = stencil_mbinx_multi[icollection][jcollection];
      mbiny = stencil_mbiny_multi[icollection][jcollection];
      mbinz = stencil_mbinz_multi[icollection][jcollection];

      bin_collection = bin_collection_multi[icollection][jcollection];
      cutsq = cutcollectionsq[icollection][jcollection];
      half_flag = flag_half_multi[icollection][jcollection];

      // Half and ortho stencils include central bin first
      // This preserves the historical order of the neighbor list
      // as the old npair classes used to separately parse the central bin first
      // This !TRI condition (and the one below) are now unnecessary
      // since triclinic only uses full stencils - kept the flags for clarity
      if (HALF && (!TRI))
        if (half_flag) stencil_multi[icollection][jcollection][ns++] = 0;

      // For half stencils, only the upper plane is needed
      int sy_min = sy;
      int sz_min = sz;
      if (HALF) {
        if (half_flag && (!DIM_3D)) sy_min = 0;
        if (half_flag && DIM_3D) sz_min = 0;
      }

      for (k = -sz_min; k <= sz; k++) {
        for (j = -sy_min; j <= sy; j++) {
          for (i = -sx; i <= sx; i++) {
            // Now only include "upper right" bins for half and ortho stencils
            if (HALF && (!TRI)) {
              if (half_flag) {
                if (DIM_3D) {
                  if (k <= 0 && j <= 0 && (j != 0 || i <= 0)) continue;
                } else {
                  if (j <= 0 && (j != 0 || i <= 0)) continue;
                }
              }
            }
            if (bin_distance_multi(i, j, k, bin_collection) < cutsq)
              stencil_multi[icollection][jcollection][ns++] = k * mbiny * mbinx + j * mbinx + i;
          }
        }
      }

      nstencil_multi[icollection][jcollection] = ns;
    }
  }
}

namespace LAMMPS_NS {
template class NStencilMulti<0,0,0>;
template class NStencilMulti<0,1,0>;
template class NStencilMulti<1,0,0>;
template class NStencilMulti<1,0,1>;
template class NStencilMulti<1,1,0>;
template class NStencilMulti<1,1,1>;
}
