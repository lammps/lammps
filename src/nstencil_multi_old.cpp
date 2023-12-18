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

#include "nstencil_multi_old.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
NStencilMultiOld<HALF, DIM_3D, TRI>::NStencilMultiOld(LAMMPS *lmp) : NStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

template<int HALF, int DIM_3D, int TRI>
void NStencilMultiOld<HALF, DIM_3D, TRI>::create()
{
  int i, j, k, n, itype;
  double rsq, typesq;
  int *s;
  double *distsq;
  int ntypes = atom->ntypes;

  // For half stencils, only the upper plane is needed
  // for triclinic, need to use full stencil in all dims
  //   not a half stencil in y
  // b/c transforming orthog -> lambda -> orthog for ghost atoms
  //   with an added PBC offset can shift both coords by epsilon
  // thus for an I/J owned/ghost pair, the xy coords
  //   and bin assignments can be different on I proc vs J proc

  int sy_min = sy;
  int sz_min = sz;
  if ((!TRI) && HALF && (!DIM_3D)) sy_min = 0;
  if ((!TRI) && HALF && DIM_3D) sz_min = 0;

  for (itype = 1; itype <= ntypes; itype++) {

    typesq = cuttypesq[itype];
    s = stencil_multi_old[itype];
    distsq = distsq_multi_old[itype];
    n = 0;

    // Half and ortho stencils include central bin first
    // This preserves the historical order of the neighbor list
    // as the old npair classes used to separately parse the central bin first
    if (HALF && (!TRI)) s[n++] = 0;

    for (k = -sz_min; k <= sz; k++) {
      for (j = -sy_min; j <= sy; j++) {
        for (i = -sx; i <= sx; i++) {

          // Now only include "upper right" bins for half and ortho stencils
          if (HALF && (!DIM_3D) && (!TRI))
            if (! (j > 0 || (j == 0 && i > 0))) continue;
          if (HALF && DIM_3D && (!TRI))
            if (! (k > 0 || j > 0 || (j == 0 && i > 0))) continue;

          rsq = bin_distance(i, j, k);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = k * mbiny * mbinx + j * mbinx + i;
          }
        }
      }
    }
    nstencil_multi_old[itype] = n;
  }
}

namespace LAMMPS_NS {
template class NStencilMultiOld<0,0,0>;
template class NStencilMultiOld<0,1,0>;
template class NStencilMultiOld<1,0,0>;
template class NStencilMultiOld<1,0,1>;
template class NStencilMultiOld<1,1,0>;
template class NStencilMultiOld<1,1,1>;
}
