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

#include "nstencil_manual.h"

#include "domain.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilManual::NStencilManual(LAMMPS *lmp) : NStencil(lmp)
{
  half = 0;
  dim_3d = domain->dimension == 3;
  tri = domain->triclinic;
  cutoff_custom = 0.0;
}

/* ----------------------------------------------------------------------
   manually assign neighbor info (only supports bin)
------------------------------------------------------------------------- */

void NStencilManual::assign_neighbor_info(double cutoff)
{
  cutneighmax = cutoff;
  cutneighmaxsq = cutoff * cutoff;
  neighstyle = Neighbor::BIN;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilManual::create()
{
  int i, j, k;

  // For half stencils, only the upper plane is needed
  // for triclinic, need to use full stencil in all dims
  //   not a half stencil in y
  // b/c transforming orthog -> lambda -> orthog for ghost atoms
  //   with an added PBC offset can shift both coords by epsilon
  // thus for an I/J owned/ghost pair, the xy coords
  //   and bin assignments can be different on I proc vs J proc

  int sy_min = sy;
  int sz_min = sz;
  if ((!tri) && half && (!dim_3d)) sy_min = 0;
  if ((!tri) && half && dim_3d) sz_min = 0;

  nstencil = 0;

  // Half and ortho stencils include central bin first
  // This preserves the historical order of the neighbor list
  // as the old npair classes used to separately parse the central bin first
  if (half && (!tri)) stencil[nstencil++] = 0;

  for (k = -sz_min; k <= sz; k++) {
    for (j = -sy_min; j <= sy; j++) {
      for (i = -sx; i <= sx; i++) {

        // Now only include "upper right" bins for half and ortho stencils
        if (half && (!dim_3d) && (!tri))
          if (j <= 0 && (j != 0 || i <= 0)) continue;
        if (half && dim_3d && (!tri))
          if (k <= 0 && j <= 0 && (j != 0 || i <= 0)) continue;

        if (bin_distance(i, j, k) < cutneighmaxsq)
          stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
      }
    }
  }
}
