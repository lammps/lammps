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

#include "nstencil_half_bin_2d_tri.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfBin2dTri::NStencilHalfBin2dTri(LAMMPS *lmp) : NStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfBin2dTri::create()
{
  int i, j;

  // for triclinic, need to use full stencil in all dims
  //   not a half stencil in y
  // b/c transforming orthog -> lambda -> orthog for ghost atoms
  //   with an added PBC offset can shift both coords by epsilon
  // thus for an I/J owned/ghost pair, the xy coords
  //   and bin assignments can be different on I proc vs J proc

  nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i, j, 0) < cutneighmaxsq)
	stencil[nstencil++] = j * mbinx + i;
}
