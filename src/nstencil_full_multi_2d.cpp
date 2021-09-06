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

#include "nstencil_full_multi_2d.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilFullMulti2d::NStencilFullMulti2d(LAMMPS *lmp) : NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilFullMulti2d::set_stencil_properties()
{
  int n = ncollections;
  int i, j;

  // Always look up neighbor using full stencil and neighbor's bin

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      flag_half_multi[i][j] = 0;
      flag_skip_multi[i][j] = 0;
      bin_collection_multi[i][j] = j;
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullMulti2d::create()
{
  int icollection, jcollection, bin_collection, i, j, ns;
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

      mbinx = stencil_mbinx_multi[icollection][jcollection];
      mbiny = stencil_mbiny_multi[icollection][jcollection];

      bin_collection = bin_collection_multi[icollection][jcollection];

      cutsq = cutcollectionsq[icollection][jcollection];

      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++)
          if (bin_distance_multi(i,j,0,bin_collection) < cutsq)
                stencil_multi[icollection][jcollection][ns++] = j*mbinx + i;

      nstencil_multi[icollection][jcollection] = ns;
    }
  }
}
