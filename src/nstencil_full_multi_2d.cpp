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
  int n = n_multi_groups;
  int i, j;
  
  // Always look up neighbor using full stencil and neighbor's bin

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      flag_half_multi[i][j] = 0;
      flag_skip_multi[i][j] = 0;
      bin_group_multi[i][j] = j;
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullMulti2d::create()
{
  int igroup, jgroup, bin_group, i, j, k, ns;
  int n = n_multi_groups;
  double cutsq;
  
  
  for (igroup = 0; igroup < n; igroup++) {
    for (jgroup = 0; jgroup < n; jgroup++) {
      if (flag_skip_multi[igroup][jgroup]) continue;
      
      ns = 0;
      
      sx = stencil_sx_multi[igroup][jgroup];
      sy = stencil_sy_multi[igroup][jgroup];
      
      mbinx = stencil_mbinx_multi[igroup][jgroup];
      mbiny = stencil_mbiny_multi[igroup][jgroup];
      
      bin_group = bin_group_multi[igroup][jgroup];
      
      cutsq = cutmultisq[igroup][jgroup];
      
      for (j = -sy; j <= sy; j++)
        for (i = -sx; i <= sx; i++)
          if (bin_distance_multi(i,j,0,bin_group) < cutsq)
	        stencil_multi[igroup][jgroup][ns++] = j*mbinx + i;
      
      nstencil_multi[igroup][jgroup] = ns;
    }
  }
}
