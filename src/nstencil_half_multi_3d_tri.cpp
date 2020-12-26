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

#include "nstencil_half_multi_3d_tri.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfMulti3dTri::NStencilHalfMulti3dTri(LAMMPS *lmp) :
  NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilHalfMulti3dTri::set_stencil_properties()
{
  int n = n_multi_groups;
  int i, j;
  
  // Cross groups: use full stencil, looking one way through hierarchy
  // smaller -> larger => use full stencil in larger bin
  // larger -> smaller => no nstencil required
  // If cut offs are same, use half stencil

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if(cutmultisq[i][i] > cutmultisq[j][j]) continue;

      flag_skip_multi[i][j] = 0;
      
      if(cutmultisq[i][i] == cutmultisq[j][j]){
        flag_half_multi[i][j] = 1;
        bin_group_multi[i][j] = i;
      } else {
        flag_half_multi[i][j] = 0;
        bin_group_multi[i][j] = j;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfMulti3dTri::create()
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
      sz = stencil_sz_multi[igroup][jgroup];
      
      mbinx = stencil_mbinx_multi[igroup][jgroup];
      mbiny = stencil_mbiny_multi[igroup][jgroup];
      mbinz = stencil_mbinz_multi[igroup][jgroup];
      
      bin_group = bin_group_multi[igroup][jgroup];
      
      cutsq = cutmultisq[igroup][jgroup];
      
      if (flag_half_multi[igroup][jgroup]) {
        for (k = 0; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
              if (bin_distance_multi(i,j,k,bin_group) < cutsq)
	            stencil_multi[igroup][jgroup][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      } else {
        for (k = -sz; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
	          if (bin_distance_multi(i,j,k,bin_group) < cutsq)
	            stencil_multi[igroup][jgroup][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      }
      
      nstencil_multi[igroup][jgroup] = ns;
    }
  }
}
