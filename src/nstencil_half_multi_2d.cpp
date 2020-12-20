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

#include "nstencil_half_multi_2d.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfMulti2d::NStencilHalfMulti2d(LAMMPS *lmp) :
  NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilHalfMulti2d::set_stencil_properties()
{
  int n = atom->ntypes;
  int i, j;

  // Cross types: use full stencil, looking one way through hierarchy
  // smaller -> larger => use full stencil in larger bin
  // larger -> smaller => no nstencil required
  // If cut offs are same, use half stencil

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if(cutneighsq[i][i] > cutneighsq[j][j]) continue;

      flag_skip_multi[i][j] = 0;
      
      if(cutneighsq[i][i] == cutneighsq[j][j]){
        flag_half_multi[i][j] = 1;
        bin_type_multi[i][j] = i;
      } else {
        flag_half_multi[i][j] = 0;
        bin_type_multi[i][j] = j;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfMulti2d::create()
{
  int itype, jtype, bin_type, i, j, ns;
  int n = atom->ntypes;
  double cutsq;
  
  
  for (itype = 1; itype <= n; itype++) {
    for (jtype = 1; jtype <= n; jtype++) {
      if (flag_skip_multi[itype][jtype]) continue;
      
      ns = 0;
      
      sx = stencil_sx_multi[itype][jtype];
      sy = stencil_sy_multi[itype][jtype];
      
      mbinx = stencil_mbinx_multi[itype][jtype];
      mbiny = stencil_mbiny_multi[itype][jtype];
      
      bin_type = bin_type_multi[itype][jtype];      
      
      cutsq = cutneighsq[itype][jtype];
      
      if (flag_half_multi[itype][jtype]) {
        for (j = 0; j <= sy; j++)
          for (i = -sx; i <= sx; i++)
            if (j > 0 || (j == 0 && i > 0)) { 
              if (bin_distance_multi(i,j,0,bin_type) < cutsq)
                  stencil_multi[itype][jtype][ns++] = j*mbinx + i;
	        }
      } else {
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
              if (bin_distance_multi(i,j,0,bin_type) < cutsq)
	            stencil_multi[itype][jtype][ns++] = j*mbinx + i;
      }
      
      nstencil_multi[itype][jtype] = ns;
    }
  }
}

