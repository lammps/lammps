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

#include "nstencil_half_multi2_3d_tri.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfMulti23dTri::NStencilHalfMulti23dTri(LAMMPS *lmp) :
  NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilHalfMulti23dTri::set_stencil_properties()
{
  int n = atom->ntypes;
  int i, j;
  
  // like -> like => use half stencil
  for (i = 1; i <= n; i++) {
    stencil_half[i][i] = 1;
    stencil_skip[i][i] = 0;
    stencil_bin_type[i][i] = i;
    stencil_cut[i][i] = sqrt(cutneighsq[i][i]);
  }

  // Cross types: use full stencil, looking one way through hierarchy
  // smaller -> larger => use full stencil in larger bin
  // larger -> smaller => no nstecil required
  // If cut offs are same, use existing type-type stencil

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if(i == j) continue;
      if(cuttypesq[i] > cuttypesq[j]) continue;

      stencil_skip[i][j] = 0;
      stencil_cut[i][j] = sqrt(cutneighsq[i][j]);          
      
      if(cuttypesq[i] == cuttypesq[j]){
        stencil_half[i][j] = 1;
        stencil_bin_type[i][j] = i;
      } else {
        stencil_half[i][j] = 0;
        stencil_bin_type[i][j] = j;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfMulti23dTri::create()
{
  int itype, jtype, i, j, k, ns;
  int n = atom->ntypes;
  double cutsq;
  
  
  for (itype = 1; itype <= n; itype++) {
    for (jtype = 1; jtype <= n; jtype++) {
      if (stencil_skip[itype][jtype]) continue;
      
      ns = 0;
      
      sx = sx_multi2[itype][jtype];
      sy = sy_multi2[itype][jtype];
      sz = sz_multi2[itype][jtype];
      
      mbinx = mbinx_multi2[itype][jtype];
      mbiny = mbiny_multi2[itype][jtype];
      mbinz = mbinz_multi2[itype][jtype];
      
      cutsq = stencil_cut[itype][jtype];
      
      if (stencil_half[itype][jtype]) {
        for (k = 0; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
              if (bin_distance(i,j,k) < cutsq)
	            stencil_multi2[itype][jtype][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      } else {
        for (k = -sz; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
	          if (bin_distance(i,j,k) < cutsq)
	            stencil_multi2[itype][jtype][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      }
      
      nstencil_multi2[itype][jtype] = ns;
    }
  }
}
