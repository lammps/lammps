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

#include <math.h>
#include "npair.h"
#include "neighbor.h"
#include "nbin.h"
#include "nstencil.h"
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPair::NPair(LAMMPS *lmp) : Pointers(lmp)
{
  last_build = -1;
  last_copy_bin_setup = last_copy_bin = last_copy_stencil = -1;

  molecular = atom->molecular;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
------------------------------------------------------------------------- */

void NPair::copy_neighbor_info()
{
  // general params

  includegroup = neighbor->includegroup;
  exclude = neighbor->exclude;
  skin = neighbor->skin;
  cutneighsq = neighbor->cutneighsq;
  cutneighghostsq = neighbor->cutneighghostsq;
  cut_inner_sq = neighbor->cut_inner_sq;
  cut_middle_sq = neighbor->cut_middle_sq;
  cut_middle_inside_sq = neighbor->cut_middle_inside_sq;
  zeroes = neighbor->zeroes;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // exclusion info

  nex_type = neighbor->nex_type;
  ex1_type = neighbor->ex1_type;
  ex2_type = neighbor->ex2_type;
  ex_type = neighbor->ex_type;

  nex_group = neighbor->nex_group;
  ex1_group = neighbor->ex1_group;
  ex2_group = neighbor->ex2_group;
  ex1_bit = neighbor->ex1_bit;
  ex2_bit = neighbor->ex2_bit;

  nex_mol = neighbor->nex_mol;
  ex_mol_group = neighbor->ex_mol_group;
  ex_mol_bit = neighbor->ex_mol_bit;

  // special info

  special_flag = neighbor->special_flag;
}

/* ----------------------------------------------------------------------
   copy bin geometry info from NBin class to this build class
------------------------------------------------------------------------- */

void NPair::copy_bin_setup_info()
{
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;
 
  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;
}

/* ----------------------------------------------------------------------
   copy per-atom and per-bin vectors from NBin class to this build class
------------------------------------------------------------------------- */

void NPair::copy_bin_info()
{
  bins = nb->bins;
  binhead = nb->binhead;
}

/* ----------------------------------------------------------------------
   copy needed info from NStencil class to this build class
------------------------------------------------------------------------- */

void NPair::copy_stencil_info()
{
  nstencil = ns->nstencil;
  stencil = ns->stencil;
  stencilxyz = ns->stencilxyz;
  nstencil_multi = ns->nstencil_multi;
  stencil_multi = ns->stencil_multi;
  distsq_multi = ns->distsq_multi;
}

/* ----------------------------------------------------------------------
   copy needed info from NStencil class to this build class
------------------------------------------------------------------------- */

void NPair::build_setup()
{
  if (nb && last_copy_bin_setup < nb->last_setup) {
    copy_bin_setup_info();
    last_copy_bin_setup = update->ntimestep;
  }
  if (nb && last_copy_bin < nb->last_bin_memory) {
    copy_bin_info();
    last_copy_bin = update->ntimestep;
  }
  if (ns && last_copy_stencil < ns->last_create) {
    copy_stencil_info();
    last_copy_stencil = update->ntimestep;
  }

  last_build = update->ntimestep;
}

/* ----------------------------------------------------------------------
   test if atom pair i,j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
------------------------------------------------------------------------- */

int NPair::exclusion(int i, int j, int itype, int jtype,
                          int *mask, tagint *molecule) const {
  int m;

  if (nex_type && ex_type[itype][jtype]) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask[i] & ex1_bit[m] && mask[j] & ex2_bit[m]) return 1;
      if (mask[i] & ex2_bit[m] && mask[j] & ex1_bit[m]) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)
      if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
          molecule[i] == molecule[j]) return 1;
  }

  return 0;
}
