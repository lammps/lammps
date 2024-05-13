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

#include "npair_copy.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairCopy::NPairCopy(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
------------------------------------------------------------------------- */

void NPairCopy::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->ipage = listcopy->ipage;
}
