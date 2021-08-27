// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_update_special_bonds.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"

#include <cstring>
#include <set>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixUpdateSpecialBonds::FixUpdateSpecialBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix update/special/bonds command");
  comm_forward = 1+atom->maxspecial;  
}

/* ---------------------------------------------------------------------- */

FixUpdateSpecialBonds::~FixUpdateSpecialBonds()
{
}

/* ---------------------------------------------------------------------- */

int FixUpdateSpecialBonds::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::setup(int /*vflag*/)
{
  // Require atoms know about all of their bonds and if they break
  if (force->newton_bond)
    error->all(FLERR,"Fix update/special/bonds requires Newton bond off");

  if (!atom->avec->bonds_allow)
    error->all(FLERR,"Fix update/special/bonds requires atom bonds");

  // special lj must be 0 1 1 to censor pair forces between bonded particles
  // special coulomb must be 1 1 1 to ensure all pairs are included in the 
  //   neighbor list and 1-3 and 1-4 special bond lists are skipped  
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0)
    error->all(FLERR,"Fix update/special/bonds requires special LJ weights = 0,1,1");
  if (force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0)
    error->all(FLERR,"Fix update/special/bonds requires special Coulomb weights = 1,1,1");
    
  broken_pairs.clear();
}

/* ----------------------------------------------------------------------
  Update special bond list and atom bond arrays, empty broken bond list
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_exchange()
{
  int i, j, key, m, n1, n3;
  tagint min_tag, max_tag;
  int nlocal = atom->nlocal;
  
  tagint *tag = atom->tag;  
  tagint *slist;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  
  for (auto const &key : broken_pairs) {
    min_tag = key.first;
    max_tag = key.second;
    
    i = atom->map(min_tag);
    j = atom->map(max_tag);
    
    // remove i from special bond list for atom j and vice versa
    if (i < nlocal) {
      slist = special[i];
      n1 = nspecial[i][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == max_tag) break;
      n3 = nspecial[i][2];
      for (; m < n3-1; m++) slist[m] = slist[m+1];
      nspecial[i][0]--;
      nspecial[i][1]--;
      nspecial[i][2]--;
    }
    
    if (j < nlocal) {
      slist = special[j];
      n1 = nspecial[j][0];
      for (m = 0; m < n1; m++)
        if (slist[m] == min_tag) break;
      n3 = nspecial[j][2];
      for (; m < n3-1; m++) slist[m] = slist[m+1];
      nspecial[j][0]--;
      nspecial[j][1]--;
      nspecial[j][2]--;
    }
  }
  
  // Forward updated special bond list
  comm->forward_comm_fix(this);
  
  broken_pairs.clear();
}

/* ----------------------------------------------------------------------
  Loop neighbor list and update special bond lists for recently broken bonds
------------------------------------------------------------------------- */

void FixUpdateSpecialBonds::pre_force(int /*vflag*/)
{
  int i,j,n,m,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  tagint min_tag, max_tag;
  std::pair <tagint, tagint> key;
  
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  
  tagint *tag = atom->tag;
  NeighList *list = force->pair->list;  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      min_tag = tag[i];
      max_tag = tag[j];
      if (max_tag < min_tag) {
        min_tag = tag[j];
        max_tag = tag[i];
      }
      key = std::make_pair(min_tag, max_tag);

      if (broken_pairs.find(key) != broken_pairs.end()) 
        jlist[jj] = j; // Clear special bond bits
    } 
  }
}

/* ---------------------------------------------------------------------- */

int FixUpdateSpecialBonds::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m,ns;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ns = nspecial[j][0];
    buf[m++] = ubuf(ns).d;
    for (k = 0; k < ns; k++)
      buf[m++] = ubuf(special[j][k]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,m,ns,last;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    ns = (int) ubuf(buf[m++]).i;
    nspecial[i][0] = ns;
    for (j = 0; j < ns; j++)
      special[i][j] = (tagint) ubuf(buf[m++]).i;
  }
}

/* ---------------------------------------------------------------------- */

void FixUpdateSpecialBonds::add_broken_bond(int i, int j)
{
  tagint *tag = atom->tag;
  
  tagint min_tag = tag[i];
  tagint max_tag = tag[j];
  if (max_tag < min_tag) {
    min_tag = tag[j];
    max_tag = tag[i];
  }
  std::pair <tagint, tagint> key = std::make_pair(min_tag, max_tag);

  broken_pairs.insert(key);
}
