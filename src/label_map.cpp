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

#include "label_map.h"

#include "force.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

LabelMap::LabelMap(LAMMPS *lmp) : Pointers(lmp)
{
  natomtypes = nbondtypes = nangletypes = 0;
  ndihedraltypes = nimpropertypes = 0;
}

/* ---------------------------------------------------------------------- */

LabelMap::~LabelMap()
{
  // delete type labels

  typelabel.clear();
  btypelabel.clear();
  atypelabel.clear();
  dtypelabel.clear();
  itypelabel.clear();

  delete [] lmap2lmap.atom;
  delete [] lmap2lmap.bond;
  delete [] lmap2lmap.angle;
  delete [] lmap2lmap.dihedral;
  delete [] lmap2lmap.improper;
}

/* ----------------------------------------------------------------------
   allocate character-based type arrays (labels) of length ntypes
   always allocated (for both numeric and character-based type modes)
   initialize label with (a string of) its numeric counterpart
------------------------------------------------------------------------- */

void LabelMap::allocate_type_labels()
{
  typelabel.resize(natomtypes);
  for (int i = 0; i < natomtypes; i++)
    typelabel[i] = fmt::format("{}",i);

  if (force->bond) {
    btypelabel.resize(nbondtypes);
    for (int i = 0; i < nbondtypes; i++)
      btypelabel[i] = fmt::format("{}",i);
  }

  if (force->angle) {
    atypelabel.resize(nangletypes);
    for (int i = 0; i < nangletypes; i++)
      atypelabel[i] = fmt::format("{}",i);
  }

  if (force->dihedral) {
    dtypelabel.resize(ndihedraltypes);
    for (int i = 0; i < ndihedraltypes; i++)
      dtypelabel[i] = fmt::format("{}",i);
  }

  if (force->improper) {
    itypelabel.resize(nimpropertypes);
    for (int i = 0; i < nimpropertypes; i++)
      itypelabel[i] = fmt::format("{}",i);
  }

  // allocate space for lmap2lmap

  lmap2lmap.atom = new int[natomtypes];
  lmap2lmap.bond = new int[nbondtypes];
  lmap2lmap.angle = new int[nangletypes];
  lmap2lmap.dihedral = new int[ndihedraltypes];
  lmap2lmap.improper = new int[nimpropertypes];
}

/* ----------------------------------------------------------------------
   copy another map (lmap2) into this one
   if label already exists, leave in place
   else, put new label in next available slot
------------------------------------------------------------------------- */

void LabelMap::merge_lmap(class LabelMap *lmap2)
{
  for (int i = 0; i < lmap2->natomtypes; i++)
    find_or_create(lmap2->typelabel[i],typelabel,natomtypes);

  if (force->bond) {
    for (int i = 0; i < lmap2->nbondtypes; i++)
      find_or_create(lmap2->btypelabel[i],btypelabel,nbondtypes);
  }

  if (force->angle) {
    for (int i = 0; i < lmap2->nangletypes; i++)
      find_or_create(lmap2->atypelabel[i],atypelabel,nangletypes);
  }

  if (force->dihedral) {
    for (int i = 0; i < lmap2->ndihedraltypes; i++)
      find_or_create(lmap2->dtypelabel[i],dtypelabel,ndihedraltypes);
  }

  if (force->improper) {
    for (int i = 0; i < lmap2->nimpropertypes; i++)
      find_or_create(lmap2->itypelabel[i],itypelabel,nimpropertypes);
  }
}

/* ----------------------------------------------------------------------
   get mapping between this label map and another (lmap2)
   values of lmap2lmap point to equivalent indices in lmap2
------------------------------------------------------------------------- */

void LabelMap::create_lmap2lmap(class LabelMap *lmap2)
{
  int type;

  for (int i = 0; i < natomtypes; i++) {
    type = find(typelabel[i],lmap2->typelabel,lmap2->natomtypes);
    lmap2lmap.atom[i] = type - 1;
  }

  if (force->bond) {
    lmap2lmap.bond = new int[nbondtypes];
    for (int i = 0; i < nbondtypes; i++) {
      type = find(btypelabel[i],lmap2->btypelabel,lmap2->nbondtypes);
      lmap2lmap.bond[i] = type - 1;
    }
  }

  if (force->angle) {
    for (int i = 0; i < nangletypes; i++) {
      type = find(atypelabel[i],lmap2->atypelabel,lmap2->nangletypes);
      lmap2lmap.angle[i] = type - 1;
    }
  }

  if (force->dihedral) {
    for (int i = 0; i < ndihedraltypes; i++) {
      type = find(dtypelabel[i],lmap2->dtypelabel,lmap2->ndihedraltypes);
      lmap2lmap.dihedral[i] = type - 1;
    }
  }

  if (force->improper) {
    for (int i = 0; i < nimpropertypes; i++) {
      type = find(itypelabel[i],lmap2->itypelabel,lmap2->nimpropertypes);
      lmap2lmap.improper[i] = type - 1;
    }
  }
}

/* ----------------------------------------------------------------------
   find type label with name or create type if it doesn't exist
   return numeric type
------------------------------------------------------------------------- */

int LabelMap::find_or_create(std::string mylabel, std::vector<std::string> labels, int ntypes)
{
  for (int i = 0; i < ntypes; i++)
    if (labels[i] == mylabel) return i+1;

  // if no match found, create new label at next available index
  // label map assumed to be intialized with numeric index
  // user labels are assumed to be alphanumeric (not a number)

  for (int i = 0; i < ntypes; i++) {
    if (stoi(labels[i]) != i+1) {
      labels[i] = mylabel;
      return i+1;
    }
  }

  // if label cannot be found or created, need more space reserved

  error->all(FLERR,"Topology type exceeds system topology type");
}

/* ----------------------------------------------------------------------
   find integer type given a type label
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::find(std::string mylabel, std::vector<std::string> labels, int ntypes)
{
  for (int i = 0; i < ntypes; i++)
    if (labels[i] == mylabel) return i+1;

  return -1;
}
