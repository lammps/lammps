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
  lmap2lmap.atom = new int[natomtypes];

  if (force->bond) {
    btypelabel.resize(nbondtypes);
    lmap2lmap.bond = new int[nbondtypes];
  }

  if (force->angle) {
    atypelabel.resize(nangletypes);
    lmap2lmap.angle = new int[nangletypes];
  }

  if (force->dihedral) {
    dtypelabel.resize(ndihedraltypes);
    lmap2lmap.dihedral = new int[ndihedraltypes];
  }

  if (force->improper) {
    itypelabel.resize(nimpropertypes);
    lmap2lmap.improper = new int[nimpropertypes];
  }
}

/* ----------------------------------------------------------------------
   copy another map (lmap2) into this one
   if label already exists, leave in place
   else, put new label in next available slot
------------------------------------------------------------------------- */

void LabelMap::merge_lmap(LabelMap *lmap2, int mode)
{
  if (mode == ATOM)
    for (int i = 0; i < lmap2->natomtypes; i++)
      find_or_create(lmap2->typelabel[i],typelabel,natomtypes);

  if (mode == BOND)
    for (int i = 0; i < lmap2->nbondtypes; i++)
      find_or_create(lmap2->btypelabel[i],btypelabel,nbondtypes);

  if (mode == ANGLE)
    for (int i = 0; i < lmap2->nangletypes; i++)
      find_or_create(lmap2->atypelabel[i],atypelabel,nangletypes);

  if (mode == DIHEDRAL)
    for (int i = 0; i < lmap2->ndihedraltypes; i++)
      find_or_create(lmap2->dtypelabel[i],dtypelabel,ndihedraltypes);

  if (mode == IMPROPER)
    for (int i = 0; i < lmap2->nimpropertypes; i++)
      find_or_create(lmap2->itypelabel[i],itypelabel,nimpropertypes);
}

/* ----------------------------------------------------------------------
   get mapping between this label map and another (lmap2)
   values of lmap2lmap point to equivalent types in lmap2
------------------------------------------------------------------------- */

void LabelMap::create_lmap2lmap(LabelMap *lmap2, int mode)
{
  if (mode == ATOM)
    for (int i = 0; i < natomtypes; i++)
      lmap2lmap.atom[i] = find(typelabel[i],lmap2->typelabel,
                               lmap2->natomtypes);

  if (mode == BOND)
    for (int i = 0; i < nbondtypes; i++)
      lmap2lmap.bond[i] = find(btypelabel[i],lmap2->btypelabel,
                               lmap2->nbondtypes);

  if (mode == ANGLE)
    for (int i = 0; i < nangletypes; i++)
      lmap2lmap.angle[i] = find(atypelabel[i],lmap2->atypelabel,
                                lmap2->nangletypes);

  if (mode == DIHEDRAL)
    for (int i = 0; i < ndihedraltypes; i++)
      lmap2lmap.dihedral[i] = find(dtypelabel[i],lmap2->dtypelabel,
                                   lmap2->ndihedraltypes);

  if (mode == IMPROPER)
    for (int i = 0; i < nimpropertypes; i++)
      lmap2lmap.improper[i] = find(itypelabel[i],lmap2->itypelabel,
                                   lmap2->nimpropertypes);
}

/* ----------------------------------------------------------------------
   find type label with name or create type if it doesn't exist
   return numeric type
------------------------------------------------------------------------- */

int LabelMap::find_or_create(std::string mylabel, std::vector<std::string> &labels, int ntypes)
{
  for (int i = 0; i < ntypes; i++)
    if (labels[i] == mylabel) return i+1;

  // if no match found, create new label at next available index
  // label map assumed to be intialized with numeric index
  // user labels are assumed to be alphanumeric (not a number)

  for (int i = 0; i < ntypes; i++) {
    if (labels[i].empty()) {
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

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void LabelMap::write_data(FILE *fp)
{
  fmt::print(fp,"\nAtom Type Labels\n\n");
  for (int i = 0; i < natomtypes; i++)
    fmt::print(fp,"{} {}\n",i+1,typelabel[i]);

  if (force->bond) {
    fmt::print(fp,"\nBond Type Labels\n\n");
    for (int i = 0; i < nbondtypes; i++)
      fmt::print(fp,"{} {}\n",i+1,btypelabel[i]);
  }

  if (force->angle) {
    fmt::print(fp,"\nAngle Type Labels\n\n");
    for (int i = 0; i < nangletypes; i++)
      fmt::print(fp,"{} {}\n",i+1,atypelabel[i]);
  }

  if (force->dihedral) {
    fmt::print(fp,"\nDihedral Type Labels\n\n");
    for (int i = 0; i < ndihedraltypes; i++)
      fmt::print(fp,"{} {}\n",i+1,dtypelabel[i]);
  }

  if (force->improper) {
    fmt::print(fp,"\nImproper Type Labels\n\n");
    for (int i = 0; i < nimpropertypes; i++)
      fmt::print(fp,"{} {}\n",i+1,itypelabel[i]);
  }
}
