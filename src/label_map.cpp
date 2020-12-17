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
}

/* ----------------------------------------------------------------------
   find integer type given a type label
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::find_type(char *mytypelabel, char **typelabelarray, int num_types)
{
  for (int i = 0; i < num_types; i++) {
    if (typelabelarray[i] && strcmp(mytypelabel,typelabelarray[i]) == 0) return i+1;
  }
  return -1;
}
