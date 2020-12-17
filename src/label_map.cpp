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

  typelabel = btypelabel = atypelabel = NULL;
  dtypelabel = itypelabel = NULL;
}

/* ---------------------------------------------------------------------- */

LabelMap::~LabelMap()
{
  // delete type labels

  for (int i = 0; i < natomtypes; i++) delete typelabel[i];
  memory->sfree(typelabel);
  for (int i = 0; i < nbondtypes; i++) delete btypelabel[i];
  memory->sfree(btypelabel);
  for (int i = 0; i < nangletypes; i++) delete atypelabel[i];
  memory->sfree(atypelabel);
  for (int i = 0; i < ndihedraltypes; i++) delete dtypelabel[i];
  memory->sfree(dtypelabel);
  for (int i = 0; i < nimpropertypes; i++) delete itypelabel[i];
  memory->sfree(itypelabel);
}

/* ----------------------------------------------------------------------
   allocate character-based type arrays (labels) of length ntypes
   always allocated (for both numeric and character-based type modes)
   initialize label with (a string of) its numeric counterpart
------------------------------------------------------------------------- */

void LabelMap::allocate_type_labels()
{
  char *char_type = new char[256];

  typelabel = (char **) memory->srealloc(typelabel,
             natomtypes*sizeof(char *),"atom:typelabel");
  for (int i = 0; i < natomtypes; i++) {
    sprintf(char_type,"%d",i+1);
    int n = strlen(char_type) + 1;
    typelabel[i] = new char[n];
    strcpy(typelabel[i],char_type);
  }
  if (force->bond) {
    btypelabel = (char **) memory->srealloc(btypelabel,
                nbondtypes*sizeof(char *),"atom:btypelabel");
    for (int i = 0; i < nbondtypes; i++) {
      sprintf(char_type,"%d",i+1);
      int n = strlen(char_type) + 1;
      btypelabel[i] = new char[n];
      strcpy(btypelabel[i],char_type);
    }
  }
  if (force->angle) {
    atypelabel = (char **) memory->srealloc(atypelabel,
                nangletypes*sizeof(char *),"atom:atypelabel");
    for (int i = 0; i < nangletypes; i++) {
      sprintf(char_type,"%d",i+1);
      int n = strlen(char_type) + 1;
      atypelabel[i] = new char[n];
      strcpy(atypelabel[i],char_type);
    }
  }
  if (force->dihedral) {
    dtypelabel = (char **) memory->srealloc(dtypelabel,
                ndihedraltypes*sizeof(char *),"atom:dtypelabel");
    for (int i = 0; i < ndihedraltypes; i++) {
      sprintf(char_type,"%d",i+1);
      int n = strlen(char_type) + 1;
      dtypelabel[i] = new char[n];
      strcpy(dtypelabel[i],char_type);
    }
  }
  if (force->improper) {
    itypelabel = (char **) memory->srealloc(itypelabel,
                nimpropertypes*sizeof(char *),"atom:itypelabel");
    for (int i = 0; i < nimpropertypes; i++) {
      sprintf(char_type,"%d",i+1);
      int n = strlen(char_type) + 1;
      itypelabel[i] = new char[n];
      strcpy(itypelabel[i],char_type);
    }
  }
  delete [] char_type;
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
