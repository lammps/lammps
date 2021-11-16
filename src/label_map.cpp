/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "label_map.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
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
  delete[] lmap2lmap.atom;
  delete[] lmap2lmap.bond;
  delete[] lmap2lmap.angle;
  delete[] lmap2lmap.dihedral;
  delete[] lmap2lmap.improper;
}

/* ----------------------------------------------------------------------
   allocate character-based type arrays (labels) of length ntypes
------------------------------------------------------------------------- */

void LabelMap::allocate_type_labels()
{
  typelabel.resize(natomtypes);
  lmap2lmap.atom = new int[natomtypes];

  btypelabel.resize(nbondtypes);
  lmap2lmap.bond = new int[nbondtypes];

  atypelabel.resize(nangletypes);
  lmap2lmap.angle = new int[nangletypes];

  dtypelabel.resize(ndihedraltypes);
  lmap2lmap.dihedral = new int[ndihedraltypes];

  itypelabel.resize(nimpropertypes);
  lmap2lmap.improper = new int[nimpropertypes];
}

/* ----------------------------------------------------------------------
   labelmap command in input script
------------------------------------------------------------------------- */

void LabelMap::modify_lmap(int narg, char **arg)
{
  if (narg < 3 || narg % 2 == 0) error->all(FLERR,"Illegal labelmap command");

  int ntypes;
  std::vector<std::string> *labels;
  const std::string tlabel(arg[0]);
  if (tlabel == "atom") {
    ntypes = natomtypes;
    labels = &typelabel;
  } else if (tlabel == "bond") {
    ntypes = nbondtypes;
    labels = &btypelabel;
  } else if (tlabel == "angle") {
    ntypes = nangletypes;
    labels = &atypelabel;
  } else if (tlabel == "dihedral") {
    ntypes = ndihedraltypes;
    labels = &dtypelabel;
  } else if (tlabel == "improper") {
    ntypes = nimpropertypes;
    labels = &itypelabel;
  } else error->all(FLERR,"Illegal labelmap command");

  int itype;
  int iarg = 1;
  while (iarg < narg) {
    itype = utils::inumeric(FLERR,arg[iarg++],false,lmp);
    if (itype > ntypes)
      error->all(FLERR,"Topology type exceeds system topology type");
    std::string slabel(arg[iarg++]);
    if (isdigit(slabel[0]))
      error->all(FLERR,"Type labels cannot start with a number");
    if (search(slabel,(*labels),ntypes) != -1)
      error->all(FLERR,"Type label already exists: types labels must be unique");
    (*labels)[itype-1] = slabel;
  }
}

/* ----------------------------------------------------------------------
   copy another map (lmap2) into this one
   if label already exists, leave in place
   else, put new label in next available slot
------------------------------------------------------------------------- */

void LabelMap::merge_lmap(LabelMap *lmap2, int mode)
{
  switch (mode)
  {
  case Atom::ATOM:
    for (int i = 0; i < lmap2->natomtypes; i++)
      find_or_create(lmap2->typelabel[i],typelabel,natomtypes);
    break;
  case Atom::BOND:
    for (int i = 0; i < lmap2->nbondtypes; i++)
      find_or_create(lmap2->btypelabel[i],btypelabel,nbondtypes);
    break;
  case Atom::ANGLE:
    for (int i = 0; i < lmap2->nangletypes; i++)
      find_or_create(lmap2->atypelabel[i],atypelabel,nangletypes);
    break;
  case Atom::DIHEDRAL:
    for (int i = 0; i < lmap2->ndihedraltypes; i++)
      find_or_create(lmap2->dtypelabel[i],dtypelabel,ndihedraltypes);
    break;
  case Atom::IMPROPER:
    for (int i = 0; i < lmap2->nimpropertypes; i++)
      find_or_create(lmap2->itypelabel[i],itypelabel,nimpropertypes);
    break;
  }
}

/* ----------------------------------------------------------------------
   get mapping between this label map and another (lmap2)
   values of lmap2lmap point to equivalent types in lmap2
------------------------------------------------------------------------- */

void LabelMap::create_lmap2lmap(LabelMap *lmap2, int mode)
{
  if (mode == Atom::ATOM)
    for (int i = 0; i < natomtypes; i++)
      lmap2lmap.atom[i] = search(typelabel[i],lmap2->typelabel,
                                 lmap2->natomtypes);

  if (mode == Atom::BOND)
    for (int i = 0; i < nbondtypes; i++)
      lmap2lmap.bond[i] = search(btypelabel[i],lmap2->btypelabel,
                                 lmap2->nbondtypes);

  if (mode == Atom::ANGLE)
    for (int i = 0; i < nangletypes; i++)
      lmap2lmap.angle[i] = search(atypelabel[i],lmap2->atypelabel,
                                  lmap2->nangletypes);

  if (mode == Atom::DIHEDRAL)
    for (int i = 0; i < ndihedraltypes; i++)
      lmap2lmap.dihedral[i] = search(dtypelabel[i],lmap2->dtypelabel,
                                     lmap2->ndihedraltypes);

  if (mode == Atom::IMPROPER)
    for (int i = 0; i < nimpropertypes; i++)
      lmap2lmap.improper[i] = search(itypelabel[i],lmap2->itypelabel,
                                     lmap2->nimpropertypes);
}

/* ----------------------------------------------------------------------
   find type label with name or create type if it doesn't exist
   return numeric type
------------------------------------------------------------------------- */

int LabelMap::find_or_create(const std::string &mylabel,
                             std::vector<std::string> &labels, int ntypes)
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

  // never reaches here, just to prevent compiler warning

  return -1;
}

/* ----------------------------------------------------------------------
   return numeric type given a type label
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::find(const std::string &mylabel, int mode)
{
  if (mode == Atom::ATOM)
    return search(mylabel,typelabel,natomtypes);

  if (mode == Atom::BOND)
    return search(mylabel,btypelabel,nbondtypes);

  if (mode == Atom::ANGLE)
    return search(mylabel,atypelabel,nangletypes);

  if (mode == Atom::DIHEDRAL)
    return search(mylabel,dtypelabel,ndihedraltypes);

  if (mode == Atom::IMPROPER)
    return search(mylabel,itypelabel,nimpropertypes);

  return -1;
}

/* ----------------------------------------------------------------------
   get index+1 given vector of strings
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::search(const std::string &mylabel,
                     const std::vector<std::string> &labels, int ntypes)
{
  for (int i = 0; i < ntypes; i++)
    if (labels[i] == mylabel) return i+1;

  return -1;
}

/* ----------------------------------------------------------------------
   check that all types have been assigned a unique type label
------------------------------------------------------------------------- */

int LabelMap::is_complete(int mode)
{
  int i,j;

  if (mode == Atom::ATOM)
    for (i = 0; i < natomtypes; i++) {
      if (typelabel[i].empty()) return 0;
      for (j = i+1; j < natomtypes; j++)
        if (typelabel[i] == typelabel[j]) return 0;
    }

  if (force->bond && mode == Atom::BOND)
    for (i = 0; i < nbondtypes; i++) {
      if (btypelabel[i].empty()) return 0;
      for (j = i+1; j < nbondtypes; j++)
        if (btypelabel[i] == btypelabel[j]) return 0;
    }

  if (force->angle && mode == Atom::ANGLE)
    for (i = 0; i < nangletypes; i++) {
      if (atypelabel[i].empty()) return 0;
      for (j = i+1; j < nangletypes; j++)
        if (atypelabel[i] == atypelabel[j]) return 0;
    }

  if (force->dihedral && mode == Atom::DIHEDRAL)
    for (i = 0; i < ndihedraltypes; i++) {
      if (dtypelabel[i].empty()) return 0;
      for (j = i+1; j < ndihedraltypes; j++)
        if (dtypelabel[i] == dtypelabel[j]) return 0;
    }

  if (force->improper && mode == Atom::IMPROPER)
    for (i = 0; i < nimpropertypes; i++) {
      if (itypelabel[i].empty()) return 0;
      for (j = i+1; j < nimpropertypes; j++)
        if (itypelabel[i] == itypelabel[j]) return 0;
    }

  return 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void LabelMap::write_data(FILE *fp)
{
  if (is_complete(Atom::ATOM)) {
    fmt::print(fp,"\nAtom Type Labels\n\n");
    for (int i = 0; i < natomtypes; i++)
      fmt::print(fp,"{} {}\n",i+1,typelabel[i]);
  }

  if (force->bond && is_complete(Atom::BOND)) {
    fmt::print(fp,"\nBond Type Labels\n\n");
    for (int i = 0; i < nbondtypes; i++)
      fmt::print(fp,"{} {}\n",i+1,btypelabel[i]);
  }

  if (force->angle && is_complete(Atom::ANGLE)) {
    fmt::print(fp,"\nAngle Type Labels\n\n");
    for (int i = 0; i < nangletypes; i++)
      fmt::print(fp,"{} {}\n",i+1,atypelabel[i]);
  }

  if (force->dihedral && is_complete(Atom::DIHEDRAL)) {
    fmt::print(fp,"\nDihedral Type Labels\n\n");
    for (int i = 0; i < ndihedraltypes; i++)
      fmt::print(fp,"{} {}\n",i+1,dtypelabel[i]);
  }

  if (force->improper && is_complete(Atom::IMPROPER)) {
    fmt::print(fp,"\nImproper Type Labels\n\n");
    for (int i = 0; i < nimpropertypes; i++)
      fmt::print(fp,"{} {}\n",i+1,itypelabel[i]);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void LabelMap::read_restart(FILE *fp)
{
  char *charlabel;

  for (int i = 0; i < natomtypes; i++) {
    charlabel = read_string(fp);
    typelabel[i] = charlabel;
    delete[] charlabel;
  }

  for (int i = 0; i < nbondtypes; i++) {
    charlabel = read_string(fp);
    btypelabel[i] = charlabel;
    delete[] charlabel;
  }

  for (int i = 0; i < nangletypes; i++) {
    charlabel = read_string(fp);
    atypelabel[i] = charlabel;
    delete[] charlabel;
  }

  for (int i = 0; i < ndihedraltypes; i++) {
    charlabel = read_string(fp);
    dtypelabel[i] = charlabel;
    delete[] charlabel;
  }

  for (int i = 0; i < nimpropertypes; i++) {
    charlabel = read_string(fp);
    itypelabel[i] = charlabel;
    delete[] charlabel;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void LabelMap::write_restart(FILE *fp)
{
  for (int i = 0; i < natomtypes; i++)
    write_string(typelabel[i],fp);

  for (int i = 0; i < nbondtypes; i++)
    write_string(btypelabel[i],fp);

  for (int i = 0; i < nangletypes; i++)
    write_string(atypelabel[i],fp);

  for (int i = 0; i < ndihedraltypes; i++)
    write_string(dtypelabel[i],fp);

  for (int i = 0; i < nimpropertypes; i++)
    write_string(itypelabel[i],fp);
}

/* ----------------------------------------------------------------------
   read a char string (including nullptr) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *LabelMap::read_string(FILE *fp)
{
  int n = read_int(fp);
  if (n < 0) error->all(FLERR,"Illegal size string or corrupt restart");
  char *value = new char[n];
  if (comm->me == 0)
    utils::sfread(FLERR,value,sizeof(char),n,fp,nullptr,error);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   write a flag and a C-style char string (including the terminating null
   byte) into the restart file
------------------------------------------------------------------------- */

void LabelMap::write_string(const std::string &str, FILE *fp)
{
  const char *cstr = str.c_str();
  int n = strlen(cstr) + 1;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(cstr,sizeof(char),n,fp);
}

/* ----------------------------------------------------------------------
   read an int from restart file and bcast it
------------------------------------------------------------------------- */

int LabelMap::read_int(FILE *fp)
{
  int value;
  if ((comm->me == 0) && (fread(&value,sizeof(int),1,fp) < 1)) value = -1;
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}
