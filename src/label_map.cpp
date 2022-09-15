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

LabelMap::LabelMap(LAMMPS *_lmp, int _natomtypes, int _nbondtypes, int _nangletypes,
                   int _ndihedraltypes, int _nimpropertypes) :
    Pointers(_lmp),
    natomtypes(_natomtypes), nbondtypes(_nbondtypes), nangletypes(_nangletypes),
    ndihedraltypes(_ndihedraltypes), nimpropertypes(_nimpropertypes)
{
  lmap2lmap.atom = lmap2lmap.bond = lmap2lmap.angle = lmap2lmap.dihedral = lmap2lmap.improper =
      nullptr;
  reset_type_labels();
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
   reset/allocate character-based type arrays (labels) of length ntypes
------------------------------------------------------------------------- */

void LabelMap::reset_type_labels()
{
  typelabel_map.clear();
  typelabel.resize(natomtypes);
  delete[] lmap2lmap.atom;
  lmap2lmap.atom = new int[natomtypes];
  for (auto &i : typelabel) i.clear();
  memset(lmap2lmap.atom, 0, natomtypes * sizeof(int));

  btypelabel_map.clear();
  btypelabel.resize(nbondtypes);
  delete[] lmap2lmap.bond;
  for (auto &i : btypelabel) i.clear();
  lmap2lmap.bond = new int[nbondtypes];
  memset(lmap2lmap.bond, 0, nbondtypes * sizeof(int));

  atypelabel_map.clear();
  atypelabel.resize(nangletypes);
  delete[] lmap2lmap.angle;
  for (auto &i : atypelabel) i.clear();
  lmap2lmap.angle = new int[nangletypes];
  memset(lmap2lmap.angle, 0, nangletypes * sizeof(int));

  dtypelabel_map.clear();
  dtypelabel.resize(ndihedraltypes);
  delete[] lmap2lmap.dihedral;
  for (auto &i : dtypelabel) i.clear();
  lmap2lmap.dihedral = new int[ndihedraltypes];
  memset(lmap2lmap.dihedral, 0, ndihedraltypes * sizeof(int));

  itypelabel_map.clear();
  itypelabel.resize(nimpropertypes);
  delete[] lmap2lmap.improper;
  for (auto &i : itypelabel) i.clear();
  lmap2lmap.improper = new int[nimpropertypes];
  memset(lmap2lmap.improper, 0, nimpropertypes * sizeof(int));
}

/* ----------------------------------------------------------------------
   labelmap command in input script
------------------------------------------------------------------------- */

void LabelMap::modify_lmap(int narg, char **arg)
{
  if ((narg < 1) || ((narg > 2) && ((narg % 2) == 0)))
    error->all(FLERR, "Incorrect number of arguments for labelmap command");

  int ntypes;
  std::vector<std::string> *labels;
  std::unordered_map<std::string, int> *labels_map;
  const std::string tlabel(arg[0]);
  if (tlabel == "atom") {
    ntypes = natomtypes;
    labels = &typelabel;
    labels_map = &typelabel_map;
  } else if (tlabel == "bond") {
    ntypes = nbondtypes;
    labels = &btypelabel;
    labels_map = &btypelabel_map;
  } else if (tlabel == "angle") {
    ntypes = nangletypes;
    labels = &atypelabel;
    labels_map = &atypelabel_map;
  } else if (tlabel == "dihedral") {
    ntypes = ndihedraltypes;
    labels = &dtypelabel;
    labels_map = &dtypelabel_map;
  } else if (tlabel == "improper") {
    ntypes = nimpropertypes;
    labels = &itypelabel;
    labels_map = &itypelabel_map;
  } else if (tlabel == "clear") {
    if (narg != 1) error->all(FLERR, "Incorrect number of arguments for labelmap clear command");
    reset_type_labels();
    return;
  } else if (tlabel == "write") {
    if (narg != 2) error->all(FLERR, "Incorrect number of arguments for labelmap write command");
    write_map(arg[1]);
    return;
  } else
    error->all(FLERR, "Unknown labelmap keyword {}", tlabel);

  int iarg = 1;
  if (narg == 1) utils::missing_cmd_args(FLERR, "labelmap " + tlabel, error);
  while (iarg < narg) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "labelmap " + tlabel, error);
    if (ntypes < 1) error->all(FLERR, "No {} types allowed with current box settings", tlabel);
    int itype = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    if ((itype < 1) || (itype > ntypes))
      error->all(FLERR, "Labelmap {} type {} must be within 1-{}", tlabel, itype, ntypes);
    std::string slabel = utils::utf8_subst(utils::trim(arg[iarg++]));
    if (utils::is_type(slabel) != 1)
      error->all(FLERR, "Type label string {} for {} type {} is invalid", slabel, tlabel, itype);
    int found = search(slabel, (*labels_map));
    if ((found != -1) && (found != itype))
      error->all(FLERR, "The {} type label {} is already in use for type {}", tlabel, slabel,
                 (*labels_map)[slabel]);
    std::string &str = (*labels)[itype - 1];
    if (!str.empty()) (*labels_map).erase(str);
    str = slabel;
    (*labels_map)[slabel] = itype;
  }
}

/* ----------------------------------------------------------------------
   copy another map (lmap2) into this one
   if label already exists, leave in place
   else, put new label in next available slot
------------------------------------------------------------------------- */

void LabelMap::merge_lmap(LabelMap *lmap2, int mode)
{
  switch (mode) {
    case Atom::ATOM:
      for (auto &it : lmap2->typelabel) find_or_create(it, typelabel, typelabel_map);
      break;
    case Atom::BOND:
      for (auto &it : lmap2->btypelabel) find_or_create(it, btypelabel, btypelabel_map);
      break;
    case Atom::ANGLE:
      for (auto &it : lmap2->atypelabel) find_or_create(it, atypelabel, atypelabel_map);
      break;
    case Atom::DIHEDRAL:
      for (auto &it : lmap2->dtypelabel) find_or_create(it, dtypelabel, dtypelabel_map);
      break;
    case Atom::IMPROPER:
      for (auto &it : lmap2->itypelabel) find_or_create(it, itypelabel, itypelabel_map);
      break;
  }
}

/* ----------------------------------------------------------------------
   get mapping between this label map and another (lmap2)
   values of lmap2lmap point to equivalent types in lmap2
------------------------------------------------------------------------- */

void LabelMap::create_lmap2lmap(LabelMap *lmap2, int mode)
{
  switch (mode) {
    case Atom::ATOM:
      for (int i = 0; i < natomtypes; ++i)
        lmap2lmap.atom[i] = search(typelabel[i], lmap2->typelabel_map);
      break;
    case Atom::BOND:
      for (int i = 0; i < nbondtypes; ++i)
        lmap2lmap.bond[i] = search(btypelabel[i], lmap2->btypelabel_map);
      break;
    case Atom::ANGLE:
      for (int i = 0; i < nangletypes; ++i)
        lmap2lmap.angle[i] = search(atypelabel[i], lmap2->atypelabel_map);
      break;
    case Atom::DIHEDRAL:
      for (int i = 0; i < ndihedraltypes; ++i)
        lmap2lmap.dihedral[i] = search(dtypelabel[i], lmap2->dtypelabel_map);
      break;
    case Atom::IMPROPER:
      for (int i = 0; i < nimpropertypes; ++i)
        lmap2lmap.improper[i] = search(itypelabel[i], lmap2->itypelabel_map);
      break;
  }
}

/* ----------------------------------------------------------------------
   find type label with name or create type if it doesn't exist
   return numeric type
------------------------------------------------------------------------- */

int LabelMap::find_or_create(const std::string &mylabel, std::vector<std::string> &labels,
                             std::unordered_map<std::string, int> &labels_map)
{
  auto search = labels_map.find(mylabel);
  if (search != labels_map.end()) return search->second;

  // if no match found, create new label at next available index
  // label map assumed to be intialized with numeric index
  // user labels are assumed to be alphanumeric (not a number)

  auto labels_map_size = labels_map.size();
  if (labels_map_size < labels.size()) {
    labels[labels_map_size] = mylabel;
    int index = static_cast<int>(labels_map_size + 1);
    labels_map[mylabel] = index;
    return index;
  }

  // if label cannot be found or created, need more space reserved

  error->all(FLERR, "Topology type exceeds system topology type");

  // never reaches here, just to prevent compiler warning

  return -1;
}

/* ----------------------------------------------------------------------
   return numeric type given a type label
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::find(const std::string &mylabel, int mode) const
{
  switch (mode) {
    case Atom::ATOM:
      return search(mylabel, typelabel_map);
      break;
    case Atom::BOND:
      return search(mylabel, btypelabel_map);
      break;
    case Atom::ANGLE:
      return search(mylabel, atypelabel_map);
      break;
    case Atom::DIHEDRAL:
      return search(mylabel, dtypelabel_map);
      break;
    case Atom::IMPROPER:
      return search(mylabel, itypelabel_map);
      break;
    default:
      return -1;
  }
}

/* ----------------------------------------------------------------------
   get type given type labels map
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::search(const std::string &mylabel,
                     const std::unordered_map<std::string, int> &labels_map) const
{
  auto search = labels_map.find(mylabel);
  if (search == labels_map.end()) return -1;
  return search->second;
}

/* ----------------------------------------------------------------------
   check that all types have been assigned a unique type label
------------------------------------------------------------------------- */

bool LabelMap::is_complete(int mode) const
{
  switch (mode) {
    case Atom::ATOM:
      return static_cast<int>(typelabel_map.size()) == natomtypes;
      break;
    case Atom::BOND:
      return static_cast<int>(btypelabel_map.size()) == nbondtypes;
      break;
    case Atom::ANGLE:
      return static_cast<int>(atypelabel_map.size()) == nangletypes;
      break;
    case Atom::DIHEDRAL:
      return static_cast<int>(dtypelabel_map.size()) == ndihedraltypes;
      break;
    case Atom::IMPROPER:
      return static_cast<int>(itypelabel_map.size()) == nimpropertypes;
      break;
  }
  return false;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void LabelMap::write_data(FILE *fp)
{
  if (is_complete(Atom::ATOM)) {
    fmt::print(fp, "\nAtom Type Labels\n\n");
    for (int i = 0; i < natomtypes; i++) fmt::print(fp, "{} {}\n", i + 1, typelabel[i]);
  }

  if (force->bond && is_complete(Atom::BOND)) {
    fmt::print(fp, "\nBond Type Labels\n\n");
    for (int i = 0; i < nbondtypes; i++) fmt::print(fp, "{} {}\n", i + 1, btypelabel[i]);
  }

  if (force->angle && is_complete(Atom::ANGLE)) {
    fmt::print(fp, "\nAngle Type Labels\n\n");
    for (int i = 0; i < nangletypes; i++) fmt::print(fp, "{} {}\n", i + 1, atypelabel[i]);
  }

  if (force->dihedral && is_complete(Atom::DIHEDRAL)) {
    fmt::print(fp, "\nDihedral Type Labels\n\n");
    for (int i = 0; i < ndihedraltypes; i++) fmt::print(fp, "{} {}\n", i + 1, dtypelabel[i]);
  }

  if (force->improper && is_complete(Atom::IMPROPER)) {
    fmt::print(fp, "\nImproper Type Labels\n\n");
    for (int i = 0; i < nimpropertypes; i++) fmt::print(fp, "{} {}\n", i + 1, itypelabel[i]);
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
    typelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nbondtypes; i++) {
    charlabel = read_string(fp);
    btypelabel[i] = charlabel;
    btypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nangletypes; i++) {
    charlabel = read_string(fp);
    atypelabel[i] = charlabel;
    atypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < ndihedraltypes; i++) {
    charlabel = read_string(fp);
    dtypelabel[i] = charlabel;
    dtypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nimpropertypes; i++) {
    charlabel = read_string(fp);
    itypelabel[i] = charlabel;
    itypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void LabelMap::write_restart(FILE *fp)
{
  for (int i = 0; i < natomtypes; i++) write_string(typelabel[i], fp);

  for (int i = 0; i < nbondtypes; i++) write_string(btypelabel[i], fp);

  for (int i = 0; i < nangletypes; i++) write_string(atypelabel[i], fp);

  for (int i = 0; i < ndihedraltypes; i++) write_string(dtypelabel[i], fp);

  for (int i = 0; i < nimpropertypes; i++) write_string(itypelabel[i], fp);
}

/* ----------------------------------------------------------------------
   read a char string (including nullptr) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *LabelMap::read_string(FILE *fp)
{
  int n = read_int(fp);
  if (n < 0) error->all(FLERR, "Illegal size string or corrupt restart");
  char *value = new char[n];
  if (comm->me == 0) utils::sfread(FLERR, value, sizeof(char), n, fp, nullptr, error);
  MPI_Bcast(value, n, MPI_CHAR, 0, world);
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
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(cstr, sizeof(char), n, fp);
}

/* ----------------------------------------------------------------------
   read an int from restart file and bcast it
------------------------------------------------------------------------- */

int LabelMap::read_int(FILE *fp)
{
  int value;
  if ((comm->me == 0) && (fread(&value, sizeof(int), 1, fp) < 1)) value = -1;
  MPI_Bcast(&value, 1, MPI_INT, 0, world);
  return value;
}

/* ----------------------------------------------------------------------
   write out all current label map values as labelmap commands
------------------------------------------------------------------------- */

void LabelMap::write_map(const std::string &filename)
{
  if (comm->me == 0) {
    FILE *fp = fopen(filename.c_str(), "w");
    if (!fp) error->one(FLERR, "Cannot open label map file {}: {}", filename, utils::getsyserror());
    if (typelabel_map.size() > 0) {
      fputs("labelmap atom", fp);
      for (int i = 0; i < natomtypes; ++i)
        if (!typelabel[i].empty()) fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, typelabel[i]);
      fputc('\n', fp);
    }
    if (btypelabel_map.size() > 0) {
      fputs("labelmap bond", fp);
      for (int i = 0; i < nbondtypes; ++i)
        if (!btypelabel[i].empty()) fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, btypelabel[i]);
      fputc('\n', fp);
    }
    if (atypelabel_map.size() > 0) {
      fputs("labelmap angle", fp);
      for (int i = 0; i < nangletypes; ++i)
        if (!atypelabel[i].empty()) fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, atypelabel[i]);
      fputc('\n', fp);
    }
    if (dtypelabel_map.size() > 0) {
      fputs("labelmap dihedral", fp);
      for (int i = 0; i < ndihedraltypes; ++i)
        if (!dtypelabel[i].empty()) fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, dtypelabel[i]);
      fputc('\n', fp);
    }
    if (itypelabel_map.size() > 0) {
      fputs("labelmap improper", fp);
      for (int i = 0; i < nimpropertypes; ++i)
        if (!itypelabel[i].empty()) fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, itypelabel[i]);
      fputc('\n', fp);
    }
    fclose(fp);
  }
}
