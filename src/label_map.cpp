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

#include "label_map.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "safe_pointers.h"
#include "tokenizer.h"

#include <algorithm>
#include <cstring>
#include <utility>

using namespace LAMMPS_NS;

static const char cite_type_label_framework[] =
    "Type Label Framework: https://doi.org/10.1021/acs.jpcb.3c08419\n\n"
    "@Article{Gissinger24,\n"
    " author = {Jacob R. Gissinger, Ilia Nikiforov, Yaser Afshar, Brendon Waters, Moon-ki Choi,"
    " Daniel S. Karls, Alexander Stukowski, Wonpil Im, Hendrik Heinz, Axel Kohlmeyer, and Ellad B. Tadmor},\n"
    " title = {Type Label Framework for Bonded Force Fields in LAMMPS},\n"
    " journal = {J. Phys. Chem. B},\n"
    " year =    2024,\n"
    " volume =  128,\n"
    " number =  13,\n"
    " pages =   {3282--3297}\n"
    "}\n\n";

namespace {
const std::string empty;
}

/* ---------------------------------------------------------------------- */

LabelMap::LabelMap(LAMMPS *_lmp, int _natomtypes, int _nbondtypes, int _nangletypes,
                   int _ndihedraltypes, int _nimpropertypes) :
    Pointers(_lmp), natomtypes(_natomtypes), nbondtypes(_nbondtypes), nangletypes(_nangletypes),
    ndihedraltypes(_ndihedraltypes), nimpropertypes(_nimpropertypes)
{
  lmap2lmap.atom = lmap2lmap.bond = lmap2lmap.angle = lmap2lmap.dihedral = lmap2lmap.improper =
      nullptr;
  checkflag = 0;
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
  for (int i = 0; i < 4; i++) check_which_labels[i] = 0;

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

  if (lmp->citeme) lmp->citeme->add(cite_type_label_framework);

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
  } else if (tlabel == "check_labels") {
    if (narg != 2)
      error->all(FLERR, "Incorrect number of arguments for labelmap check_labels command");
    for (int j = 0; j < 4; j++) check_which_labels[j] = 0;
    int i = 0;
    char option;
    while ((option = arg[1][i++]) != '\0') {
      switch (option) {
        case 'b':
          check_which_labels[0] = 1;
          break;
        case 'a':
          check_which_labels[1] = 1;
          break;
        case 'd':
          check_which_labels[2] = 1;
          break;
        case 'i':
          check_which_labels[3] = 1;
          break;
        default:
          error->all(FLERR, "Labelmap command: Illegal check_labels option {}", option);
          break;
      }
    }
    checkflag = 1;
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

  if (lmp->citeme) lmp->citeme->add(cite_type_label_framework);

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

  error->all(FLERR, "Topology type exceeds system topology type" + utils::errorurl(25));

  // never reaches here, just to prevent compiler warning

  return -1;
}

/* ----------------------------------------------------------------------
   return numeric type given a type label
   return -1 if type not yet defined
------------------------------------------------------------------------- */

int LabelMap::find_type(const std::string &mylabel, int mode) const
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
   return type label given a numeric type
   return "" if type label does not exist
------------------------------------------------------------------------- */

const std::string &LabelMap::find_label(int i, int mode) const
{
  switch (mode) {
    case Atom::ATOM:
      if ((i > 0) && (i <= atom->ntypes)) {
        if (is_complete(mode)) return typelabel[i - 1];
      }
      break;
    case Atom::BOND:
      if ((i > 0) && (i <= atom->nbondtypes)) {
        if (is_complete(mode)) return btypelabel[i - 1];
      }
      break;
    case Atom::ANGLE:
      if ((i > 0) && (i <= atom->nangletypes)) {
        if (is_complete(mode)) return atypelabel[i - 1];
      }
      break;
    case Atom::DIHEDRAL:
      if ((i > 0) && (i <= atom->ndihedraltypes)) {
        if (is_complete(mode)) return dtypelabel[i - 1];
      }
      break;
    case Atom::IMPROPER:
      if ((i > 0) && (i <= atom->nimpropertypes)) {
        if (is_complete(mode)) return itypelabel[i - 1];
      }
      break;
    default:
      return empty;
  }
  return empty;
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
   infer bond type from two atom types
   input/output is numeric types, uses type labels internally
   assumes bond type labels are of the form "a-b" for atom types 'a' and 'b'
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_bondtype(int type1, int type2)
{
  // check for out of range input
  if ((type1 < 1) || (type1 > natomtypes) || (type2 < 1) || (type2 > natomtypes)) return 0;

  // convert numeric atom types to type label
  std::vector<std::string> mytypes(2);
  mytypes[0] = typelabel[type1 - 1];
  mytypes[1] = typelabel[type2 - 1];
  if (mytypes[0].empty() || mytypes[1].empty()) return 0;

  return infer_bondtype(mytypes);
}

/* ----------------------------------------------------------------------
   infer numeric type from two atom type labels
   assumes bond types are of the form "a-b" for atom types 'a' and 'b'
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_bondtype(const std::vector<std::string> &mytypes)
{
  // search for matching bond type label with symmetry considerations
  int out = 0;
  std::vector<std::string> btypes(2);
  for (int i = 0; i < nbondtypes; i++) {
    int status = parse_typelabel(2, btypelabel[i], btypes);
    if ((status != -1) && (btypes.size() == 2)) {
      if (mytypes[0] == btypes[0] && mytypes[1] == btypes[1]) return i + 1;
      if (mytypes[0] == btypes[1] && mytypes[1] == btypes[0]) out = -(i + 1);
    }
  }
  return out;
}

/* ----------------------------------------------------------------------
   infer angle type from three atom types
   input/output is numeric types, uses type labels internally
   assumes angle types of the form "a-b-c" for atom types 'a', 'b', 'c'
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_angletype(int type1, int type2, int type3)
{
  // check for out of range input
  if ((type1 < 1) || (type1 > natomtypes) || (type2 < 1) || (type2 > natomtypes) || (type3 < 1) ||
      (type3 > natomtypes))
    return 0;

  // convert numeric atom types to type label
  std::vector<std::string> mytypes(3);
  mytypes[0] = typelabel[type1 - 1];
  mytypes[1] = typelabel[type2 - 1];
  mytypes[2] = typelabel[type3 - 1];
  for (size_t i = 0; i < 3; i++)
    if (mytypes[i].empty()) return 0;

  return infer_angletype(mytypes);
}

/* ----------------------------------------------------------------------
   infer angle type from three atom types
   input/output is numeric types, uses type labels internally
   assumes angle types of the form "a-b-c" for atom types 'a', 'b', 'c'
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_angletype(const std::vector<std::string> &mytypes)
{
  // search for matching angle type label, with symmetry considerations

  int out = 0;
  int status;
  std::vector<std::string> atypes(3);
  for (int i = 0; i < nangletypes; i++) {
    status = parse_typelabel(3, atypelabel[i], atypes);
    if (status != -1 && mytypes[1] == atypes[1]) {
      if (mytypes[0] == atypes[0] && mytypes[2] == atypes[2]) return i + 1;
      if (mytypes[0] == atypes[2] && mytypes[2] == atypes[0]) out = -(i + 1);
    }
  }
  return out;
}

/* ----------------------------------------------------------------------
   infer dihedral type from four atom types
   input/output is numeric types, uses type labels internally
   assumes dihedral types of the form "a-b-c-d"
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_dihedraltype(int type1, int type2, int type3, int type4)
{
  // check for out of range input
  if ((type1 < 1) || (type1 > natomtypes) || (type2 < 1) || (type2 > natomtypes) || (type3 < 1) ||
      (type3 > natomtypes) || (type4 < 1) || (type4 > natomtypes))
    return 0;

  // convert numeric atom types to type label
  std::vector<std::string> mytypes(4);
  mytypes[0] = typelabel[type1 - 1];
  mytypes[1] = typelabel[type2 - 1];
  mytypes[2] = typelabel[type3 - 1];
  mytypes[3] = typelabel[type4 - 1];
  for (size_t i = 0; i < 4; i++)
    if (mytypes[i].empty()) return 0;

  return infer_dihedraltype(mytypes);
}

/* ----------------------------------------------------------------------
   infer dihedral type from four atom types
   input/output is numeric types, uses type labels internally
   assumes dihedral types of the form "a-b-c-d"
   returns negative of numeric type if constituent atoms types in reverse order
------------------------------------------------------------------------- */

int LabelMap::infer_dihedraltype(const std::vector<std::string> &mytypes)
{
  // search for matching dihedral type label

  int out = 0;
  int status;
  std::vector<std::string> dtypes(4);
  for (int i = 0; i < ndihedraltypes; i++) {
    status = parse_typelabel(4, dtypelabel[i], dtypes);
    if (status != -1) {
      if (mytypes[0] == dtypes[0] && mytypes[1] == dtypes[1] && mytypes[2] == dtypes[2] &&
          mytypes[3] == dtypes[3])
        return i + 1;
      if (mytypes[3] == dtypes[0] && mytypes[2] == dtypes[1] && mytypes[1] == dtypes[2] &&
          mytypes[0] == dtypes[3])
        out = -(i + 1);
    }
  }
  return out;
}

/* ----------------------------------------------------------------------
   infer improper type from four atom types
   input/output is numeric types, uses type labels internally
   assumes improper types of the form "a-b-c-d"
   the symmetry of the improper is encoded in improper.symmatoms
------------------------------------------------------------------------- */

int LabelMap::infer_impropertype(int type1, int type2, int type3, int type4, std::array<int, 4> *iorder)
{
  // check for out of range input
  if ((type1 < 1) || (type1 > natomtypes) || (type2 < 1) || (type2 > natomtypes) || (type3 < 1) ||
      (type3 > natomtypes) || (type4 < 1) || (type4 > natomtypes))
    return 0;

  // convert numeric atom types to type label
  std::vector<std::string> mytypes(4);
  mytypes[0] = typelabel[type1 - 1];
  mytypes[1] = typelabel[type2 - 1];
  mytypes[2] = typelabel[type3 - 1];
  mytypes[3] = typelabel[type4 - 1];
  for (int i = 0; i < 4; i++)
    if (mytypes[i].empty()) return 0;

  return infer_impropertype(mytypes, iorder);
}

/* ----------------------------------------------------------------------
   infer improper type from four atom types
   input/output is numeric types, uses type labels internally
   assumes improper types of the form "a-b-c-d"
   the symmetry of the improper is encoded in improper.symmatoms
------------------------------------------------------------------------- */

int LabelMap::infer_impropertype(const std::vector<std::string> &mytypes, std::array<int, 4> *iorder)
{
  // search for matching improper type label
  int out = 0;
  int status, navail_types;
  std::vector<std::string> itypes(4);
  std::vector<std::string> avail_types;
  for (int i = 0; i < nimpropertypes; i++) {
    status = parse_typelabel(4, itypelabel[i], itypes);
    if (status != -1) {
      if (mytypes[0] == itypes[0] && mytypes[1] == itypes[1] && mytypes[2] == itypes[2] &&
          mytypes[3] == itypes[3])
        return i + 1;
      navail_types = 4;
      avail_types = mytypes;
      for (int j = 0; j < 4; j++) {
        if (force->improper && force->improper->symmatoms[j] == 1) {
          if (mytypes[j] != itypes[j]) {
            status = -1;
            break;
          }
          avail_types[j] = "";
          navail_types--;
        }
      }
      if (status == -1) continue;

      if (iorder) *iorder = {0, 1, 2, 3};
      for (int j = 0; j < 4; j++) {
        if (std::string(force->improper_style) == "none" || force->improper->symmatoms[j] != 1) {
          for (int k = 0; k < 4; k++) {
            if (itypes[j] == avail_types[k]) {
              avail_types[k] = "";
              navail_types--;
              if (iorder) (*iorder)[j] = k;
              break;
            }
          }
        }
      }
      if (navail_types == 0) out = -(i + 1);
    }
  }
  return out;
}

/* ----------------------------------------------------------------------
   return -1 if number of parsed strings is not equal to ntypes input
------------------------------------------------------------------------- */

int LabelMap::parse_typelabel(int ntypes, const std::string &label, std::vector<std::string> &types)
{
  auto out = Tokenizer(label, "-").as_vector();
  if ((int) out.size() != ntypes) return -1;
  types = std::move(out);
  return 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void LabelMap::write_data(FILE *fp)
{
  if (is_complete(Atom::ATOM)) {
    utils::print(fp, "\nAtom Type Labels\n\n");
    for (int i = 0; i < natomtypes; i++) utils::print(fp, "{} {}\n", i + 1, typelabel[i]);
  }

  if (force->bond && is_complete(Atom::BOND)) {
    utils::print(fp, "\nBond Type Labels\n\n");
    for (int i = 0; i < nbondtypes; i++) utils::print(fp, "{} {}\n", i + 1, btypelabel[i]);
  }

  if (force->angle && is_complete(Atom::ANGLE)) {
    utils::print(fp, "\nAngle Type Labels\n\n");
    for (int i = 0; i < nangletypes; i++) utils::print(fp, "{} {}\n", i + 1, atypelabel[i]);
  }

  if (force->dihedral && is_complete(Atom::DIHEDRAL)) {
    utils::print(fp, "\nDihedral Type Labels\n\n");
    for (int i = 0; i < ndihedraltypes; i++) utils::print(fp, "{} {}\n", i + 1, dtypelabel[i]);
  }

  if (force->improper && is_complete(Atom::IMPROPER)) {
    utils::print(fp, "\nImproper Type Labels\n\n");
    for (int i = 0; i < nimpropertypes; i++) utils::print(fp, "{} {}\n", i + 1, itypelabel[i]);
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
    if (strlen(charlabel) > 0) typelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nbondtypes; i++) {
    charlabel = read_string(fp);
    btypelabel[i] = charlabel;
    if (strlen(charlabel) > 0) btypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nangletypes; i++) {
    charlabel = read_string(fp);
    atypelabel[i] = charlabel;
    if (strlen(charlabel) > 0) atypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < ndihedraltypes; i++) {
    charlabel = read_string(fp);
    dtypelabel[i] = charlabel;
    if (strlen(charlabel) > 0) dtypelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }

  for (int i = 0; i < nimpropertypes; i++) {
    charlabel = read_string(fp);
    itypelabel[i] = charlabel;
    if (strlen(charlabel) > 0) itypelabel_map[charlabel] = i + 1;
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
    SafeFilePtr fp = fopen(filename.c_str(), "w");
    if (!fp) error->one(FLERR, "Cannot open label map file {}: {}", filename, utils::getsyserror());
    if (typelabel_map.size() > 0) {
      fputs("labelmap atom", fp);
      for (int i = 0; i < natomtypes; ++i)
        if (!typelabel[i].empty()) utils::print(fp, R"( {} """ {} """)", i + 1, typelabel[i]);
      fputc('\n', fp);
    }
    if (btypelabel_map.size() > 0) {
      fputs("labelmap bond", fp);
      for (int i = 0; i < nbondtypes; ++i)
        if (!btypelabel[i].empty()) utils::print(fp, R"( {} """ {} """)", i + 1, btypelabel[i]);
      fputc('\n', fp);
    }
    if (atypelabel_map.size() > 0) {
      fputs("labelmap angle", fp);
      for (int i = 0; i < nangletypes; ++i)
        if (!atypelabel[i].empty()) utils::print(fp, R"( {} """ {} """)", i + 1, atypelabel[i]);
      fputc('\n', fp);
    }
    if (dtypelabel_map.size() > 0) {
      fputs("labelmap dihedral", fp);
      for (int i = 0; i < ndihedraltypes; ++i)
        if (!dtypelabel[i].empty()) utils::print(fp, R"( {} """ {} """)", i + 1, dtypelabel[i]);
      fputc('\n', fp);
    }
    if (itypelabel_map.size() > 0) {
      fputs("labelmap improper", fp);
      for (int i = 0; i < nimpropertypes; ++i)
        if (!itypelabel[i].empty()) utils::print(fp, R"( {} """ {} """)", i + 1, itypelabel[i]);
      fputc('\n', fp);
    }
  }
}

/* ----------------------------------------------------------------------
   check type label self-consistency
------------------------------------------------------------------------- */

void LabelMap::check_labels()
{
  int *type = atom->type;
  tagint *tag = atom->tag;
  // in rare cases, bonds are not symmetric. only check if newton on for bonds
  int globally_perfect_labels;
  int perfect_labels = 1;
  if (force->newton_bond && check_which_labels[0]) {
    for (int i = 0; i < atom->nlocal; i++) {
      int atom1 = i;
      for (int j = 0; j < atom->num_bond[i]; j++) {
        int btype = atom->bond_type[i][j];
        if (btype < 1) continue;
        int atom2 = atom->map(atom->bond_atom[i][j]);
        int inferred_type = atom->lmap->infer_bondtype(type[atom1], type[atom2]);
        if (inferred_type != btype) {
          perfect_labels = 0;
          std::string atom1_label = atom->lmap->find_label(type[atom1], Atom::ATOM);
          std::string atom2_label = atom->lmap->find_label(type[atom2], Atom::ATOM);
          std::string blabel = atom->lmap->find_label(btype, Atom::BOND);
          if (inferred_type == -btype)
            error->warning(FLERR,
                           "Bond between atoms {}, {} has constituent atom types ({}, {}) in "
                           "reverse order compared "
                           "to its bond type label ({})",
                           tag[atom1], tag[atom2], atom1_label, atom2_label, blabel);
          else
            error->warning(
                FLERR,
                "Bond between atoms {}, {} has constituent atom types ({}, {}) that do not match "
                "its type label ({})",
                tag[atom1], tag[atom2], atom1_label, atom2_label, blabel);
        }
      }
    }
    MPI_Reduce(&perfect_labels, &globally_perfect_labels, 1, MPI_INT, MPI_SUM, 0, world);
    if (comm->me == 0 && globally_perfect_labels == comm->nprocs)
      utils::logmesg(lmp, "All bonds in the simulation have self-consistent type labels\n");
  }

  // some angles are not symmetric, like class2
  perfect_labels = 1;
  if (check_which_labels[1]) {
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = 0; j < atom->num_angle[i]; j++) {
        int atype = atom->angle_type[i][j];
        if (atype < 1) continue;
        int atom1 = atom->map(atom->angle_atom1[i][j]);
        int atom2 = atom->map(atom->angle_atom2[i][j]);
        int atom3 = atom->map(atom->angle_atom3[i][j]);
        int inferred_type = atom->lmap->infer_angletype(type[atom1], type[atom2], type[atom3]);
        if (inferred_type != atype) {
          perfect_labels = 0;
          std::string atom1_label = atom->lmap->find_label(type[atom1], Atom::ATOM);
          std::string atom2_label = atom->lmap->find_label(type[atom2], Atom::ATOM);
          std::string atom3_label = atom->lmap->find_label(type[atom3], Atom::ATOM);
          std::string alabel = atom->lmap->find_label(atype, Atom::ANGLE);
          if (inferred_type == -atype)
            error->warning(FLERR,
                           "Angle between atoms {}, {}, {} has constituent atom types ({}, {}, {}) "
                           "in reverse order compared "
                           "to its angle type label ({})",
                           tag[atom1], tag[atom2], tag[atom3], atom1_label, atom2_label,
                           atom3_label, alabel);
          else
            error->warning(FLERR,
                           "Angle between atoms {}, {}, {} has constituent atom types ({}, {}, {}) "
                           "that do not match its "
                           "type label ({})",
                           tag[atom1], tag[atom2], tag[atom3], atom1_label, atom2_label,
                           atom3_label, alabel);
        }
      }
    }
    MPI_Reduce(&perfect_labels, &globally_perfect_labels, 1, MPI_INT, MPI_SUM, 0, world);
    if (comm->me == 0 && globally_perfect_labels == comm->nprocs)
      utils::logmesg(lmp, "All angles in the simulation have self-consistent type labels\n");
  }

  // some dihedrals are not symmetric, like class2
  perfect_labels = 1;
  if (check_which_labels[2]) {
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = 0; j < atom->num_dihedral[i]; j++) {
        int dtype = atom->dihedral_type[i][j];
        if (dtype < 1) continue;
        int atom1 = atom->map(atom->dihedral_atom1[i][j]);
        int atom2 = atom->map(atom->dihedral_atom2[i][j]);
        int atom3 = atom->map(atom->dihedral_atom3[i][j]);
        int atom4 = atom->map(atom->dihedral_atom4[i][j]);
        int inferred_type =
            atom->lmap->infer_dihedraltype(type[atom1], type[atom2], type[atom3], type[atom4]);
        if (inferred_type != dtype) {
          perfect_labels = 0;
          std::string atom1_label = atom->lmap->find_label(type[atom1], Atom::ATOM);
          std::string atom2_label = atom->lmap->find_label(type[atom2], Atom::ATOM);
          std::string atom3_label = atom->lmap->find_label(type[atom3], Atom::ATOM);
          std::string atom4_label = atom->lmap->find_label(type[atom4], Atom::ATOM);
          std::string dlabel = atom->lmap->find_label(dtype, Atom::DIHEDRAL);
          if (inferred_type == -dtype)
            error->warning(FLERR,
                           "Dihedral between atoms {}, {}, {}, {} has constituent atom types ({}, "
                           "{}, {}, {}) in reverse order compared to its "
                           "dihedral type label ({})",
                           tag[atom1], tag[atom2], tag[atom3], tag[atom4], atom1_label, atom2_label,
                           atom3_label, atom4_label, dlabel);
          else
            error->warning(FLERR,
                           "Dihedral between atoms {}, {}, {}, {} has constituent atom types ({}, "
                           "{}, {}, {}) that do not match its "
                           "dihedral label ({})",
                           tag[atom1], tag[atom2], tag[atom3], tag[atom4], atom1_label, atom2_label,
                           atom3_label, atom4_label, dlabel);
        }
      }
    }
    MPI_Reduce(&perfect_labels, &globally_perfect_labels, 1, MPI_INT, MPI_SUM, 0, world);
    if (comm->me == 0 && globally_perfect_labels == comm->nprocs)
      utils::logmesg(lmp, "All dihedrals in the simulation have self-consistent type labels\n");
  }

  // some impropers are not symmetric, like class2
  perfect_labels = 1;
  if (check_which_labels[3]) {
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = 0; j < atom->num_improper[i]; j++) {
        int itype = atom->improper_type[i][j];
        if (itype < 1) continue;
        int atom1 = atom->map(atom->improper_atom1[i][j]);
        int atom2 = atom->map(atom->improper_atom2[i][j]);
        int atom3 = atom->map(atom->improper_atom3[i][j]);
        int atom4 = atom->map(atom->improper_atom4[i][j]);
        int inferred_type =
            atom->lmap->infer_impropertype(type[atom1], type[atom2], type[atom3], type[atom4]);
        if (inferred_type != itype) {
          perfect_labels = 0;
          std::string atom1_label = atom->lmap->find_label(type[atom1], Atom::ATOM);
          std::string atom2_label = atom->lmap->find_label(type[atom2], Atom::ATOM);
          std::string atom3_label = atom->lmap->find_label(type[atom3], Atom::ATOM);
          std::string atom4_label = atom->lmap->find_label(type[atom4], Atom::ATOM);
          std::string ilabel = atom->lmap->find_label(itype, Atom::IMPROPER);
          if (inferred_type == -itype)
            error->warning(FLERR,
                           "Improper containing atoms {}, {}, {}, {} has constituent atom types "
                           "({}, {}, {}, {}) in a different order compared to its "
                           "improper type label ({})",
                           tag[atom1], tag[atom2], tag[atom3], tag[atom4], atom1_label, atom2_label,
                           atom3_label, atom4_label, ilabel);
          else
            error->warning(FLERR,
                           "Improper containing atoms {}, {}, {}, {} has constituent atom types "
                           "({}, {}, {}, {}) that do not match its "
                           "improper label ({})",
                           tag[atom1], tag[atom2], tag[atom3], tag[atom4], atom1_label, atom2_label,
                           atom3_label, atom4_label, ilabel);
        }
      }
    }
    MPI_Reduce(&perfect_labels, &globally_perfect_labels, 1, MPI_INT, MPI_SUM, 0, world);
    if (comm->me == 0 && globally_perfect_labels == comm->nprocs)
      utils::logmesg(lmp, "All impropers in the simulation have self-consistent type labels\n");
  }
  checkflag = 0;
}
