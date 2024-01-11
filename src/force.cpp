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

#include "force.h"

#include "style_angle.h"       // IWYU pragma: keep
#include "style_bond.h"        // IWYU pragma: keep
#include "style_dihedral.h"    // IWYU pragma: keep
#include "style_improper.h"    // IWYU pragma: keep
#include "style_kspace.h"      // IWYU pragma: keep
#include "style_pair.h"        // IWYU pragma: keep

#include "angle_hybrid.h"
#include "bond_hybrid.h"
#include "dihedral_hybrid.h"
#include "improper_hybrid.h"
#include "kspace.h"
#include "pair_hybrid.h"

#include "atom.h"
#include "comm.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

// template for factory functions:
// there will be one instance for each style keyword in the respective style_xxx.h files

template <typename S, typename T> static S *style_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

Force::Force(LAMMPS *lmp) :
    Pointers(lmp), pair(nullptr), pair_style(nullptr), pair_restart(nullptr), bond(nullptr),
    bond_style(nullptr), angle(nullptr), angle_style(nullptr), dihedral(nullptr),
    dihedral_style(nullptr), improper(nullptr), improper_style(nullptr), kspace(nullptr),
    kspace_style(nullptr)
{
  newton = newton_pair = newton_bond = 1;

  special_lj[0] = special_coul[0] = 1.0;
  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;
  special_onefive = 0;
  special_extra = 0;

  dielectric = 1.0;
  qqr2e_lammps_real = 332.06371;    // these constants are toggled
  qqr2e_charmm_real = 332.0716;     // by new CHARMM pair styles

  pair = nullptr;
  bond = nullptr;
  angle = nullptr;
  dihedral = nullptr;
  improper = nullptr;
  kspace = nullptr;

  pair_style = utils::strdup("none");
  bond_style = utils::strdup("none");
  angle_style = utils::strdup("none");
  dihedral_style = utils::strdup("none");
  improper_style = utils::strdup("none");
  kspace_style = utils::strdup("none");

  pair_restart = nullptr;
  create_factories();
}

void _noopt Force::create_factories()
{
  // fill pair map with pair styles listed in style_pair.h

  pair_map = new PairCreatorMap();

#define PAIR_CLASS
#define PairStyle(key, Class) (*pair_map)[#key] = &style_creator<Pair, Class>;
#include "style_pair.h"    // IWYU pragma: keep
#undef PairStyle
#undef PAIR_CLASS

  bond_map = new BondCreatorMap();

#define BOND_CLASS
#define BondStyle(key, Class) (*bond_map)[#key] = &style_creator<Bond, Class>;
#include "style_bond.h"    // IWYU pragma: keep
#undef BondStyle
#undef BOND_CLASS

  angle_map = new AngleCreatorMap();

#define ANGLE_CLASS
#define AngleStyle(key, Class) (*angle_map)[#key] = &style_creator<Angle, Class>;
#include "style_angle.h"    // IWYU pragma: keep
#undef AngleStyle
#undef ANGLE_CLASS

  dihedral_map = new DihedralCreatorMap();

#define DIHEDRAL_CLASS
#define DihedralStyle(key, Class) (*dihedral_map)[#key] = &style_creator<Dihedral, Class>;
#include "style_dihedral.h"    // IWYU pragma: keep
#undef DihedralStyle
#undef DIHEDRAL_CLASS

  improper_map = new ImproperCreatorMap();

#define IMPROPER_CLASS
#define ImproperStyle(key, Class) (*improper_map)[#key] = &style_creator<Improper, Class>;
#include "style_improper.h"    // IWYU pragma: keep
#undef ImproperStyle
#undef IMPROPER_CLASS

  kspace_map = new KSpaceCreatorMap();

#define KSPACE_CLASS
#define KSpaceStyle(key, Class) (*kspace_map)[#key] = &style_creator<KSpace, Class>;
#include "style_kspace.h"    // IWYU pragma: keep
#undef KSpaceStyle
#undef KSPACE_CLASS
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
  delete[] pair_style;
  delete[] bond_style;
  delete[] angle_style;
  delete[] dihedral_style;
  delete[] improper_style;
  delete[] kspace_style;

  delete[] pair_restart;

  if (pair) delete pair;
  if (bond) delete bond;
  if (angle) delete angle;
  if (dihedral) delete dihedral;
  if (improper) delete improper;
  if (kspace) delete kspace;

  pair = nullptr;
  bond = nullptr;
  angle = nullptr;
  dihedral = nullptr;
  improper = nullptr;
  kspace = nullptr;

  delete pair_map;
  delete bond_map;
  delete angle_map;
  delete dihedral_map;
  delete improper_map;
  delete kspace_map;
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
  qqrd2e = qqr2e / dielectric;

  // check if pair style must be specified after restart
  if (pair_restart) {
    if (!pair)
      error->all(FLERR, "Must re-specify non-restarted pair style ({}) after read_restart",
                 pair_restart);
  }

  if (kspace) kspace->init();    // kspace must come before pair
  if (pair) pair->init();        // so g_ewald is defined
  if (bond) bond->init();
  if (angle) angle->init();
  if (dihedral) dihedral->init();
  if (improper) improper->init();

  // print warnings if topology and force field are inconsistent

  if (comm->me == 0) {
    if (!bond && (atom->nbonds > 0)) {
      error->warning(FLERR, "Bonds are defined but no bond style is set");
      if ((special_lj[1] != 1.0) || (special_coul[1] != 1.0))
        error->warning(FLERR, "Likewise 1-2 special neighbor interactions != 1.0");
    }
    if (!angle && (atom->nangles > 0)) {
      error->warning(FLERR, "Angles are defined but no angle style is set");
      if ((special_lj[2] != 1.0) || (special_coul[2] != 1.0))
        error->warning(FLERR, "Likewise 1-3 special neighbor interactions != 1.0");
    }
    if (!dihedral && (atom->ndihedrals > 0)) {
      error->warning(FLERR, "Dihedrals are defined but no dihedral style is set");
      if ((special_lj[3] != 1.0) || (special_coul[3] != 1.0))
        error->warning(FLERR, "Likewise 1-4 special neighbor interactions != 1.0");
    }
    if (!improper && (atom->nimpropers > 0))
      error->warning(FLERR, "Impropers are defined but no improper style is set");
  }
}

/* ---------------------------------------------------------------------- */

void Force::setup()
{
  if (pair) pair->setup();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_pair(const std::string &style, int trysuffix)
{
  delete[] pair_style;
  if (pair) delete pair;
  delete[] pair_restart;
  pair_style = nullptr;
  pair = nullptr;
  pair_restart = nullptr;

  int sflag;
  pair = new_pair(style, trysuffix, sflag);
  pair_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a pair class
   if trysuffix = 1, try first with suffix1/2 appended
   return sflag = 0 for no suffix added, 1 or 2 for suffix1/2 added
------------------------------------------------------------------------- */

Pair *Force::new_pair(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator &pair_creator = (*pair_map)[estyle];
        return pair_creator(lmp);
      }
    }
    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator &pair_creator = (*pair_map)[estyle];
        return pair_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator &pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("pair", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to Pair class if matches word or matches hybrid sub-style
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   if nsub > 0, match Nth hybrid sub-style
   return nullptr if no match or if nsub=0 and multiple sub-styles match
------------------------------------------------------------------------- */

Pair *Force::pair_match(const std::string &word, int exact, int nsub)
{
  int iwhich, count;

  if (exact && (word == pair_style))
    return pair;
  else if (!exact && utils::strmatch(pair_style, word))
    return pair;
  else if (utils::strmatch(pair_style, "^hybrid")) {
    auto hybrid = dynamic_cast<PairHybrid *>(pair);
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++)
      if ((exact && (word == hybrid->keywords[i])) ||
          (!exact && utils::strmatch(hybrid->keywords[i], word))) {
        iwhich = i;
        count++;
        if (nsub == count) return hybrid->styles[iwhich];
      }
    if (count == 1) return hybrid->styles[iwhich];
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return style name of Pair class that matches Pair ptr
   called by Neighbor::print_neigh_info()
   return nullptr if no match
------------------------------------------------------------------------- */

char *Force::pair_match_ptr(Pair *ptr)
{
  if (ptr == pair) return pair_style;

  if (utils::strmatch(pair_style, "^hybrid")) {
    auto hybrid = dynamic_cast<PairHybrid *>(pair);
    for (int i = 0; i < hybrid->nstyles; i++)
      if (ptr == hybrid->styles[i]) return hybrid->keywords[i];
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   create a bond style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_bond(const std::string &style, int trysuffix)
{
  delete[] bond_style;
  if (bond) delete bond;
  bond_style = nullptr;
  bond = nullptr;

  int sflag;
  bond = new_bond(style, trysuffix, sflag);
  bond_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a bond class, fist with suffix appended
------------------------------------------------------------------------- */

Bond *Force::new_bond(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (bond_map->find(estyle) != bond_map->end()) {
        BondCreator &bond_creator = (*bond_map)[estyle];
        return bond_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (bond_map->find(estyle) != bond_map->end()) {
        BondCreator &bond_creator = (*bond_map)[estyle];
        return bond_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (bond_map->find(style) != bond_map->end()) {
    BondCreator &bond_creator = (*bond_map)[style];
    return bond_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("bond", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to current bond class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Bond *Force::bond_match(const std::string &style)
{
  if (style == bond_style)
    return bond;
  else if (strcmp(bond_style, "hybrid") == 0) {
    auto hybrid = dynamic_cast<BondHybrid *>(bond);
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   create an angle style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_angle(const std::string &style, int trysuffix)
{
  delete[] angle_style;
  if (angle) delete angle;
  angle_style = nullptr;
  angle = nullptr;

  int sflag;
  angle = new_angle(style, trysuffix, sflag);
  angle_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate an angle class
------------------------------------------------------------------------- */

Angle *Force::new_angle(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (angle_map->find(estyle) != angle_map->end()) {
        AngleCreator &angle_creator = (*angle_map)[estyle];
        return angle_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (angle_map->find(estyle) != angle_map->end()) {
        AngleCreator &angle_creator = (*angle_map)[estyle];
        return angle_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (angle_map->find(style) != angle_map->end()) {
    AngleCreator &angle_creator = (*angle_map)[style];
    return angle_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("angle", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to current angle class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Angle *Force::angle_match(const std::string &style)
{
  if (style == angle_style)
    return angle;
  else if (utils::strmatch(angle_style, "^hybrid")) {
    auto hybrid = dynamic_cast<AngleHybrid *>(angle);
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   create a dihedral style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_dihedral(const std::string &style, int trysuffix)
{
  delete[] dihedral_style;
  if (dihedral) delete dihedral;
  dihedral_style = nullptr;
  dihedral = nullptr;

  int sflag;
  dihedral = new_dihedral(style, trysuffix, sflag);
  dihedral_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a dihedral class
------------------------------------------------------------------------- */

Dihedral *Force::new_dihedral(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (dihedral_map->find(estyle) != dihedral_map->end()) {
        DihedralCreator &dihedral_creator = (*dihedral_map)[estyle];
        return dihedral_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (dihedral_map->find(estyle) != dihedral_map->end()) {
        DihedralCreator &dihedral_creator = (*dihedral_map)[estyle];
        return dihedral_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (dihedral_map->find(style) != dihedral_map->end()) {
    DihedralCreator &dihedral_creator = (*dihedral_map)[style];
    return dihedral_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("dihedral", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to current angle class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Dihedral *Force::dihedral_match(const std::string &style)
{
  if (style == dihedral_style)
    return dihedral;
  else if (utils::strmatch(dihedral_style, "^hybrid")) {
    auto hybrid = dynamic_cast<DihedralHybrid *>(dihedral);
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   create an improper style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_improper(const std::string &style, int trysuffix)
{
  delete[] improper_style;
  if (improper) delete improper;
  improper_style = nullptr;
  improper = nullptr;

  int sflag;
  improper = new_improper(style, trysuffix, sflag);
  improper_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a improper class
------------------------------------------------------------------------- */

Improper *Force::new_improper(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (improper_map->find(estyle) != improper_map->end()) {
        ImproperCreator &improper_creator = (*improper_map)[estyle];
        return improper_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (improper_map->find(estyle) != improper_map->end()) {
        ImproperCreator &improper_creator = (*improper_map)[estyle];
        return improper_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (improper_map->find(style) != improper_map->end()) {
    ImproperCreator &improper_creator = (*improper_map)[style];
    return improper_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("improper", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to current improper class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Improper *Force::improper_match(const std::string &style)
{
  if (style == improper_style)
    return improper;
  else if (utils::strmatch(improper_style, "^hybrid")) {
    auto hybrid = dynamic_cast<ImproperHybrid *>(improper);
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   new kspace style
------------------------------------------------------------------------- */

void Force::create_kspace(const std::string &style, int trysuffix)
{
  delete[] kspace_style;
  if (kspace) delete kspace;
  kspace_style = nullptr;
  kspace = nullptr;

  int sflag;
  kspace = new_kspace(style, trysuffix, sflag);
  kspace_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a kspace class
------------------------------------------------------------------------- */

KSpace *Force::new_kspace(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2 * lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (kspace_map->find(estyle) != kspace_map->end()) {
        KSpaceCreator &kspace_creator = (*kspace_map)[estyle];
        return kspace_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (kspace_map->find(estyle) != kspace_map->end()) {
        KSpaceCreator &kspace_creator = (*kspace_map)[estyle];
        return kspace_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return nullptr;
  if (kspace_map->find(style) != kspace_map->end()) {
    KSpaceCreator &kspace_creator = (*kspace_map)[style];
    return kspace_creator(lmp);
  }

  error->all(FLERR, utils::check_packages_for_style("kspace", style, lmp));

  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to Kspace class if matches word
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   return nullptr if no match
------------------------------------------------------------------------- */

KSpace *Force::kspace_match(const std::string &word, int exact)
{
  if (exact && (word == kspace_style))
    return kspace;
  else if (!exact && utils::strmatch(kspace_style, word))
    return kspace;
  return nullptr;
}

/* ----------------------------------------------------------------------
   store style name in str allocated here
   if sflag = 0, no suffix
   if sflag = 1/2, append suffix or suffix2 to style
------------------------------------------------------------------------- */

char *Force::store_style(const std::string &style, int sflag)
{
  std::string estyle = style;

  if (sflag == 1)
    estyle += std::string("/") + lmp->suffix;
  else if (sflag == 2)
    estyle += std::string("/") + lmp->suffix2;
  else if ((sflag == 3) && lmp->non_pair_suffix())
    estyle += std::string("/") + lmp->non_pair_suffix();
  return utils::strdup(estyle);
}

/* ----------------------------------------------------------------------
   set special bond values
------------------------------------------------------------------------- */

void Force::set_special(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR, "Illegal special_bonds command");

  // defaults, but do not reset special_extra

  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;
  special_onefive = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "amber") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.5;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 5.0 / 6.0;
      iarg += 1;
    } else if (strcmp(arg[iarg], "charmm") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 0.0;
      iarg += 1;
    } else if (strcmp(arg[iarg], "dreiding") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg], "fene") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 1.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 1.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg], "lj/coul") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = special_coul[1] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      special_lj[2] = special_coul[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      special_lj[3] = special_coul[3] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "lj") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_lj[1] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      special_lj[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      special_lj[3] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "coul") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_coul[1] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      special_coul[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      special_coul[3] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "angle") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_angle = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dihedral") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal special_bonds command");
      special_dihedral = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "one/five") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal special_bonds command");
      if (strcmp(arg[iarg + 1], "no") == 0)
        special_onefive = 0;
      else if (strcmp(arg[iarg + 1], "yes") == 0)
        special_onefive = 1;
      else
        error->all(FLERR, "Illegal special_bonds command");
      if (special_onefive && atom->nspecial15_flag == 0)
        error->all(FLERR,
                   "Cannot set special_bonds one/five if "
                   "atom style does not support it");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal special_bonds command");
  }

  for (int i = 1; i <= 3; i++)
    if (special_lj[i] < 0.0 || special_lj[i] > 1.0 || special_coul[i] < 0.0 ||
        special_coul[i] > 1.0)
      error->all(FLERR, "Illegal special_bonds command");
}

/* ----------------------------------------------------------------------
   memory usage of force classes
------------------------------------------------------------------------- */

double Force::memory_usage()
{
  double bytes = 0;
  if (pair) bytes += pair->memory_usage();
  if (bond) bytes += bond->memory_usage();
  if (angle) bytes += angle->memory_usage();
  if (dihedral) bytes += dihedral->memory_usage();
  if (improper) bytes += improper->memory_usage();
  if (kspace) bytes += kspace->memory_usage();
  return bytes;
}
