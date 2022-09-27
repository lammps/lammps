/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FORCE_H
#define LMP_FORCE_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {
class Angle;
class Bond;
class Dihedral;
class Improper;
class KSpace;
class Pair;

enum { ENERGY_NONE = 0x00, ENERGY_GLOBAL = 0x01, ENERGY_ATOM = 0x02 };

// clang-format off
enum {
  VIRIAL_NONE     = 0x00,
  VIRIAL_PAIR     = 0x01,
  VIRIAL_FDOTR    = 0x02,
  VIRIAL_ATOM     = 0x04,
  VIRIAL_CENTROID = 0x08
};
// clang-format on

enum { CENTROID_SAME = 0, CENTROID_AVAIL = 1, CENTROID_NOTAVAIL = 2 };

class Force : protected Pointers {
 public:
  double boltz;          // Boltzmann constant (eng/degree-K)
  double hplanck;        // Planck's constant (energy-time)
  double mvv2e;          // conversion of mv^2 to energy
  double ftm2v;          // conversion of ft/m to velocity
  double mv2d;           // conversion of mass/volume to density
  double nktv2p;         // conversion of NkT/V to pressure
  double qqr2e;          // conversion of q^2/r to energy
  double qe2f;           // conversion of qE to force
  double vxmu2f;         // conversion of vx dynamic-visc to force
  double xxt2kmu;        // conversion of xx/t to kinematic-visc
  double dielectric;     // dielectric constant
  double qqrd2e;         // q^2/r to energy w/ dielectric constant
  double e_mass;         // electron mass
  double hhmrr2e;        // conversion of (hbar)^2/(mr^2) to energy
  double mvh2r;          // conversion of mv/hbar to distance
                         // hbar = h/(2*pi)
  double angstrom;       // 1 angstrom in native units
  double femtosecond;    // 1 femtosecond in native units
  double qelectron;      // 1 electron charge abs() in native units

  double qqr2e_lammps_real;    // different versions of this constant
  double qqr2e_charmm_real;    // used by new CHARMM pair styles

  int newton, newton_pair, newton_bond;    // Newton's 3rd law settings

  Pair *pair;
  char *pair_style;
  char *pair_restart;

  Bond *bond;
  char *bond_style;

  Angle *angle;
  char *angle_style;

  Dihedral *dihedral;
  char *dihedral_style;

  Improper *improper;
  char *improper_style;

  KSpace *kspace;
  char *kspace_style;

  typedef Pair *(*PairCreator)(LAMMPS *);
  typedef Bond *(*BondCreator)(LAMMPS *);
  typedef Angle *(*AngleCreator)(LAMMPS *);
  typedef Dihedral *(*DihedralCreator)(LAMMPS *);
  typedef Improper *(*ImproperCreator)(LAMMPS *);
  typedef KSpace *(*KSpaceCreator)(LAMMPS *);

  typedef std::map<std::string, PairCreator> PairCreatorMap;
  typedef std::map<std::string, BondCreator> BondCreatorMap;
  typedef std::map<std::string, AngleCreator> AngleCreatorMap;
  typedef std::map<std::string, DihedralCreator> DihedralCreatorMap;
  typedef std::map<std::string, ImproperCreator> ImproperCreatorMap;
  typedef std::map<std::string, KSpaceCreator> KSpaceCreatorMap;

  PairCreatorMap *pair_map;
  BondCreatorMap *bond_map;
  AngleCreatorMap *angle_map;
  DihedralCreatorMap *dihedral_map;
  ImproperCreatorMap *improper_map;
  KSpaceCreatorMap *kspace_map;

  // index [0] is not used in these arrays
  double special_lj[4];      // 1-2, 1-3, 1-4 prefactors for LJ
  double special_coul[4];    // 1-2, 1-3, 1-4 prefactors for Coulombics
  int special_angle;         // 0 if defined angles are ignored
                             // 1 if only weight 1,3 atoms if in an angle
  int special_dihedral;      // 0 if defined dihedrals are ignored
                             // 1 if only weight 1,4 atoms if in a dihedral
  int special_extra;         // extra space for added bonds
  int special_onefive;       // 0 if 1-5 neighbors are not stored, 1 if yes

  Force(class LAMMPS *);
  ~Force() override;
  void init();
  void setup();

  void create_pair(const std::string &, int);
  Pair *new_pair(const std::string &, int, int &);
  Pair *pair_match(const std::string &, int, int nsub = 0);
  char *pair_match_ptr(Pair *);

  void create_bond(const std::string &, int);
  Bond *new_bond(const std::string &, int, int &);
  Bond *bond_match(const std::string &);

  void create_angle(const std::string &, int);
  Angle *new_angle(const std::string &, int, int &);
  Angle *angle_match(const std::string &);

  void create_dihedral(const std::string &, int);
  Dihedral *new_dihedral(const std::string &, int, int &);
  Dihedral *dihedral_match(const std::string &);

  void create_improper(const std::string &, int);
  Improper *new_improper(const std::string &, int, int &);
  Improper *improper_match(const std::string &);

  void create_kspace(const std::string &, int);
  KSpace *new_kspace(const std::string &, int, int &);
  KSpace *kspace_match(const std::string &, int);

  char *store_style(const std::string &, int);
  void set_special(int, char **);

  double memory_usage();

 private:
  void create_factories();
};

}    // namespace LAMMPS_NS

#endif
