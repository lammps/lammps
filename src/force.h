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

  Force(class LAMMPS *);
  ~Force();
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

  void store_style(char *&, const std::string &, int);
  void set_special(int, char **);

  double memory_usage();

 private:
  void create_factories();
  template <typename T> static Pair *pair_creator(LAMMPS *);
  template <typename T> static Bond *bond_creator(LAMMPS *);
  template <typename T> static Angle *angle_creator(LAMMPS *);
  template <typename T> static Dihedral *dihedral_creator(LAMMPS *);
  template <typename T> static Improper *improper_creator(LAMMPS *);
  template <typename T> static KSpace *kspace_creator(LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Must re-specify non-restarted pair style (%s) after read_restart

For pair styles, that do not store their settings in a restart file,
it must be defined with a new 'pair_style' command after read_restart.

E: Unrecognized pair style %s

The choice of pair style is unknown.

E: Unrecognized bond style %s

The choice of bond style is unknown.

E: Unrecognized angle style %s

The choice of angle style is unknown.

E: Unrecognized dihedral style %s

The choice of dihedral style is unknown.

E: Unrecognized improper style %s

The choice of improper style is unknown.

E: Cannot yet use KSpace solver with grid with comm style tiled

This is current restriction in LAMMPS.

E: Unrecognized kspace style %s

The choice of kspace style is unknown.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Unknown pair style

The choice of pair style is unknown.

U: Unknown bond style

The choice of bond style is unknown.

U: Unknown angle style

The choice of angle style is unknown.

U: Unknown dihedral style

The choice of dihedral style is unknown.

U: Unknown improper style

The choice of improper style is unknown.

U: Unknown kspace style

The choice of kspace style is unknown.

U: Numeric index is out of bounds

A command with an argument that specifies an integer or range of
integers is using a value that is less than 1 or greater than the
maximum allowed limit.

W: Bonds are defined but no bond style is set

The topology contains bonds, but there are no bond forces computed
since there was no bond_style command.

W: Angles are defined but no angle style is set

The topology contains angles, but there are no angle forces computed
since there was no angle_style command.

W: Dihedrals are defined but no dihedral style is set

The topology contains dihedrals, but there are no dihedral forces computed
since there was no dihedral_style command.

W: Impropers are defined but no improper style is set

The topology contains impropers, but there are no improper forces computed
since there was no improper_style command.

W: Likewise 1-2 special neighbor interactions != 1.0

The topology contains bonds, but there is no bond style defined
and a 1-2 special neighbor scaling factor was not 1.0. This
means that pair style interactions may have scaled or missing
pairs in the neighbor list in expectation of interactions for
those pairs being computed from the bond style.

W: Likewise 1-3 special neighbor interactions != 1.0

The topology contains angles, but there is no angle style defined
and a 1-3 special neighbor scaling factor was not 1.0. This
means that pair style interactions may have scaled or missing
pairs in the neighbor list in expectation of interactions for
those pairs being computed from the angle style.

W: Likewise 1-4 special neighbor interactions != 1.0

The topology contains dihedrals, but there is no dihedral style defined
and a 1-4 special neighbor scaling factor was not 1.0. This
means that pair style interactions may have scaled or missing
pairs in the neighbor list in expectation of interactions for
those pairs being computed from the dihedral style.

*/
