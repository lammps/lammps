/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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
#include <cstdio>
#include <map>
#include <string>

namespace LAMMPS_NS {

class Force : protected Pointers {
 public:
  double boltz;                      // Boltzmann constant (eng/degree-K)
  double hplanck;                    // Planck's constant (energy-time)
  double mvv2e;                      // conversion of mv^2 to energy
  double ftm2v;                      // conversion of ft/m to velocity
  double mv2d;                       // conversion of mass/volume to density
  double nktv2p;                     // conversion of NkT/V to pressure
  double qqr2e;                      // conversion of q^2/r to energy
  double qe2f;                       // conversion of qE to force
  double vxmu2f;                     // conversion of vx dynamic-visc to force
  double xxt2kmu;                    // conversion of xx/t to kinematic-visc
  double dielectric;                 // dielectric constant
  double qqrd2e;                     // q^2/r to energy w/ dielectric constant
  double e_mass;                     // electron mass
  double hhmrr2e;                    // conversion of (hbar)^2/(mr^2) to energy
  double mvh2r;                      // conversion of mv/hbar to distance
                                     // hbar = h/(2*pi)
  double angstrom;                   // 1 angstrom in native units
  double femtosecond;                // 1 femtosecond in native units
  double qelectron;                  // 1 electron charge abs() in native units

  double qqr2e_lammps_real;          // different versions of this constant
  double qqr2e_charmm_real;          // used by new CHARMM pair styles

  int newton,newton_pair,newton_bond;   // Newton's 3rd law settings

  class Pair *pair;
  char *pair_style;
  char *pair_restart;

  class Bond *bond;
  char *bond_style;

  class Angle *angle;
  char *angle_style;

  class Dihedral *dihedral;
  char *dihedral_style;

  class Improper *improper;
  char *improper_style;

  class KSpace *kspace;
  char *kspace_style;

  typedef Pair *(*PairCreator)(LAMMPS *);
  typedef Bond *(*BondCreator)(LAMMPS *);
  typedef Angle *(*AngleCreator)(LAMMPS *);
  typedef Dihedral *(*DihedralCreator)(LAMMPS *);
  typedef Improper *(*ImproperCreator)(LAMMPS *);
  typedef KSpace *(*KSpaceCreator)(LAMMPS *);

  typedef std::map<std::string,PairCreator> PairCreatorMap;
  typedef std::map<std::string,BondCreator> BondCreatorMap;
  typedef std::map<std::string,AngleCreator> AngleCreatorMap;
  typedef std::map<std::string,DihedralCreator> DihedralCreatorMap;
  typedef std::map<std::string,ImproperCreator> ImproperCreatorMap;
  typedef std::map<std::string,KSpaceCreator> KSpaceCreatorMap;

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

  void create_pair(const char *, int);
  class Pair *new_pair(const char *, int, int &);
  class Pair *pair_match(const char *, int, int nsub=0);
  char *pair_match_ptr(Pair *);

  void create_bond(const char *, int);
  class Bond *new_bond(const char *, int, int &);
  class Bond *bond_match(const char *);

  void create_angle(const char *, int);
  class Angle *new_angle(const char *, int, int &);
  class Angle *angle_match(const char *);

  void create_dihedral(const char *, int);
  class Dihedral *new_dihedral(const char *, int, int &);
  class Dihedral *dihedral_match(const char *);

  void create_improper(const char *, int);
  class Improper *new_improper(const char *, int, int &);
  class Improper *improper_match(const char *);

  void create_kspace(const char *, int);
  class KSpace *new_kspace(const char *, int, int &);
  class KSpace *kspace_match(const char *, int);

  void store_style(char *&, const char *, int);
  void set_special(int, char **);
  void bounds(const char *, int, char *, int, int &, int &, int nmin=1);
  void boundsbig(const char *, int, char *, bigint, bigint &, bigint &, bigint nmin=1);
  double numeric(const char *, int, char *);
  int inumeric(const char *, int, char *);
  bigint bnumeric(const char *, int, char *);
  tagint tnumeric(const char *, int, char *);

  FILE *open_potential(const char *);
  const char *potential_name(const char *);
  void potential_date(FILE *, const char *);

  bigint memory_usage();

 private:
  template <typename T> static Pair *pair_creator(LAMMPS *);
  template <typename T> static Bond *bond_creator(LAMMPS *);
  template <typename T> static Angle *angle_creator(LAMMPS *);
  template <typename T> static Dihedral *dihedral_creator(LAMMPS *);
  template <typename T> static Improper *improper_creator(LAMMPS *);
  template <typename T> static KSpace *kspace_creator(LAMMPS *);
};

}

#endif

/* ERROR/WARNING messages:

E: Must re-specify non-restarted pair style (%s) after read_restart

UNDOCUMENTED

E: Unknown pair style %s

UNDOCUMENTED

E: Unknown bond style %s

UNDOCUMENTED

E: Unknown angle style %s

UNDOCUMENTED

E: Unknown dihedral style %s

UNDOCUMENTED

E: Unknown improper style %s

UNDOCUMENTED

E: Cannot yet use KSpace solver with grid with comm style tiled

This is current restriction in LAMMPS.

E: Unknown kspace style %s

UNDOCUMENTED

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

*/
