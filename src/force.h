/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

/** Flags for energy computation modes */
enum { ENERGY_NONE = 0x00, ENERGY_GLOBAL = 0x01, ENERGY_ATOM = 0x02 };

// clang-format off
/** Flags for virial computation modes */
enum {
  VIRIAL_NONE     = 0x00,    /**< No virial computation */
  VIRIAL_PAIR     = 0x01,    /**< Virial from pair interactions */
  VIRIAL_FDOTR    = 0x02,    /**< Virial from f dot r */
  VIRIAL_ATOM     = 0x04,    /**< Per-atom virial */
  VIRIAL_CENTROID = 0x08     /**< Centroid virial for rigid bodies */
};
// clang-format on

/** Flags for centroid stress availability */
enum { CENTROID_SAME = 0, CENTROID_AVAIL = 1, CENTROID_NOTAVAIL = 2 };

/** \class LAMMPS_NS::Force
 * \brief Manager for interatomic force field styles in LAMMPS
 *
 * The Force class manages all force field components in a LAMMPS simulation,
 * including pair potentials, bond/angle/dihedral/improper interactions, and
 * long-range electrostatics (KSpace). It serves as a container and factory
 * for these interaction styles.
 *
 * Key responsibilities include:
 *
 * - Creating and managing pair, bond, angle, dihedral, improper, and kspace styles
 * - Storing physical constants for unit conversions
 * - Managing Newton's third law settings for pairwise and bonded interactions
 * - Storing special bonds settings for 1-2, 1-3, 1-4 neighbor exclusions
 * - Providing style factories for dynamic creation of force field components
 *
 * The class maintains pointers to the currently active force field styles
 * and provides methods to create, match, and query these styles.
 *
 * \sa LAMMPS_NS::Pair, LAMMPS_NS::Bond, LAMMPS_NS::Angle, LAMMPS_NS::KSpace */
class Force : protected Pointers {
 public:
  double boltz;          /**< Boltzmann constant (energy/degree-K) */
  double hplanck;        /**< Planck's constant (energy-time) */
  double mvv2e;          /**< Conversion factor: mv^2 to energy */
  double ftm2v;          /**< Conversion factor: ft/m to velocity */
  double mv2d;           /**< Conversion factor: mass/volume to density */
  double nktv2p;         /**< Conversion factor: NkT/V to pressure */
  double qqr2e;          /**< Conversion factor: q^2/r to energy */
  double qe2f;           /**< Conversion factor: qE to force */
  double vxmu2f;         /**< Conversion factor: v*dynamic-viscosity to force */
  double xxt2kmu;        /**< Conversion factor: xx/t to kinematic-viscosity */
  double dielectric;     /**< Dielectric constant */
  double qqrd2e;         /**< Conversion factor: q^2/r to energy with dielectric */
  double e_mass;         /**< Electron mass */
  double hhmrr2e;        /**< Conversion factor: (hbar)^2/(mr^2) to energy */
  double mvh2r;          /**< Conversion factor: mv/hbar to distance */
  double angstrom;       /**< 1 angstrom in native units */
  double femtosecond;    /**< 1 femtosecond in native units */
  double qelectron;      /**< 1 electron charge (absolute) in native units */

  double qqr2e_lammps_real;    /**< LAMMPS version of qqr2e for real units */
  double qqr2e_charmm_real;    /**< CHARMM version of qqr2e for real units */

  int newton;          /**< Newton's 3rd law setting: 0 = off, 1 = on */
  int newton_pair;     /**< Newton's 3rd law for pair: 0 = off, 1 = on */
  int newton_bond;     /**< Newton's 3rd law for bonds: 0 = off, 1 = on */

  Pair *pair;              /**< Pointer to current pair style instance */
  char *pair_style;        /**< Current pair style name */
  char *pair_restart;      /**< Pair style name from restart file */

  Bond *bond;              /**< Pointer to current bond style instance */
  char *bond_style;        /**< Current bond style name */

  Angle *angle;            /**< Pointer to current angle style instance */
  char *angle_style;       /**< Current angle style name */

  Dihedral *dihedral;      /**< Pointer to current dihedral style instance */
  char *dihedral_style;    /**< Current dihedral style name */

  Improper *improper;      /**< Pointer to current improper style instance */
  char *improper_style;    /**< Current improper style name */

  KSpace *kspace;          /**< Pointer to current kspace style instance */
  char *kspace_style;      /**< Current kspace style name */

  using PairCreator = Pair *(*)(LAMMPS *);            /**< Function pointer type for pair creation */
  using BondCreator = Bond *(*)(LAMMPS *);            /**< Function pointer type for bond creation */
  using AngleCreator = Angle *(*)(LAMMPS *);          /**< Function pointer type for angle creation */
  using DihedralCreator = Dihedral *(*)(LAMMPS *);    /**< Function pointer type for dihedral creation */
  using ImproperCreator = Improper *(*)(LAMMPS *);    /**< Function pointer type for improper creation */
  using KSpaceCreator = KSpace *(*)(LAMMPS *);        /**< Function pointer type for kspace creation */

  using PairCreatorMap = std::map<std::string, PairCreator>;          /**< Map type for pair style factory */
  using BondCreatorMap = std::map<std::string, BondCreator>;          /**< Map type for bond style factory */
  using AngleCreatorMap = std::map<std::string, AngleCreator>;        /**< Map type for angle style factory */
  using DihedralCreatorMap = std::map<std::string, DihedralCreator>;  /**< Map type for dihedral style factory */
  using ImproperCreatorMap = std::map<std::string, ImproperCreator>;  /**< Map type for improper style factory */
  using KSpaceCreatorMap = std::map<std::string, KSpaceCreator>;      /**< Map type for kspace style factory */

  PairCreatorMap *pair_map;          /**< Factory map for pair styles */
  BondCreatorMap *bond_map;          /**< Factory map for bond styles */
  AngleCreatorMap *angle_map;        /**< Factory map for angle styles */
  DihedralCreatorMap *dihedral_map;  /**< Factory map for dihedral styles */
  ImproperCreatorMap *improper_map;  /**< Factory map for improper styles */
  KSpaceCreatorMap *kspace_map;      /**< Factory map for kspace styles */

  // index [0] is not used in these arrays
  double special_lj[4];      /**< 1-2, 1-3, 1-4 weighting factors for LJ interactions */
  double special_coul[4];    /**< 1-2, 1-3, 1-4 weighting factors for Coulombic interactions */
  int special_angle;         /**< 0 = ignore defined angles, 1 = only weight 1,3 if in angle */
  int special_dihedral;      /**< 0 = ignore defined dihedrals, 1 = only weight 1,4 if in dihedral */
  int special_extra;         /**< Extra space for dynamically added bonds */
  int special_onefive;       /**< 0 = 1-5 neighbors not stored, 1 = stored */

  /** Constructor for Force class
   * \param lmp Pointer to the main LAMMPS instance */
  Force(class LAMMPS *);

  /** Destructor - cleans up all force field styles and maps */
  ~Force() override;

  /** Initialize all force field components before a run */
  void init();

  /** Setup all force field components at start of run */
  void setup();

  /** Create a pair style
   * \param style Pair style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_pair(const std::string &, int);

  /** Create a new pair style instance
   * \param style Pair style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new Pair instance */
  Pair *new_pair(const std::string &, int, int &);

  /** Find a pair style matching a pattern
   * \param style Style name pattern to match
   * \param exact 1 for exact match, 0 for partial match
   * \param nsub Sub-style index for hybrid (default 0)
   * \return Pointer to matching Pair, or nullptr */
  Pair *pair_match(const std::string &, int, int nsub = 0);

  /** Get the style name of a pair style
   * \param ptr Pointer to Pair instance
   * \return Style name string */
  char *pair_match_ptr(Pair *);

  /** Create a bond style
   * \param style Bond style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_bond(const std::string &, int);

  /** Create a new bond style instance
   * \param style Bond style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new Bond instance */
  Bond *new_bond(const std::string &, int, int &);

  /** Find a bond style matching a pattern
   * \param style Style name pattern to match
   * \return Pointer to matching Bond, or nullptr */
  Bond *bond_match(const std::string &);

  /** Create an angle style
   * \param style Angle style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_angle(const std::string &, int);

  /** Create a new angle style instance
   * \param style Angle style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new Angle instance */
  Angle *new_angle(const std::string &, int, int &);

  /** Find an angle style matching a pattern
   * \param style Style name pattern to match
   * \return Pointer to matching Angle, or nullptr */
  Angle *angle_match(const std::string &);

  /** Create a dihedral style
   * \param style Dihedral style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_dihedral(const std::string &, int);

  /** Create a new dihedral style instance
   * \param style Dihedral style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new Dihedral instance */
  Dihedral *new_dihedral(const std::string &, int, int &);

  /** Find a dihedral style matching a pattern
   * \param style Style name pattern to match
   * \return Pointer to matching Dihedral, or nullptr */
  Dihedral *dihedral_match(const std::string &);

  /** Create an improper style
   * \param style Improper style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_improper(const std::string &, int);

  /** Create a new improper style instance
   * \param style Improper style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new Improper instance */
  Improper *new_improper(const std::string &, int, int &);

  /** Find an improper style matching a pattern
   * \param style Style name pattern to match
   * \return Pointer to matching Improper, or nullptr */
  Improper *improper_match(const std::string &);

  /** Create a kspace style
   * \param style KSpace style name
   * \param trysuffix Flag to try accelerated suffix versions */
  void create_kspace(const std::string &, int);

  /** Create a new kspace style instance
   * \param style KSpace style name
   * \param trysuffix Flag to try accelerated suffix versions
   * \param sflag Output: 1 if suffix was used, 0 otherwise
   * \return Pointer to new KSpace instance */
  KSpace *new_kspace(const std::string &, int, int &);

  /** Find a kspace style matching a pattern
   * \param style Style name pattern to match
   * \param exact 1 for exact match, 0 for partial match
   * \return Pointer to matching KSpace, or nullptr */
  KSpace *kspace_match(const std::string &, int);

  /** Store a style name with optional suffix
   * \param style Style name
   * \param sflag Suffix flag
   * \return Allocated style name string */
  char *store_style(const std::string &, int);

  /** Set special bonds settings
   * \param narg Number of arguments
   * \param arg Argument array */
  void set_special(int, char **);

  /** Calculate memory usage
   * \return Memory usage in bytes */
  double memory_usage();

 private:
  void create_factories();
};

}    // namespace LAMMPS_NS

#endif
