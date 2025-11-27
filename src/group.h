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

#ifndef LMP_GROUP_H
#define LMP_GROUP_H

#include "pointers.h"

namespace LAMMPS_NS {
class Region;

/** \class LAMMPS_NS::Group
 * \brief Atom group management for LAMMPS
 *
 * The Group class manages named groups of atoms in LAMMPS. Groups provide a
 * way to select subsets of atoms for various operations such as applying
 * fixes, computing properties, or output. Each atom can belong to multiple
 * groups simultaneously.
 *
 * Key features include:
 *
 * - Definition of groups by atom type, region, molecule ID, or other criteria
 * - Dynamic groups that update membership during simulation
 * - Efficient bitmask-based membership testing
 * - Calculation of group properties (mass, charge, COM, velocity, etc.)
 * - Up to 32 groups (limited by 32-bit bitmask)
 *
 * The special group "all" (index 0) always contains all atoms. Group
 * membership is stored as a bitmask in each atom's mask property, allowing
 * efficient O(1) membership testing.
 *
 * \sa LAMMPS_NS::Pointers */
class Group : protected Pointers {
 public:
  int ngroup;          /**< Number of defined groups */
  char **names;        /**< Array of group names [ngroup] */
  int *bitmask;        /**< One-bit mask for each group [ngroup] */
  int *inversemask;    /**< Inverse mask (all bits except group's bit) [ngroup] */
  int *dynamic;        /**< 1 if group is dynamic, 0 if static [ngroup] */

  /** Constructor for Group class
   * \param lmp Pointer to the main LAMMPS instance */
  Group(class LAMMPS *);

  /** Destructor - cleans up group names and arrays */
  ~Group() override;

  /** Assign atoms to a group based on criteria
   * \param narg Number of arguments
   * \param arg Argument array specifying group and criteria */
  void assign(int, char **);

  /** Assign atoms to a group from a single command string
   * \param str Command string for group assignment */
  void assign(const std::string &);

  /** Create a group from flagged atoms
   * \param name Name of the group to create
   * \param flag Array of flags (1 = include, 0 = exclude) for each atom */
  void create(const std::string &, int *);

  /** Find a group by name
   * \param name Group name to find
   * \return Group index (0 to ngroup-1), or -1 if not found */
  int find(const std::string &);

  /** Find a group by name, creating it if it doesn't exist
   * \param name Group name to find or create
   * \return Group index */
  int find_or_create(const char *);

  /** Get group bitmask by ID with error checking
   * \param id Group name
   * \param which Flag indicating context
   * \param file Source file for error message
   * \param line Error context string
   * \return Group bitmask */
  int get_bitmask_by_id(const std::string &, int, const std::string &, const std::string &);

  /** Write group info to restart file
   * \param fp File pointer for restart file */
  void write_restart(FILE *);

  /** Read group info from restart file
   * \param fp File pointer for restart file */
  void read_restart(FILE *);

  /** Count all atoms (fast path for group "all")
   * \return Total number of atoms across all processors */
  bigint count_all();

  /** Count atoms in a group
   * \param igroup Group index
   * \return Number of atoms in the group */
  bigint count(int);

  /** Count atoms in a group within a region
   * \param igroup Group index
   * \param region Pointer to Region
   * \return Number of atoms in group and region */
  bigint count(int, Region *);

  /** Calculate total mass of atoms in a group
   * \param igroup Group index
   * \return Total mass */
  double mass(int);

  /** Calculate total mass of atoms in a group within a region
   * \param igroup Group index
   * \param region Pointer to Region
   * \return Total mass in group and region */
  double mass(int, Region *);

  /** Calculate total charge of atoms in a group
   * \param igroup Group index
   * \return Total charge */
  double charge(int);

  /** Calculate total charge of atoms in a group within a region
   * \param igroup Group index
   * \param region Pointer to Region
   * \return Total charge in group and region */
  double charge(int, Region *);

  /** Calculate bounding box of atoms in a group
   * \param igroup Group index
   * \param minmax Output array [6]: xlo,xhi,ylo,yhi,zlo,zhi */
  void bounds(int, double *);

  /** Calculate bounding box of atoms in a group within a region
   * \param igroup Group index
   * \param minmax Output array [6]
   * \param region Pointer to Region */
  void bounds(int, double *, Region *);

  /** Calculate center of mass of a group
   * \param igroup Group index
   * \param masstotal Total mass of group
   * \param cm Output: center of mass coordinates [3] */
  void xcm(int, double, double *);

  /** Calculate center of mass within a region
   * \param igroup Group index
   * \param masstotal Total mass
   * \param cm Output: center of mass [3]
   * \param region Pointer to Region */
  void xcm(int, double, double *, Region *);

  /** Calculate center of mass velocity of a group
   * \param igroup Group index
   * \param masstotal Total mass
   * \param vcm Output: COM velocity [3] */
  void vcm(int, double, double *);

  /** Calculate center of mass velocity within a region
   * \param igroup Group index
   * \param masstotal Total mass
   * \param vcm Output: COM velocity [3]
   * \param region Pointer to Region */
  void vcm(int, double, double *, Region *);

  /** Calculate total force on a group
   * \param igroup Group index
   * \param fcm Output: total force [3] */
  void fcm(int, double *);

  /** Calculate total force on a group within a region
   * \param igroup Group index
   * \param fcm Output: total force [3]
   * \param region Pointer to Region */
  void fcm(int, double *, Region *);

  /** Calculate kinetic energy of a group
   * \param igroup Group index
   * \return Kinetic energy */
  double ke(int);

  /** Calculate kinetic energy of a group within a region
   * \param igroup Group index
   * \param region Pointer to Region
   * \return Kinetic energy */
  double ke(int, Region *);

  /** Calculate radius of gyration of a group
   * \param igroup Group index
   * \param masstotal Total mass
   * \param cm Center of mass [3]
   * \return Radius of gyration */
  double gyration(int, double, double *);

  /** Calculate radius of gyration within a region
   * \param igroup Group index
   * \param masstotal Total mass
   * \param cm Center of mass [3]
   * \param region Pointer to Region
   * \return Radius of gyration */
  double gyration(int, double, double *, Region *);

  /** Calculate angular momentum of a group
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param lmom Output: angular momentum [3] */
  void angmom(int, double *, double *);

  /** Calculate angular momentum within a region
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param lmom Output: angular momentum [3]
   * \param region Pointer to Region */
  void angmom(int, double *, double *, Region *);

  /** Calculate torque on a group
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param tq Output: torque [3] */
  void torque(int, double *, double *);

  /** Calculate torque on a group within a region
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param tq Output: torque [3]
   * \param region Pointer to Region */
  void torque(int, double *, double *, Region *);

  /** Calculate moment of inertia tensor of a group
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param inertia Output: inertia tensor [3][3] */
  void inertia(int, double *, double[3][3]);

  /** Calculate moment of inertia tensor within a region
   * \param igroup Group index
   * \param cm Center of mass [3]
   * \param inertia Output: inertia tensor [3][3]
   * \param region Pointer to Region */
  void inertia(int, double *, double[3][3], Region *);

  /** Calculate angular velocity from angular momentum and inertia
   * \param angmom Angular momentum [3]
   * \param inertia Inertia tensor [3][3]
   * \param omega Output: angular velocity [3] */
  void omega(double *, double[3][3], double *);

 private:
  int me;    /**< MPI rank of this processor */

  int find_unused();
  void add_molecules(int, int);

  // callback functions for ring communication

  static void molring(int, char *, void *);
  int molbit;    /**< Bitmask for molecule-based group operations */
};

}    // namespace LAMMPS_NS

#endif
