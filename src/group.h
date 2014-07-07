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

#ifndef LMP_GROUP_H
#define LMP_GROUP_H

#include "stdio.h"
#include "pointers.h"
#include <map>

namespace LAMMPS_NS {

class Group : protected Pointers {
 public:
  int ngroup;                  // # of defined groups
  char **names;                // name of each group
  int *bitmask;                // one-bit mask for each group
  int *inversemask;            // inverse mask for each group
  int *dynamic;                // 1 if dynamic, 0 if not

  Group(class LAMMPS *);
  ~Group();
  void assign(int, char **);         // assign atoms to a group
  void create(char *, int *);        // add flagged atoms to a group
  int find(const char *);            // lookup name in list of groups
  void write_restart(FILE *);
  void read_restart(FILE *);

  bigint count(int);                       // count atoms in group
  bigint count(int,int);                   // count atoms in group & region
  double mass(int);                        // total mass of atoms in group
  double mass(int,int);
  double charge(int);                      // total charge of atoms in group
  double charge(int,int);
  void bounds(int, double *);              // bounds of atoms in group
  void bounds(int, double *, int);
  void xcm(int, double, double *);         // center-of-mass coords of group
  void xcm(int, double, double *, int);
  void vcm(int, double, double *);         // center-of-mass velocity of group
  void vcm(int, double, double *, int);
  void fcm(int, double *);                 // total force on group
  void fcm(int, double *, int);
  double ke(int);                          // kinetic energy of group
  double ke(int, int);
  double gyration(int, double, double *);  // radius-of-gyration of group
  double gyration(int, double, double *, int);
  void angmom(int, double *, double *);    // angular momentum of group
  void angmom(int, double *, double *, int);
  void torque(int, double *, double *);    // torque on group
  void torque(int, double *, double *, int);
  void inertia(int, double *, double [3][3]);     // inertia tensor
  void inertia(int, double *, double [3][3], int);
  void omega(double *, double [3][3], double *);  // angular velocity

 private:
  int me;
  std::map<tagint,int> *hash;

  int find_unused();
  void add_molecules(int, int);

  // static variable for ring communication callback to access class data
  // callback functions for ring communication

  static Group *cptr;
  static void molring(int, char *);
  int molbit;
};

}

#endif

/* ERROR/WARNING messages:

E: Group command before simulation box is defined

The group command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find group delete group ID

Self-explanatory.

E: Cannot delete group all

Self-explanatory.

E: Cannot delete group currently used by a fix

Self-explanatory.

E: Cannot delete group currently used by a compute

Self-explanatory.

E: Cannot delete group currently used by a dump

Self-explanatory.

E: Cannot delete group currently used by atom_modify first

Self-explanatory.

E: Too many groups

The maximum number of atom groups (including the "all" group) is
given by MAX_GROUP in group.cpp and is 32.

E: Group region ID does not exist

A region ID used in the group command does not exist.

E: Variable name for group does not exist

Self-explanatory.

E: Variable for group is invalid style

Only atom-style variables can be used.

E: Group ID does not exist

A group ID used in the group command does not exist.

*/
