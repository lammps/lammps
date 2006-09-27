/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FORCE_H
#define FORCE_H

#include "lammps.h"

class Pair;
class Bond;
class Angle;
class Dihedral;
class Improper;
class KSpace;
class Temperature;
class Pressure;

class Force : public LAMMPS {
 public:
  double boltz;                      // Boltzmann constant (eng/degree-K)
  double mvv2e;                      // conversion of mv^2 to energy
  double ftm2v;                      // conversion of ft/m to velocity
  double nktv2p;                     // conversion of NkT/V to pressure
  double qqr2e;                      // conversion of q^2/r to energy
  double qe2f;                       // conversion of qE to force
  double dielectric;                 // dielectric constant
  double qqrd2e;                     // q^2/r to energy w/ dielectric constant

  int dimension;                        // 2 = 2d, 3 = 3d
  int newton,newton_pair,newton_bond;   // Newton's 3rd law settings

  Pair *pair;
  char *pair_style;

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
                             // index [0] is not used in these arrays
  double special_lj[4];      // 1-2, 1-3, 1-4 prefactors for LJ
  double special_coul[4];    // 1-2, 1-3, 1-4 prefactors for Coulombics

  int ntemp;                   // # of defined Temperatures
  int maxtemp;                 // max number list can hold
  Temperature **templist;      // list of defined Temperatures

  Pressure *pressure;          // single defined Pressure

  Force();
  ~Force();
  void init();

  void create_pair(char *);
  Pair *new_pair(char *);
  Pair *pair_match(char *);

  void create_bond(char *);
  Bond *new_bond(char *);
  Bond *bond_match(char *);

  void create_angle(char *);
  Angle *new_angle(char *);

  void create_dihedral(char *);
  Dihedral *new_dihedral(char *);

  void create_improper(char *);
  Improper *new_improper(char *);

  void create_kspace(int, char **);

  void set_special(int, char **);
  void add_temp(int, char **, int);
  Temperature *find_temp(char *);
  void modify_temp(int, char **);
  void bounds(char *, int, int &, int &);
  int memory_usage();
};

#endif
