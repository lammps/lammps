/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/long/tip4p/long,PairLJLongTIP4PLong)

#else

#ifndef LMP_PAIR_LJ_LONG_TIP4P_LONG_H
#define LMP_PAIR_LJ_LONG_TIP4P_LONG_H

#include "pair_lj_long_coul_long.h"

namespace LAMMPS_NS {

class PairLJLongTIP4PLong : public PairLJLongCoulLong {
 public:
  PairLJLongTIP4PLong(class LAMMPS *);
  ~PairLJLongTIP4PLong();
  virtual void compute(int, int);
  void settings(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart_settings(FILE *fp);
  void read_restart_settings(FILE *fp);
  void *extract(const char *, int &);
  double memory_usage();


 protected:
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  int typeA,typeB;             // angle and bond types of TIP4P water
  double alpha;                // geometric constraint parameter for TIP4P

  int nmax;                    // info on off-oxygen charge sites
  int **hneigh;                // 0,1 = indices of 2 H associated with O
                               // 2 = 0 if site loc not yet computed, 1 if yes
  double **newsite;            // locations of charge sites

  void compute_newsite(double *, double *, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: TIP4P hydrogen is missing

The TIP4P pairwise computation failed to find the correct H atom
within a water molecule.

E: TIP4P hydrogen has incorrect atom type

The TIP4P pairwise computation found an H atom whose type does not
agree with the specified H type.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Mixing forced for lj coefficients

UNDOCUMENTED

W: Using largest cut-off for lj/coul long long

UNDOCUMENTED

E: Coulombic cut not supported in pair_style lj/coul

UNDOCUMENTED

E: Pair style lj/coul/tip4p requires atom IDs

UNDOCUMENTED

E: Pair style lj/coul/tip4p requires newton pair on

UNDOCUMENTED

E: Pair style lj/coul/tip4p requires atom attribute q

UNDOCUMENTED

E: Must use a bond style with TIP4P potential

TIP4P potentials assume bond lengths in water are constrained
by a fix shake command.

E: Must use an angle style with TIP4P potential

TIP4P potentials assume angles in water are constrained by a fix shake
command.

E: Water H epsilon must be 0.0 for pair style lj/coul/tip4p

UNDOCUMENTED

U: Pair style lj/cut/coul/long/tip4p requires atom IDs

There are no atom IDs defined in the system and the TIP4P potential
requires them to find O,H atoms with a water molecule.

U: Pair style lj/cut/coul/long/tip4p requires newton pair on

This is because the computation of constraint forces within a water
molecule adds forces to atoms owned by other processors.

U: Pair style lj/cut/coul/long/tip4p requires atom attribute q

The atom style defined does not have these attributes.

U: Pair style is incompatible with KSpace style

If a pair style with a long-range Coulombic component is selected,
then a kspace style must also be used.

*/
