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

#ifdef DIHEDRAL_CLASS

DihedralStyle(hybrid,DihedralHybrid)

#else

#ifndef LMP_DIHEDRAL_HYBRID_H
#define LMP_DIHEDRAL_HYBRID_H

#include <cstdio>
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHybrid : public Dihedral {
 public:
  int nstyles;                  // # of different dihedral styles
  Dihedral **styles;            // class list for each Dihedral style
  char **keywords;              // keyword for each dihedral style

  DihedralHybrid(class LAMMPS *);
  ~DihedralHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 private:
  int *map;                     // which style each dihedral type points to

  int *ndihedrallist;           // # of dihedrals in sub-style dihedrallists
  int *maxdihedral;             // max # of dihedrals sub-style lists can store
  int ***dihedrallist;          // dihedrallist for each sub-style

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dihedral style hybrid cannot use same dihedral style twice

Self-explanatory.

E: Dihedral style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Dihedral style hybrid cannot have none as an argument

Self-explanatory.

E: Dihedral coeff for hybrid has invalid style

Dihedral style hybrid uses another dihedral style as one of its
coefficients.  The dihedral style used in the dihedral_coeff command
or read from a restart file is not recognized.

*/
