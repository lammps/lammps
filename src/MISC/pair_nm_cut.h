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

PairStyle(nm/cut,PairNMCut)

#else

#ifndef LMP_PAIR_NM_CUT_H
#define LMP_PAIR_NM_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNMCut : public Pair {
 public:
  PairNMCut(class LAMMPS *);
  virtual ~PairNMCut();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut;
  double **e0,**r0,**nn, **mm;
  double **nm,**e0nm,**r0n,**r0m,**offset;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_atoms command before simulation box is defined

The create_atoms command cannot be used before a read_data,
read_restart, or create_box command.

E: Cannot create_atoms after reading restart file with per-atom info

The per-atom info was stored to be used when by a fix that you
may re-define.  If you add atoms before re-defining the fix, then
there will not be a correct amount of per-atom info.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Create_atoms region ID does not exist

A region ID used in the create_atoms command does not exist.

E: Invalid basis setting in create_atoms command

UNDOCUMENTED

E: Molecule template ID for create_atoms does not exist

UNDOCUMENTED

W: Molecule template for create_atoms has multiple molecules

UNDOCUMENTED

E: Invalid atom type in create_atoms command

The create_box command specified the range of valid atom types.
An invalid type is being requested.

E: Create_atoms molecule must have coordinates

UNDOCUMENTED

E: Create_atoms molecule must have atom types

UNDOCUMENTED

E: Invalid atom type in create_atoms mol command

UNDOCUMENTED

E: Create_atoms molecule has atom IDs, but system does not

UNDOCUMENTED

E: Cannot create atoms with undefined lattice

Must use the lattice command before using the create_atoms
command.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: No overlap of box and region for create_atoms

Self-explanatory.

*/
