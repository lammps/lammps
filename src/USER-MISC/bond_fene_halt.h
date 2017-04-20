/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Identical to bond fene except it takes an additional argument which is
   the maximal distance of fene bonds. Simulation will be terminated if
   a bond larger than the max distance is encountered.

   E.G.
     bond_style fenehalt
     bond_coeff <bondtype> k R_0  \epsilon  \sigma  r_{halt}

------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(fenehalt,BondFENEHALT)

#else

#ifndef LMP_BOND_FENEHALT_H
#define LMP_BOND_FENEHALT_H

#include <stdio.h>
#include "bond.h"

namespace LAMMPS_NS {

class BondFENEHALT : public Bond {
 public:
  BondFENEHALT(class LAMMPS *);
  virtual ~BondFENEHALT();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  double TWO_1_3;
  double *k,*r0,*epsilon,*sigma;
  double *rhalt2;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: FENE bond too long: %ld %d %d %g

A FENE bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

E: FENE bond longer than specified halt distance

A FENE bond extended beyond the specific halt distance. Simulation terminated.

E: Bad FENE bond

Two atoms in a FENE bond have become so far apart that the bond cannot
be computed.

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

W: Use special bonds = 0,1,1 with bond style fene

Most FENE models need this setting for the special_bonds command.

W: FENE bond too long: %ld %g

A FENE bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

*/
