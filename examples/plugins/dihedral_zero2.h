/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Identical to dihedral harmonic, except if all k's are zero the
   force loop is skipped.

------------------------------------------------------------------------- */

#ifndef LMP_DIHEDRAL_ZERO2_H
#define LMP_DIHEDRAL_ZERO2_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralZero2 : public Dihedral {
 public:
  DihedralZero2(class LAMMPS *);
  virtual ~DihedralZero2();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  virtual void settings(int, char **);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  int coeffflag;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Incorrect args for dihedral coefficients

UNDOCUMENTED

*/
