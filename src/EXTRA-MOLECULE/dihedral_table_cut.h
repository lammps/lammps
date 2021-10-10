/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(table/cut,DihedralTableCut);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_TABLE_CUT_H
#define LMP_DIHEDRAL_TABLE_CUT_H

#include "dihedral_table.h"

namespace LAMMPS_NS {

class DihedralTableCut : public DihedralTable {
 public:
  DihedralTableCut(class LAMMPS *);
  virtual ~DihedralTableCut();
  virtual void compute(int, int);
  virtual void coeff(int, char **);

 protected:
  double *aat_k, *aat_theta0_1, *aat_theta0_2;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

*/
