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
DihedralStyle(opls,DihedralOPLS);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_OPLS_H
#define LMP_DIHEDRAL_OPLS_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralOPLS : public Dihedral {
 public:
  DihedralOPLS(class LAMMPS *);
  ~DihedralOPLS() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k1, *k2, *k3, *k4;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
