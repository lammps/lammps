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
DihedralStyle(multi/harmonic,DihedralMultiHarmonic);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_MULTI_HARMONIC_H
#define LMP_DIHEDRAL_MULTI_HARMONIC_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralMultiHarmonic : public Dihedral {
 public:
  DihedralMultiHarmonic(class LAMMPS *);
  ~DihedralMultiHarmonic() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *a1, *a2, *a3, *a4, *a5;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
