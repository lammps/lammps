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
DihedralStyle(nharmonic,DihedralNHarmonic);
// clang-format on
#else

#ifndef DIHEDRAL_NHARMONIC_H
#define DIHEDRAL_NHARMONIC_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralNHarmonic : public Dihedral {
 public:
  DihedralNHarmonic(class LAMMPS *);
  ~DihedralNHarmonic() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  void born_matrix(int /*dtype*/, int, int, int, int, double &, double &) override;

 protected:
  int *nterms;
  double **a;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
