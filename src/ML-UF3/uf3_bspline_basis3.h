// clang-format off
/* ----------------------------------------------------------------------
   lammps - large-scale atomic/molecular massively parallel simulator
   https://www.lammps.org/, sandia national laboratories
   lammps development team: developers@lammps.org

   copyright (2003) sandia corporation.  under the terms of contract
   de-ac04-94al85000 with sandia corporation, the u.s. government retains
   certain rights in this software.  this software is distributed under
   the gnu general public license.

   see the readme file in the top-level lammps directory.
------------------------------------------------------------------------- */

#include "pointers.h"

#include <vector>

#ifndef UF3_BSPLINE_BASIS3_H
#define UF3_BSPLINE_BASIS3_H

namespace LAMMPS_NS {

class uf3_bspline_basis3 {
 private:
  LAMMPS *lmp;
  //std::vector<double> constants;

 public:
  uf3_bspline_basis3(LAMMPS *ulmp, const double *knots, double coefficient);
  ~uf3_bspline_basis3();
  double constants[16] = {};
  double eval0(double, double, double);
  double eval1(double, double, double);
  double eval2(double, double, double);
  double eval3(double, double, double);

  double memory_usage();
};

}    // namespace LAMMPS_NS
#endif
