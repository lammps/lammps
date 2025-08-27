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

#ifndef UF3_BSPLINE_BASIS2_H
#define UF3_BSPLINE_BASIS2_H

namespace LAMMPS_NS {

class uf3_bspline_basis2 {
 private:
  class LAMMPS *lmp;

 public:
  uf3_bspline_basis2(LAMMPS *ulmp, const double *knots, double coefficient);
  ~uf3_bspline_basis2() = default;
  double constants[9];

  // Evaluate outer-left part of spline
  double eval0(double rsq, double r) const
  {
    return rsq * constants[2] + r * constants[1] + constants[0];
  }

  // Evaluate center-left part of spline
  double eval1(double rsq, double r) const
  {
    return rsq * constants[5] + r * constants[4] + constants[3];
  }

  // Evaluate center-right part of spline
  double eval2(double rsq, double r) const
  {
    return rsq * constants[8] + r * constants[7] + constants[6];
  }

  double memory_usage();
};

}    // namespace LAMMPS_NS
#endif
