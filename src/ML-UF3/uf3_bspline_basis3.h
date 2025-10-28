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

#ifndef UF3_BSPLINE_BASIS3_H
#define UF3_BSPLINE_BASIS3_H

namespace LAMMPS_NS {

class uf3_bspline_basis3 {
 private:
  class LAMMPS *lmp;

 public:
  uf3_bspline_basis3(LAMMPS *ulmp, const double *knots, double coefficient);
  ~uf3_bspline_basis3() = default;

  double constants[16];

  // Evaluate outer-left part of spline
  double eval0(double rth, double rsq, double r) const
  {
    return rth * constants[3] + rsq * constants[2] + r * constants[1] + constants[0];
  }

  // Evaluate center-left part of spline
  double eval1(double rth, double rsq, double r) const
  {
    return rth * constants[7] + rsq * constants[6] + r * constants[5] + constants[4];
  }

  // Evaluate center-right part of spline
  double eval2(double rth, double rsq, double r) const
  {
    return rth * constants[11] + rsq * constants[10] + r * constants[9] + constants[8];
  }

  // Evaluate outer-right part of spline
  double eval3(double rth, double rsq, double r) const
  {
    return rth * constants[15] + rsq * constants[14] + r * constants[13] + constants[12];
  }

  double memory_usage();
};

}    // namespace LAMMPS_NS
#endif
