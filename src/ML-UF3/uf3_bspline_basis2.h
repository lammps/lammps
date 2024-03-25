//De Boor's algorithm @
//https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/de-Boor.html
//For values outside the domain,
//extrapoltaes the left(right) hand side piece of the curve
//Only works for bspline degree upto 3 becuase of definiation of P
//
#include "pointers.h"

#include <vector>

#ifndef UF3_BSPLINE_BASIS2_H
#define UF3_BSPLINE_BASIS2_H

namespace LAMMPS_NS {

class uf3_bspline_basis2 {
 private:
  LAMMPS *lmp;
  std::vector<double> constants;

 public:
  uf3_bspline_basis2(LAMMPS *ulmp, const double *knots, double coefficient);
  ~uf3_bspline_basis2();
  double eval0(double, double);
  double eval1(double, double);
  double eval2(double, double);

  double memory_usage();
};

}    // namespace LAMMPS_NS
#endif
