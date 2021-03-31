#include "manifold_cylinder_dent.h"
#include "math_const.h"

#include <cmath>


using namespace LAMMPS_NS;

using namespace user_manifold;

manifold_cylinder_dent::manifold_cylinder_dent( LAMMPS *lmp, int /*argc*/,
                                                char **/*argv*/ ) : manifold(lmp)
{}


double manifold_cylinder_dent::g( const double *x )
{
  double l = params[1], R = params[0], a = params[2];
  double r2 = x[1]*x[1] + x[0]*x[0];
  if (fabs(x[2]) < 0.5*l) {
    double k = MathConst::MY_2PI / l;
    double c = R - 0.5*a*( 1.0 + cos(k*x[2]) );
    return c*c - r2;
  } else {
    return R*R - r2;
  }
}


void manifold_cylinder_dent::n( const double *x, double *n )
{
  double l = params[1], R = params[0], a = params[2];
  if (fabs(x[2]) < 0.5*l) {
    double k = MathConst::MY_2PI / l;
    double c = R - 0.5*a*(1.0 + cos(k*x[2]));
    n[0] = -2*x[0];
    n[1] = -2*x[1];
    n[2] = c*a*k*sin(k*x[2]);
  } else {
    n[0] = -2*x[0];
    n[1] = -2*x[1];
    n[2] = 0.0;
  }
}
