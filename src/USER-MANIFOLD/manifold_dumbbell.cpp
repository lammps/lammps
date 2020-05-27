#include "manifold_dumbbell.h"

#include <cmath>

using namespace LAMMPS_NS;

using namespace user_manifold;

manifold_dumbbell::manifold_dumbbell( LAMMPS *lmp, int /*argc*/, char **/*argv*/ ) : manifold(lmp)
{}


/*
 * Equation for dumbbell is:
 *
 * -x^2 - y^2 + (a^2 - (z/c)^2)*( 1 + A*sin(const(z) *z^2) )^4
 * const(z) = (z < 0) ? B2 : B
 * params[4] = { a, A, B, c }
 */


double manifold_dumbbell::g( const double *x )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double c  = params[3];
  return -1.0*(x[0]*x[0]+x[1]*x[1])+(a*a-x[2]*x[2]/(c*c))*(1.0+pow(A*sin(B *x[2]*x[2]),4));
}

void manifold_dumbbell::n( const double *x, double *nn )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double c  = params[3];
  nn[0] = -2.0*x[0];
  nn[1] = -2.0*x[1];
  nn[2] = 8.0*A*A*A*A*B*x[2]*(a*a-(x[2]*x[2]/(c*c)))*cos(B*x[2]*x[2]) *
          pow(sin(B*x[2]*x[2]),3)-2*x[2]*(1+pow(A*sin(B*x[2]*x[2]),4))/(c*c);
}


