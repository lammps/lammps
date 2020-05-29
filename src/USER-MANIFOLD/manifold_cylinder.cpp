#include "manifold_cylinder.h"

using namespace LAMMPS_NS;

using namespace user_manifold;


manifold_cylinder::manifold_cylinder( LAMMPS *lmp, int /*argc*/,
                                      char **/*argv*/ ) : manifold(lmp)
{}


double manifold_cylinder::g( const double *x )
{
  double R = params[0];
  double r2 = x[0]*x[0] + x[1]*x[1];
  return R*R - r2;
}

void manifold_cylinder::n( const double *x, double *n )
{
  n[0] = -2*x[0];
  n[1] = -2*x[1];
  n[2] = 0.0;
}


