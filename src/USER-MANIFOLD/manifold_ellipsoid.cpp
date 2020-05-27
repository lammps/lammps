#include "manifold_ellipsoid.h"

using namespace LAMMPS_NS;

using namespace user_manifold;

manifold_ellipsoid::manifold_ellipsoid( LAMMPS *lmp, int /*narg*/, char **/*argv*/ ) : manifold(lmp)
{}


double manifold_ellipsoid::g( const double *x )
{
  const double ai = 1.0 / params[0];
  const double bi = 1.0 / params[1];
  const double ci = 1.0 / params[2];
  return x[0]*x[0]*ai*ai + x[1]*x[1]*bi*bi + x[2]*x[2]*ci*ci - 1.0;
}

void manifold_ellipsoid::n( const double *x, double * n )
{
  const double ai = 1.0 / params[0];
  const double bi = 1.0 / params[1];
  const double ci = 1.0 / params[2];
  n[0] = 2*x[0]*ai*ai;
  n[1] = 2*x[1]*bi*bi;
  n[2] = 2*x[2]*ci*ci;
}

