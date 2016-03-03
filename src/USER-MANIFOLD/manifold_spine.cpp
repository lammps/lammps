#include "manifold_spine.h"

#include <math.h>

using namespace LAMMPS_NS;
using namespace user_manifold;

manifold_spine::manifold_spine( LAMMPS *lmp, int argc, char **argv ) : manifold(lmp)
{}


/*
 * Equation for spine is:
 * 
 * -x^2 - y^2 + (a^2 - (z/c)^2)*( 1 + A*sin(const(z) *z^2) )^4
 * const(z) = (z < 0) ? B2 : B
 * params[5] = { a, A, B, B2, c }
 */


double manifold_spine::g( const double *x )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double B2 = params[3];
  double c  = params[4];
  if (x[2]>0){
    return -1.0*(x[0]*x[0]+x[1]*x[1])+(a*a-x[2]*x[2]/(c*c))*(1.0+pow(A*sin(B *x[2]*x[2]),4));
  }else{
    return -1.0*(x[0]*x[0]+x[1]*x[1])+(a*a-x[2]*x[2])      *(1.0+pow(A*sin(B2*x[2]*x[2]),4));
  }
}

void manifold_spine::n( const double *x, double *nn )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double B2 = params[3];
  double c  = params[4];
  if (x[2]>0){
    nn[0] = -2.0*x[0];
    nn[1] = -2.0*x[1];
    nn[2] = 8.0*A*A*A*A*B*x[2]*(a*a-(x[2]*x[2]/(c*c)))*cos(B*x[2]*x[2]) *
	   pow(sin(B*x[2]*x[2]),3)-2*x[2]*(1+pow(A*sin(B*x[2]*x[2]),4))/(c*c);
  }else{
    nn[0] = -2.0*x[0];
    nn[1] = -2.0*x[1];
    nn[2] = 8.0*A*A*A*A*B2*x[2]*(a*a-x[2]*x[2])*cos(B2*x[2]*x[2]) *
           pow(sin(B2*x[2]*x[2]),3)-2*x[2]*(1+pow(A*sin(B2*x[2]*x[2]),4));
  }
}


