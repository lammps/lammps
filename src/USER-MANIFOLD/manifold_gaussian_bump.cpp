#include "manifold_gaussian_bump.h"

using namespace LAMMPS_NS;
using namespace user_manifold;

// The constraint is z = f(r = sqrt( x^2 + y^2) )
// f(x) = A*exp( -x^2 / 2 l^2 )       if x < rc1
//      = a + b*x + c*x**2 + d*x**3   if rc1 <= x < rc2
//      = 0                           if x >= rc2
//
double manifold_gaussian_bump::g( const double *x )
{
	double xf[3];
	xf[0] = x[0];
	xf[1] = x[1];
	xf[2] = 0.0;
	
	double x2 = dot(xf,xf);
	if( x2 < rc12 ){
		return x[2] - gaussian_bump_x2( x2 );
	}else if( x2 < rc22 ){
		double rr = sqrt(x2);
		double xi = rr - rc1;
		xi *= inv_dr;
		double xi2 = x2 * inv_dr*inv_dr;
		double xi3 = xi*xi2;
		return x[2] - ( aa + bb*xi + cc*xi2 + dd*xi3 );
		
	}else{
		return x[2];
	}
}

void   manifold_gaussian_bump::n( const double *x, double *nn )
{
	double xf[3];
	xf[0] = x[0];
	xf[1] = x[1];
	xf[2] = 0.0;
	nn[2] = 1.0;
	
	double x2 = dot(xf,xf);
	
	if( x2 < rc12 ){
		double factor = gaussian_bump_x2(x2);
		factor /= (ll*ll);
		nn[0] = factor * x[0];
		nn[1] = factor * x[1];
	}else if( x2 < rc22 ){
		double rr = sqrt(x2);
		double xi = rr - rc1;
		xi *= inv_dr;
		double xi2 = x2 * inv_dr*inv_dr;
		double der = bb + 2*cc*xi + 3*dd*xi2;
		
		nn[0] = -der * x[0] / rr;
		nn[1] = -der * x[1] / rr;
	}else{
		nn[0] = nn[1] = 0.0;
	}
}

double manifold_gaussian_bump::g_and_n( const double *x, double *nn )
{
	double xf[3];
	xf[0] = x[0];
	xf[1] = x[1];
	xf[2] = 0.0;
	nn[2] = 1.0;
	
	double x2 = dot(xf,xf);
	if( x2 < rc12 ){
		double gb = gaussian_bump_x2(x2);
		double factor = gb / (ll*ll);
		nn[0] = factor * x[0];
		nn[1] = factor * x[1];
		
		return x[2] - gb;
	}else if( x2 < rc22 ){
		
		double rr = sqrt(x2);
		double xi = rr - rc1;
		xi *= inv_dr;
		double xi2 = x2 * inv_dr*inv_dr;
		double xi3 = xi*xi2;

		double der = bb + 2*cc*xi + 3*dd*xi2;
		
		nn[0] = -der * x[0] / rr;
		nn[1] = -der * x[1] / rr;

		
		return x[2] - ( aa + bb*xi + cc*xi2 + dd*xi3 );
		
	}else{
		nn[0] = nn[1] = 0.0;
		return x[2];
	}
	
}


void manifold_gaussian_bump::post_param_init()
{
	// Read in the params:
	AA  = params[0];
	ll  = params[1];
	rc1 = params[2];
	rc2 = params[3];

	ll2 = 2.0*ll*ll;

	f_at_rc  = gaussian_bump_x2 ( rc12 );
	fp_at_rc = gaussian_bump_der( rc12 );

	rc12 = rc1*rc1;
	rc22 = rc2*rc2;
	dr = rc2 - rc1;
	inv_dr = 1.0 / dr;
}


double manifold_gaussian_bump::gaussian_bump( double x )
{
	double x2 = x*x;
	return gaussian_bump_x2( x2 );
}

double manifold_gaussian_bump::gaussian_bump_x2( double x2 )
{
	return AA*exp( -x2 / ll2 );
}

double manifold_gaussian_bump::gaussian_bump_der( double x )
{
	double x2 = x*x;
	return gaussian_bump_x2( x2 )*( -x/(ll*ll) );
}
