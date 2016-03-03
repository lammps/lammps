#include "sincos_lookup.h"

#include <cmath>
#include <iostream>

#define SINCOS_LOOKUP_NO_SAFETY


#if __cplusplus < 201103L
  #define my_pi2 1.57079632679489661923
  #define my_pi  3.14159265358979323846
  #define my_2pi 6.28318530717958647692

  #define nullptr NULL

#else
  static constexpr const double my_pi2 = 1.57079632679489661923;
  static constexpr const double my_pi  = 3.14159265358979323846;
  static constexpr const double my_2pi = 6.28318530717958647692;
#endif // __cplusplus


using namespace LAMMPS_NS;
using namespace user_manifold;


sincos_lookup::sincos_lookup( int N ) : cos_table(nullptr), size(N)
{
	std::cout << "Constructing table of size " << size
	          << ", this might take a while... ";
	make_cos_table();
	std::cout << "Done!\n";
}

sincos_lookup::~sincos_lookup()
{
	if( cos_table ) delete [] cos_table;
}

double sincos_lookup::sin( double x )
{
	return this->cos( my_pi2 - x );
}

double sincos_lookup::cos( double x )
{
	double x0 = x;

	// Wrap to first quadrant:
	x = fabs(x);
	
	const double rounded_n_over_2pi = static_cast<int>(x / my_2pi) * my_2pi;
	x -= rounded_n_over_2pi;
	if( x > my_pi ){
		// make use of cos( x ) = -cos(pi-x)
		x -= 2.0*( x - my_pi );
	}
	// Get index:
	int idx = static_cast<int>( x / dx_tab );
#ifdef SINCOS_LOOKUP_NO_SAFETY
	// Assume safe:
#else	
	if( idx > size || idx < 0 ){
		std::cerr << "Illegal table lookup for index " << idx << ", x="
		          << x << ", x0=" << x0 << "!\n";
	}
#endif // SINCOS_LOOKUP_NO_SAFETY
	
	// static_cast always truncates to left, so lookup with "delta to right":
	double delta = x - idx*dx_tab;
	double cos_idx = cos_table[idx];
	double diff_cos = cos_table[idx+1] - cos_idx;
	return cos_idx + delta*diff_cos*inv_dx_tab;
}


void sincos_lookup::sincos(double x, double *s, double *c)
{
	*s = this->sin(x);
	*c = this->cos(x);
}

void sincos_lookup::test_lookup( )
{
	double dx = 1e-5;
	double max_err_c = 0, max_err_s = 0;
  
	for( double x = 0; x < 2*my_2pi; x += dx ){
		double c = this->cos(x);
		double s = this->sin(x);
		// Compare with math.h:
		double c_err = c - cos(x);
		double s_err = s - sin(x);

		if( fabs(max_err_c) < fabs(c_err) ) max_err_c = c_err;
		if( fabs(max_err_s) < fabs(s_err) ) max_err_s = s_err;
	}
	std::cout << "Largest lookup table errors (sin/cos): ( " << max_err_s
	          << ", " << max_err_c << " ).\n";
}


void sincos_lookup::make_cos_table()
{
	cos_table = new double [size+1];
	if( !cos_table ){
		// error.
	}
	dx_tab = my_pi/static_cast<double>(size-1);
	inv_dx_tab = 1.0 / dx_tab;

	for( int i = 0; i <= size; ++i ){
		double f = i*dx_tab;
		cos_table[i] = std::cos(f);
	}	
}
