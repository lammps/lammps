//This code was written by Philip Nicoletti
//http://www.codeguru.com/forum/archive/index.php/t-129990.html
//
//Modified by Jin Ma, Oklahoma State University for LAMMPS
//erfc() is defined in GNU libraries. This code is a simplified
//version for implementation with Visual C++.
//
//Warning: these functions are not fully tested.
//
#include "erfc.h"
#include "math.h"

double erf(double x)
{
    //
    // Computation of the error function erf(x).
    //
    return (1-erfc(x));
}

//
//
double erfc(double x)
{
    //
    // Computation of the complementary error function erfc(x).
    //
    // The algorithm is based on a Chebyshev fit as denoted in
    // Numerical Recipes 2nd ed. on p. 214 (W.H.Press et al.).
    //
    // The fractional error is always less than 1.2e-7.
    //
    //
    // The parameters of the Chebyshev fit
    //
    const double a1 = -1.26551223, a2 = 1.00002368,
    a3 = 0.37409196, a4 = 0.09678418,
    a5 = -0.18628806, a6 = 0.27886807,
    a7 = -1.13520398, a8 = 1.48851587,
    a9 = -0.82215223, a10 = 0.17087277;
    //
    double v = 1; // The return value
    double z = fabs(x);
    //
    if (z == 0) return v; // erfc(0)=1
    double t = 1/(1+0.5*z);
    v = t*exp((-z*z) +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+
				t*(a7+t*(a8+t*(a9+t*a10)))))))));
    if (x < 0) v = 2-v;	  // erfc(-x)=2-erfc(x)
    return v;
}