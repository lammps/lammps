# ifndef CERF_H
# define CERF_H

# include <complex>
# include <cmath>
/*

Copyright (C) 1998, 1999 John Smith

This file is part of Octave.
Or it might be one day....

Octave is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

Octave is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

// Put together by John Smith john at arrows dot demon dot co dot uk, 
// using ideas by others.
//
// Calculate erf(z) for complex z.
// Three methods are implemented; which one is used depends on z.
//
// The code includes some hard coded constants that are intended to
// give about 14 decimal places of accuracy. This is appropriate for
// 64-bit floating point numbers. 
//
// Oct 1999: Fixed a typo that in
//     const Complex cerf_continued_fraction( const Complex z )
// that caused erroneous answers for erf(z) where real(z) negative
//


//
// Abramowitz and Stegun: (eqn: 7.1.14) gives this continued
// fraction for erfc(z)
//
// erfc(z) = sqrt(pi).exp(-z^2).  1   1/2   1   3/2   2   5/2  
//                               ---  ---  ---  ---  ---  --- ...
//                               z +  z +  z +  z +  z +  z +
//
// This is evaluated using Lentz's method, as described in the narative
// of Numerical Recipes in C.
//
// The continued fraction is true providing real(z)>0. In practice we 
// like real(z) to be significantly greater than 0, say greater than 0.5.
//
template< class Complex>
const Complex cerfc_continued_fraction( const Complex z )
{
  double tiny = 1e-20 ;     // a small number, large enough to calculate 1/tiny
  double eps = 1e-15 ;      // large enough so that 1.0+eps > 1.0, when using
                            // the floating point arithmetic
  //
  // first calculate z+ 1/2   1 
  //                    ---  --- ...
  //                    z +  z + 
  Complex f(z) ;
  Complex C(f) ;
  Complex D(0.0) ;
  Complex delta ;
  double a ;

  a = 0.0 ;
  do
    {
      a = a + 0.5 ;
      D = z + a*D ;
      C = z + a/C ;

      if (D.real() == 0.0 && D.imag() ==  0.0)
	D = tiny ;

      D = 1.0 / D ;

      delta = (C * D) ;

      f = f * delta ;

    } while (abs(1.0-delta) > eps ) ;

  //
  // Do the first term of the continued fraction
  //
  f = 1.0 / f ;

  //
  // and do the final scaling
  //
  f = f * exp(-z*z)/ sqrt(M_PI) ;

  return f ;
}

template< class Complex>
const Complex cerf_continued_fraction( const Complex z )
{
  //  warning("cerf_continued_fraction:");
  if (z.real() > 0)
    return 1.0 - cerfc_continued_fraction( z ) ;
  else
    return -1.0 + cerfc_continued_fraction( -z ) ;
}

//
// Abramawitz and Stegun, Eqn. 7.1.5 gives a series for erf(z)
// good for all z, but converges faster for smallish abs(z), say abs(z)<2.
//
template< class Complex>
const Complex cerf_series( const Complex z )
{
  double tiny = 1e-20 ;       // a small number compared with 1.
  // warning("cerf_series:");
  Complex sum(0.0) ;
  Complex term(z) ;
  Complex z2(z*z) ;

  for (int n=0; n<3 || abs(term) > abs(sum)*tiny; n++)
    {
      sum = sum + term / (2*n+1) ;
      term = -term * z2 / (n+1) ;
    }

  return sum * 2.0 / sqrt(M_PI) ;
}
  
//
// Numerical Recipes quotes a formula due to Rybicki for evaluating 
// Dawson's Integral:
//
// exp(-x^2) integral  exp(t^2).dt = 1/sqrt(pi) lim   sum  exp(-(z-n.h)^2) / n
//            0 to x                            h->0 n odd
//
// This can be adapted to erf(z).
//
template< class Complex>
const Complex cerf_rybicki( const Complex z )
{
  // warning("cerf_rybicki:");
  double h = 0.2 ;        // numerical experiment suggests this is small enough

  //
  // choose an even n0, and then shift z->z-n0.h and n->n-h. 
  // n0 is chosen so that real((z-n0.h)^2) is as small as possible. 
  // 
  int n0 = 2*(int) (floor( z.imag()/(2*h) + 0.5 )) ;

  Complex z0( 0.0, n0*h ) ;
  Complex zp(z-z0) ;
  Complex sum(0.0,0.0) ;
  //
  // limits of sum chosen so that the end sums of the sum are
  // fairly small. In this case exp(-(35.h)^2)=5e-22 
  //
  //
  for (int np=-35; np<=35; np+=2)
    {
      Complex t( zp.real(), zp.imag()-np*h) ;
      Complex b( exp(t*t) / (np+n0) ) ;
      sum += b ; 
    }

  sum = sum * 2 * exp(-z*z) / M_PI ;

  return Complex(-sum.imag(), sum.real()) ;
}

template< class Complex>
const Complex cerf( const Complex z )
{
  //
  // Use the method appropriate to size of z - 
  // there probably ought to be an extra option for NaN z, or infinite z
  //
  //
  if (abs(z) < 2.0)
    return cerf_series( z ) ;
  else if (abs(z.real()) < 0.5)
    return cerf_rybicki( z ) ;
  else
    return cerf_continued_fraction( z ) ;
}

//
// Footnote:
// 
// Using the definitions from Abramowitz and Stegun (7.3.1, 7.3.2)
// The fresnel intgerals defined as:
//
//         / t=x
// C(x) = |      cos(pi/2 t^2) dt
//        /
//         t=0 
//
// and
//         / t=x
// S(x) = |      sin(pi/2 t^2) dt
//        /
//         t=0 
//
// These can be derived from erf(x) using 7.3.22
//
// C(z) +iS(z) = (1+i)  erf( sqrt(pi)/2 (1-i) z )
//               -----
//                 2
//
// --------------------------------------------------------------------------
// Some test examples - 
// comparative data taken from Abramowitz and Stegun table 7.9. 
// Table 7.9 tabulates w(z), where w(z) = exp(-z*z) erfc(iz)
// I have copied twelve values of w(z) from the table, and separately
// calculated them using this code. The results are identical.
//
//  x    y   Abramowitz & Stegun |             Octave Calculations
//                  w(x+iy)      |          w(x+iy)          cerf ( i.(x+iy))
// 0.2  0.2   0.783538+0.157403i |   0.783538 +0.157403    0.23154672 -0.219516
// 0.2  0.7   0.515991+0.077275i |   0.515991 +0.077275    0.69741968 -0.138277
// 0.2  1.7   0.289309+0.027154i |   0.289309 +0.027154    0.98797507 -0.011744
// 0.2  2.7   0.196050+0.013002i |   0.196050 +0.013002    0.99994252 -0.000127
// 1.2  0.2   0.270928+0.469488i |   0.270928 +0.469488    0.90465623 -2.196064
// 1.2  0.7   0.280740+0.291851i |   0.280740 +0.291851    1.82926135 -0.639343
// 1.2  1.7   0.222436+0.129684i |   0.222436 +0.129684    1.00630308 +0.060067
// 1.2  2.7   0.170538+0.068617i |   0.170538 +0.068617    0.99955699 -0.000290
// 2.2  0.2   0.041927+0.287771i |   0.041927 +0.287771   24.70460755-26.205981
// 2.2  0.7   0.099943+0.242947i |   0.099943 +0.242947    9.88734713+18.310797
// 2.2  1.7   0.135021+0.153161i |   0.135021 +0.153161    1.65541359 -1.276707
// 2.2  2.7   0.127900+0.096330i |   0.127900 +0.096330    0.98619434 +0.000564

# endif
