// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */

#include "manifold_gaussian_bump.h"

#include "comm.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace user_manifold;

// This manifold comes with some baggage;
// it needs a taper function that is sufficiently smooth
// so we use a cubic Hermite interpolation and then for
// efficiency we construct a lookup table for that....

class cubic_hermite
{
public:
  // Represent x-polynomial as  a * t^3 + b * t^2 + c*t + d.
  double a, b, c, d;
  // And y-polynomial as s * t^3 + u * t^2 + v * t + w.
  double s, u, v, w;

  double x0, x1, y0, y1, yp0, yp1;

  LAMMPS_NS::Error *err;

  cubic_hermite( double x0, double x1, double y0, double y1,
                 double yp0, double yp1, LAMMPS_NS::Error *err ) :
    a(  2*x0 + 2 - 2*x1 ),
    b( -3*x0 - 3 + 3*x1 ),
    c( 1.0 ),
    d( x0 ),
    s(  2*y0 - 2*y1 +   yp0 + yp1 ),
    u( -3*y0 + 3*y1 - 2*yp0 - yp1  ),
    v(  yp0  ),
    w(  y0 ),
    x0(x0), x1(x1), y0(y0), y1(y1), yp0(yp0), yp1(yp1),
    err(err)
  {
    test();
  }


  void test()
  {
    if (fabs( x(0) - x0 ) > 1e-8 ) err->one(FLERR, "x0 wrong");
    if (fabs( x(1) - x1 ) > 1e-8 ) err->one(FLERR, "x1 wrong");
    if (fabs( y(0) - y0 ) > 1e-8 ) err->one(FLERR, "y0 wrong");
    if (fabs( y(1) - y1 ) > 1e-8 ) err->one(FLERR, "y1 wrong");
  }

  double get_t_from_x( double xx ) const
  {
    if (xx < x0 || xx > x1) {
      char msg[2048];
      sprintf(msg, "x ( %g ) out of bounds [%g, %g]", xx, x0, x1 );
      err->one(FLERR, msg);
    }

    // Newton iterate to get right t.
    double t  = xx - x0;
    double dx = x1 - x0;
    // Reasonable initial guess I hope:
    t /= dx;

    double tol = 1e-8;
    int maxit = 500;
    double ff  = x(t) - xx;
    double ffp = xp(t);
    // double l   = 1.0 / ( 1 + res*res );
    for (int i = 0; i < maxit; ++i) {
      t -= ff / ffp;
      ff  = x(t) - xx;
      ffp = xp(t);
      double res = ff;
      if (fabs( res ) < tol) {
        return t;
      }
    }
    err->warning(FLERR, "Convergence failed");
    return t;
  }

  double x( double t ) const
  {
    double t2 = t*t;
    double t3 = t2*t;
    return a*t3 + b*t2 + c*t + d;
  }

  double y_from_x( double x ) const
  {
    double t = get_t_from_x( x );
    return y(t);
  }

  double yp_from_x( double x ) const
  {
    double t = get_t_from_x( x );
    return yp(t);
  }

  double y( double t ) const
  {
    double t2 = t*t;
    double t3 = t2*t;
    return s*t3 + u*t2 + v*t + w;
  }

  void xy( double t, double &xx, double &yy ) const
  {
    xx = x(t);
    yy = y(t);
  }

  double xp( double t ) const
  {
    double t2 = t*t;
    return 3*a*t2 + 2*b*t + c;
  }

  double yp( double t ) const
  {
    double t2 = t*t;
    return 3*t2*s + 2*u*t + v;
  }

  double xpp( double t ) const
  {
    return 6*a*t + 2*b;
  }

};

// Manifold itself:
manifold_gaussian_bump::manifold_gaussian_bump(class LAMMPS* lmp,
                                               int /*narg*/, char **/*arg*/)
        : manifold(lmp), lut_z(nullptr), lut_zp(nullptr) {}


manifold_gaussian_bump::~manifold_gaussian_bump()
{
  if (lut_z ) delete lut_z;
  if (lut_zp) delete lut_zp;
}


// The constraint is z = f(r = sqrt( x^2 + y^2) )
// f(x) = A*exp( -x^2 / 2 l^2 )       if x < rc1
//      = Some interpolation          if rc1 <= rc2
//      = 0                           if x >= rc2
//
double manifold_gaussian_bump::g( const double *x )
{
  double xf[3];
  xf[0] = x[0];
  xf[1] = x[1];
  xf[2] = 0.0;

  double x2 = dot(xf,xf);
  if (x2 < rc12) {
    return x[2] - gaussian_bump_x2( x2 );
  } else if (x2 < rc22) {
    double rr = sqrt( x2 );
    double z_taper_func = lut_get_z( rr );
    return x[2] - z_taper_func;
  } else {
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

  if (x2 < rc12) {
    double factor = gaussian_bump_x2(x2);
    factor /= (ll*ll);
    nn[0] = factor * x[0];
    nn[1] = factor * x[1];
  } else if (x2 < rc22) {
    double rr = sqrt( x2 );
    double zp_taper_func = lut_get_zp( rr );

    double inv_r = 1.0 / rr;
    double der_part = zp_taper_func * inv_r;

    nn[0] = der_part * x[0];
    nn[1] = der_part * x[1];

  } else {
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
  if (x2 < rc12) {
    double gb = gaussian_bump_x2(x2);
    double factor = gb / (ll*ll);
    nn[0] = factor * x[0];
    nn[1] = factor * x[1];

    return x[2] - gb;
  } else if (x2 < rc22) {
    double z_taper_func, zp_taper_func;
    double rr = sqrt( x2 );
    lut_get_z_and_zp( rr, z_taper_func, zp_taper_func );

    double inv_r = 1.0 / rr;
    double der_part = zp_taper_func * inv_r;

    nn[0] = der_part * x[0];
    nn[1] = der_part * x[1];

    return x[2] - z_taper_func;
  } else {
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

  rc12 = rc1*rc1;
  rc22 = rc2*rc2;
  dr = rc2 - rc1;
  inv_dr = 1.0 / dr;

  f_at_rc  = gaussian_bump_x2 ( rc12 );
  fp_at_rc = gaussian_bump_der( rc1 );



  make_lut();

  // test_lut();
}


double manifold_gaussian_bump::gaussian_bump( double x ) const
{
  double x2 = x*x;
  return gaussian_bump_x2( x2 );
}

double manifold_gaussian_bump::gaussian_bump_x2( double x2 ) const
{
  return AA*exp( -x2 / ll2 );
}

double manifold_gaussian_bump::gaussian_bump_der( double x ) const
{
  double x2 = x*x;
  return gaussian_bump_x2( x2 )*( -x/(ll*ll) );
}


void manifold_gaussian_bump::make_lut()
{

  lut_x0 = rc1;
  lut_x1 = rc2;
  lut_Nbins = 1024; // Don't make it too big, it should fit in cache.
  lut_z  = new double[lut_Nbins+1];
  lut_zp = new double[lut_Nbins+1];

  lut_dx = (lut_x1 - lut_x0) / lut_Nbins;

  cubic_hermite pchip( lut_x0, lut_x1, f_at_rc, 0.0, fp_at_rc, 0.0, error );

  double xx = lut_x0;
  for (int i = 0; i <= lut_Nbins; ++i) {
    lut_z[i]  = pchip.y_from_x( xx );
    lut_zp[i] = pchip.yp_from_x( xx );
    xx += lut_dx;
  }
}


double manifold_gaussian_bump::lut_get_z ( double rr ) const
{
  double xs   = rr - lut_x0;
  double xss  = xs / lut_dx;
  int    bin  = static_cast<int>(xss);
  double frac = xss - bin;

  double zleft  = lut_z[bin];
  double zright = lut_z[bin+1];

  return zleft * ( 1 - frac ) + frac * zright;
}

double manifold_gaussian_bump::lut_get_zp( double rr ) const
{
  double xs   = rr - lut_x0;
  double xss  = xs / lut_dx;
  int    bin  = static_cast<int>(xss);
  double frac = xss - bin;

  double zleft  = lut_zp[bin];
  double zright = lut_zp[bin+1];

  return zleft * ( 1 - frac) + frac * zright;
}


void manifold_gaussian_bump::lut_get_z_and_zp( double rr, double &zz,
                                               double &zzp ) const
{
  double xs   = rr - lut_x0;
  double xss  = xs / lut_dx;
  int    bin  = static_cast<int>(xss);
  double frac = xss - bin;
  double fmin = 1 - frac;

  double zleft   = lut_z[bin];
  double zright  = lut_z[bin+1];
  double zpleft  = lut_zp[bin];
  double zpright = lut_zp[bin+1];

  zz  =  zleft * fmin +  zright * frac;
  zzp = zpleft * fmin + zpright * frac;
}


void manifold_gaussian_bump::test_lut()
{
  double x[3], nn[3];
  if (comm->me != 0) return;

  FILE *fp = fopen( "test_lut_gaussian.dat", "w" );
  double dx = 0.1;
  for (double xx = 0; xx < 20; xx += dx) {
    x[0] = xx;
    x[1] = 0.0;
    x[2] = 0.0;
    double gg = g( x );
    n( x, nn );
    double taper_z;
    if (xx <= rc1) {
            taper_z = gaussian_bump(xx);
    } else if (xx < rc2) {
            taper_z = lut_get_z( xx );
    } else {
            taper_z = 0.0;
    }
    fprintf( fp, "%g %g %g %g %g %g %g\n", xx, gaussian_bump(xx), taper_z,
             gg, nn[0], nn[1], nn[2] );
  }
  fclose(fp);
}
