#include "PolynomialSolver.h"
#include <limits>
#include <cmath>
#include <iostream>
#include "ATC_Error.h"

namespace ATC {
  // Utility functions used by solvers, but not globally accessible.
  static const double PI_OVER_3 = acos(-1.0)*(1.0/3.0);
  static bool is_zero(double x)
  {
    static double GT_ZERO = 1.0e2*std::numeric_limits<double>::epsilon();
    static double LT_ZERO = -GT_ZERO;
    return x>LT_ZERO && x<GT_ZERO;
  }
  static double sign(double x)
  {
    static double s[] = {-1.0,1.0};
    return s[x>0];
  }

  // Linear solver
  int solve_linear(double c[2], double x0[1])
  {
    if (c[1] == 0) return 0;  // constant function
    *x0 = -c[0] / c[1];
    return 1;
  }

  // Quadratic solver
  int solve_quadratic(double c[3], double x0[2])
  {
    if (is_zero(c[2])) return solve_linear(c, x0);
    const double ainv = 1.0/c[2];       // ax^2 + bx + c = 0
    const double p = 0.5 * c[1] * ainv; // -b/2a
    const double q = c[0] * ainv;       // c/a
    double D = p*p-q;

    if (is_zero(D))  { // quadratic has one repeated root
      x0[0] = -p;
      return 1;
    }
    if (D > 0) {       // quadratic has two real roots
      D = sqrt(D);
      x0[0] =  D - p;
      x0[1] = -D - p;
      return 2;
    }
    return 0;          // quadratic has no real roots
  }

  // Cubic solver
  int solve_cubic(double c[4], double x0[3])
  {
    int num_roots;
    if (is_zero(c[3])) return solve_quadratic(c, x0);
    // normalize to  x^3 + Ax^2 + Bx + C = 0
    const double c3inv = 1.0/c[3];
    const double A = c[2] * c3inv;
    const double B = c[1] * c3inv;
    const double C = c[0] * c3inv;

    // substitute x = t - A/3 so t^3 + pt + q = 0
    const double A2 = A*A;
    const double p = (1.0/3.0)*((-1.0/3.0)*A2 + B);
    const double q = 0.5*((2.0/27.0)*A*A2 - (1.0/3.0)*A*B + C);

    // Cardano's fomula
    const double p3 = p*p*p;
    const double D  = q*q + p3;
    if (is_zero(D)) {
      if (is_zero(q)) { // one triple soln
        x0[0] = 0.0;
        num_roots = 1;
      }
      else {            // one single and one double soln
        const double u  = pow(fabs(q), 1.0/3.0)*sign(q);
        x0[0] = -2.0*u;
        x0[1] = u;
        num_roots = 2;
      }
    }
    else {
      if (D < 0.0) {    // three real roots
        const double phi = 1.0/3.0 * acos(-q/sqrt(-p3));
        const double t   = 2.0 * sqrt(-p);
        x0[0] =  t * cos(phi);
        x0[1] = -t * cos(phi + PI_OVER_3);
        x0[2] = -t * cos(phi - PI_OVER_3);
        num_roots = 3;
      }
      else {            // one real root
        const double sqrt_D = sqrt(D);
        const double u      = pow(sqrt_D + fabs(q), 1.0/3.0);
        if (q > 0) x0[0] = -u + p / u;
        else       x0[0] =  u - p / u;
        num_roots = 1;
      }
    }
    double sub = (1.0/3.0)*A;
    for (int i=0; i<num_roots; i++) x0[i] -= sub;
    return num_roots;
  }

  // solve ode with polynomial source : y'n + a_n-1 y'n-1 + ... = b_n x^n +...
  void integrate_ode(double x,
                     int na, double * a, double * y0, double * y, int nb, double * /* b */ )
  {
    if (na == 2) {
      // particular
      if ( a[1] == 0) {
        if ( a[0] == 0) {
          y[0] = y0[0]+y0[1]*x;
          y[1] =       y0[1];
        }
        else {
          double c = sqrt(a[0]);
          y[0] =    y0[0]*cos(c*x)+y0[1]/c*sin(c*x);
          y[1] = -c*y0[0]*cos(c*x)+y0[1]  *sin(c*x);
        }
      }
      else {
        // use solve_quadratic
        throw ATC_Error("not yet supported");
      }
      // homogenous
      double c = 1.;
      double z = x;
      int j = 2;
      for (int i = 0; i < nb; i++,j++) {
        y[1] += j*c*z;
        c /= j;
        z *= x;
        y[0] += c*z;
      }
    }
    else throw ATC_Error("can only integrate 2nd order ODEs currently");
  }
}
