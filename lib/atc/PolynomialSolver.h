#ifndef POLYNOMIAL_SOLVER_H
#define POLYNOMIAL_SOLVER_H
namespace ATC {
//* Solves a linear system, returns the number of roots found.
int solve_linear(double c[2], double x0[1]);
//* Solves a quadratic system, returns the number of roots found.
int solve_quadratic(double c[3], double x0[2]);
//* Solves a cubic system, returns the number of roots found.
int solve_cubic(double c[4], double x0[3]);
// solve y" = b1 x + b0, ics y0, y0'
int integrate_second_order_ode_with_linear_source(double b[2], double y0[2]);
}
#endif

