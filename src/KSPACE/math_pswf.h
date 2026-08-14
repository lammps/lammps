#ifndef LMP_MATH_PSWF_H
#define LMP_MATH_PSWF_H

#include <unordered_map>
#include <vector>

namespace LAMMPS_NS::MathPSWF {

// prolate functions
void prolc180(double eps, double &c);

// prolate0 functor
struct Prolate0Fun;

double prolate0_eval_derivative(double c, double x);
/*
evaluate prolate0c at x, i.e., \psi_0^c(x)
*/
double prolate0_eval(double c, double x);

/*
evaluate prolate0c function integral of \int_0^r \psi_0^c(x) dx
*/
double prolate0_int_eval(double c, double r);

// approximation functions
void force_poly(double r_tol, const double &c, std::vector<double> &coeffs);
void energy_poly(double r_tol, const double &c, std::vector<double> &coeffs);
void fourier_poly(double r_tol, const double &c, double &lambda, std::vector<double> &coeffs);
void spread_fourier_poly(double r_tol, const double &c, double &lambda,
                         std::vector<double> &coeffs);
void spread_real_poly(int P, double r_tol, const double &c, std::vector<double> &coeffs);
}    // namespace LAMMPS_NS::MathPSWF
#endif    // LMP_MATH_PSWF_H
