/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

#ifndef SMD_KERNEL_FUNCTIONS_H_
#define SMD_KERNEL_FUNCTIONS_H_

namespace SMD_Kernels {
static inline double Kernel_Wendland_Quintic_NotNormalized(const double r, const double h)
{
  if (r < h) {
    double q = 2.0 * r / h;
    return pow(1.0 - 0.5 * q, 4) * (2.0 * q + 1.0);
  } else {
    return 0.0;
  }
}

static inline double Kernel_Cubic_Spline(const double r, const double h)
{
  double q = 2.0 * r / h;
  if (q > 2.0) {
    return 0.0;
  } else if ((q <= 2.0) && (q > 1.0)) {
    return pow(2.0 - q, 3.0) / 6.0;
  } else if ((q >= 0.0) && (q <= 1.0)) {
    return 2. / 3. - q * q + 0.5 * q * q * q;
  } else {
    return 0.0;
  }
}

static inline double Kernel_Barbara(const double r, const double h)
{
  double arg = (1.570796327 * (r + h)) / h;
  double hsq = h * h;
  //wf = (1.680351548 * (cos(arg) + 1.)) / hsq;
  return -2.639490040 * sin(arg) / (hsq * h);
}

static inline void spiky_kernel_and_derivative(const double h, const double r, const int dimension,
                                               double &wf, double &wfd)
{

  /*
         * Spiky kernel
         */

  if (r > h) {
    printf("r=%f > h=%f in Spiky kernel\n", r, h);
    wf = wfd = 0.0;
    return;
  }

  double hr = h - r;    // [m]
  if (dimension == 2) {
    double n = 0.3141592654e0 * h * h * h * h * h;    // [m^5]
    wfd = -3.0e0 * hr * hr / n;           // [m*m/m^5] = [1/m^3] ==> correct for dW/dr in 2D
    wf = -0.333333333333e0 * hr * wfd;    // [m/m^3] ==> [1/m^2] correct for W in 2D
  } else {
    wfd = -14.0323944878e0 * hr * hr /
        (h * h * h * h * h * h);          // [1/m^4] ==> correct for dW/dr in 3D
    wf = -0.333333333333e0 * hr * wfd;    // [m/m^4] ==> [1/m^3] correct for W in 3D
  }

  // alternative formulation
  //              double hr = h - r;
  //
  //              /*
  //               * Spiky kernel
  //               */
  //
  //              if (domain->dimension == 2) {
  //                      double h5 = h * h * h * h * h;
  //                      wf = 3.183098861e0 * hr * hr * hr / h5;
  //                      wfd = -9.549296583 * hr * hr / h5;
  //
  //              } else {
  //                      double h6 = h * h * h * h * h * h;
  //                      wf = 4.774648292 * hr * hr * hr / h6;
  //                      wfd = -14.32394487 * hr * hr / h6;
  //              }
  //      }
}

static inline void barbara_kernel_and_derivative(const double h, const double r,
                                                 const int dimension, double &wf, double &wfd)
{

  /*
         * Barbara kernel
         */

  double arg = (1.570796327 * (r + h)) / h;
  double hsq = h * h;

  if (r > h) {
    printf("r = %f > h = %f in barbara kernel function\n", r, h);
    exit(1);
    //wf = wfd = 0.0;
    //return;
  }

  if (dimension == 2) {
    wf = (1.680351548 * (cos(arg) + 1.)) / hsq;
    wfd = -2.639490040 * sin(arg) / (hsq * h);
  } else {
    wf = 2.051578323 * (cos(arg) + 1.) / (hsq * h);
    wfd = -3.222611694 * sin(arg) / (hsq * hsq);
  }
}

/*
 * compute a normalized smoothing kernel only
 */
static inline void Poly6Kernel(const double hsq, const double h, const double rsq,
                               const int dimension, double &wf)
{

  double tmp = hsq - rsq;
  if (dimension == 2) {
    wf = tmp * tmp * tmp / (0.7853981635e0 * hsq * hsq * hsq * hsq);
  } else {
    wf = tmp * tmp * tmp / (0.6382918409e0 * hsq * hsq * hsq * hsq * h);
  }
}

/*
 * M4 Prime Kernel
 */

static inline void M4PrimeKernel(const double s, double &wf)
{
  if (s < 1.0) {
    //wf = 1.0 - 2.5 * s * s + (3./2.) * s * s * s;
    wf = 1.0 - s * s * (2.5 - 1.5 * s);
  } else if (s < 2.0) {
    //wf = 0.5 * (1.0 - s) * ((2.0 - s) * (2.0 - s));
    wf = 2.0 + (-4.0 + (2.5 - 0.5 * s) * s) * s;
  } else {
    wf = 0.0;
  }
}

}    // namespace SMD_Kernels

#endif /* SMD_KERNEL_FUNCTIONS_H_ */
