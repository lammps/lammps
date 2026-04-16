/* -*- c++ -*- ----------------------------------------------------------
   Raiteri2011 coupling for MSEVB: V = lambda * exp(-zeta * Q^2)

   Q = |r_HY| - |r_HX|, where H is the transferring proton,
   X is the donor oxygen, Y is the acceptor oxygen.

   Forces computed via chain rule with PBC minimum image.

   Header-only, no external dependencies. All arrays are raw pointers
   for Kokkos compatibility.

   Reference: Raiteri, Gale, Bussi (2011), J. Phys.: Condens. Matter 23, 334213
   Python reference: reference/mustard/src/mustard/coupling.py class Raiteri2011
------------------------------------------------------------------------- */

#ifndef LMP_COUPLING_RAITERI2011_H
#define LMP_COUPLING_RAITERI2011_H

#include "coupling_mdf.h"
#include "domain.h"
#include <cmath>

namespace LAMMPS_NS {
namespace Raiteri2011 {

  /* ----------------------------------------------------------------------
   Compute Raiteri2011 coupling value and forces on atoms H, X, Y.

   lambda, zeta: coupling parameters
   xH, xX, xY: atom positions (double[3])
   domain: LAMMPS domain for PBC minimum image
   V: output coupling value
   fH, fX, fY: output forces on atoms H, X, Y (force convention, i.e. -dV/dr)
---------------------------------------------------------------------- */

  inline void compute(double lambda, double zeta, double *xH, double *xX, double *xY,
                      Domain *domain, double &V, double *fH, double *fX, double *fY)
  {
    // displacement vectors with PBC: dHX = H - X, dHY = H - Y
    double dHX[3], dHY[3];

    dHX[0] = xH[0] - xX[0];
    dHX[1] = xH[1] - xX[1];
    dHX[2] = xH[2] - xX[2];
    domain->minimum_image(FLERR, dHX[0], dHX[1], dHX[2]);

    dHY[0] = xH[0] - xY[0];
    dHY[1] = xH[1] - xY[1];
    dHY[2] = xH[2] - xY[2];
    domain->minimum_image(FLERR, dHY[0], dHY[1], dHY[2]);

    // distances
    double rHX = std::sqrt(dHX[0] * dHX[0] + dHX[1] * dHX[1] + dHX[2] * dHX[2]);
    double rHY = std::sqrt(dHY[0] * dHY[0] + dHY[1] * dHY[1] + dHY[2] * dHY[2]);

    // Q = rHY - rHX (positive when proton is closer to donor X)
    double Q = rHY - rHX;

    // coupling value
    V = lambda * std::exp(-zeta * Q * Q);

    // dV/dQ = -2 * zeta * Q * V
    double prefactor = -2.0 * zeta * Q * V;

    // Unit vectors: uHX = dHX / rHX, uHY = dHY / rHY
    double inv_rHX = 1.0 / rHX;
    double inv_rHY = 1.0 / rHY;

    // Chain rule:
    //   dV/dr_H(via HY) = prefactor * (dHY / rHY)    [from d(rHY)/dr_H = dHY/rHY]
    //   dV/dr_H(via HX) = prefactor * -(dHX / rHX)   [from dQ/d(rHX) = -1]
    //
    // Force = -dV/dr, so:
    //   f_H(via HY) = -prefactor * (dHY / rHY)
    //   f_H(via HX) = +prefactor * (dHX / rHX)
    //   f_Y = -f_H(via HY) = +prefactor * (dHY / rHY)   [Newton's 3rd]
    //   f_X = -f_H(via HX) = -prefactor * (dHX / rHX)   [Newton's 3rd]

    for (int d = 0; d < 3; d++) {
      double derivHY = prefactor * dHY[d] * inv_rHY;
      double derivHX = -prefactor * dHX[d] * inv_rHX;

      fH[d] = -derivHY - derivHX;
      fY[d] = derivHY;
      fX[d] = derivHX;
    }
  }

  /* ----------------------------------------------------------------------
   Tapered coupling: V_tapered = V_raw * MDF(rHY, rm, rc)
   Forces get product-rule contribution from taper derivative.
   rm = taper start radius, rc = cutoff radius.
---------------------------------------------------------------------- */

  inline void compute_tapered(double lambda, double zeta, double *xH, double *xX, double *xY,
                              Domain *domain, double rm, double rc, double &V, double *fH,
                              double *fX, double *fY)
  {
    compute(lambda, zeta, xH, xX, xY, domain, V, fH, fX, fY);
    CouplingMDF::apply_taper(xH, xY, domain, rm, rc, V, fH, fX, fY);
  }

}    // namespace Raiteri2011
}    // namespace LAMMPS_NS

#endif
