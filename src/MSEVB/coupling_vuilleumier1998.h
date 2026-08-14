/* -*- c++ -*- ----------------------------------------------------------
   Vuilleumier1998 coupling for MSEVB: V = v12 * exp(-alpha * Q - gamma * q^2)

   Q = |r_XY|  (O-O distance)
   q = |r_H - midpoint(X,Y)|  (H to midpoint of donor/acceptor)

   Forces computed via chain rule with PBC minimum image.

   Header-only, no external dependencies. All arrays are raw pointers
   for Kokkos compatibility.

   Reference: Vuilleumier, Borgis (1998), Chem. Phys. Lett. 284, 71
   Python reference: reference/mustard/src/mustard/coupling.py class Vuilleumier1998
------------------------------------------------------------------------- */

#ifndef LMP_COUPLING_VUILLEUMIER1998_H
#define LMP_COUPLING_VUILLEUMIER1998_H

#include "coupling_mdf.h"
#include "domain.h"
#include <cmath>

namespace LAMMPS_NS {
namespace Vuilleumier1998 {

  /* ----------------------------------------------------------------------
   Compute Vuilleumier1998 coupling value and forces on atoms H, X, Y.

   v12, alpha, gamma: coupling parameters
   xH, xX, xY: atom positions (double[3])
   domain: LAMMPS domain for PBC minimum image
   V: output coupling value
   fH, fX, fY: output forces on atoms H, X, Y (force convention, i.e. -dV/dr)
---------------------------------------------------------------------- */

  inline void compute(double v12, double alpha, double gamma, double *xH, double *xX, double *xY,
                      Domain *domain, double &V, double *fH, double *fX, double *fY)
  {
    // displacement vector dXY = X - Y (with PBC)
    double dXY[3];
    dXY[0] = xX[0] - xY[0];
    dXY[1] = xX[1] - xY[1];
    dXY[2] = xX[2] - xY[2];
    domain->minimum_image(FLERR, dXY[0], dXY[1], dXY[2]);

    double Q = std::sqrt(dXY[0] * dXY[0] + dXY[1] * dXY[1] + dXY[2] * dXY[2]);

    // midpoint of X and Y (in PBC-aware sense: Y + 0.5 * dXY)
    // displacement from midpoint to H: dq = H - mid = H - Y - 0.5*dXY
    double dHY[3];
    dHY[0] = xH[0] - xY[0];
    dHY[1] = xH[1] - xY[1];
    dHY[2] = xH[2] - xY[2];
    domain->minimum_image(FLERR, dHY[0], dHY[1], dHY[2]);

    double dq[3];
    dq[0] = dHY[0] - 0.5 * dXY[0];
    dq[1] = dHY[1] - 0.5 * dXY[1];
    dq[2] = dHY[2] - 0.5 * dXY[2];

    double q_sq = dq[0] * dq[0] + dq[1] * dq[1] + dq[2] * dq[2];

    // Coupling value
    V = v12 * std::exp(-alpha * Q - gamma * q_sq);

    // Force from Q term: dV/dQ = -alpha * V
    // dQ/dr_X = dXY / Q, dQ/dr_Y = -dXY / Q
    double inv_Q = (Q > 0.0) ? 1.0 / Q : 0.0;

    // Python: derivOO = -alpha * cpl * dQpos / Q  (dQpos = X - Y direction)
    //         fOO = -derivOO
    //         f_X += fOO, f_Y -= fOO

    // Force from q term: dV/dq_i = -gamma * 2 * V * dq_i (component-wise, NOT normalized)
    // Python: derivHOO = -gamma * 2 * cpl * dqpos  (dqpos = H - midpoint direction)
    //         fHOO = derivHOO
    //         f_H -= fHOO, f_X += 0.5*fHOO, f_Y += 0.5*fHOO

    for (int d = 0; d < 3; d++) {
      double derivOO = -alpha * V * dXY[d] * inv_Q;
      double derivHOO = -gamma * 2.0 * V * dq[d];

      double fOO = -derivOO;
      double fHOO = derivHOO;

      fX[d] = fOO + 0.5 * fHOO;
      fY[d] = -fOO + 0.5 * fHOO;
      fH[d] = -fHOO;
    }
  }

  /* ----------------------------------------------------------------------
   Tapered coupling: V_tapered = V_raw * MDF(rHY, rm, rc)
   Forces get product-rule contribution from taper derivative.
   rm = taper start radius, rc = cutoff radius.
   Taper is based on H-Y distance (same as Raiteri2011).
---------------------------------------------------------------------- */

  inline void compute_tapered(double v12, double alpha, double gamma, double *xH, double *xX,
                              double *xY, Domain *domain, double rm, double rc, double &V,
                              double *fH, double *fX, double *fY)
  {
    compute(v12, alpha, gamma, xH, xX, xY, domain, V, fH, fX, fY);
    CouplingMDF::apply_taper(xH, xY, domain, rm, rc, V, fH, fX, fY);
  }

}    // namespace Vuilleumier1998
}    // namespace LAMMPS_NS

#endif
