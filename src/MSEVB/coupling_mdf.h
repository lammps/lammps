/* -*- c++ -*- ----------------------------------------------------------
   Modified Damping Function (MDF) smooth cutoff utilities.
   Shared by all coupling models that support tapering.

   Reference: reference/mustard/src/mustard/utils.py MDF()
------------------------------------------------------------------------- */

#ifndef LMP_COUPLING_MDF_H
#define LMP_COUPLING_MDF_H

#include "domain.h"
#include <cmath>

namespace LAMMPS_NS {
namespace CouplingMDF {

  /* ----------------------------------------------------------------------
   Modified Damping Function (MDF) smooth cutoff.
   Returns 1 for r < rm, 0 for r > rc, smooth transition in between.
---------------------------------------------------------------------- */

  inline double MDF(double r, double rm, double rc)
  {
    if (r <= rm) return 1.0;
    if (r >= rc) return 0.0;
    double x = (r - rm) / (rc - rm);
    double omx = 1.0 - x;
    return omx * omx * omx * (1.0 + 3.0 * x + 6.0 * x * x);
  }

  /* ----------------------------------------------------------------------
   Scalar derivative dMDF/dr.
   Returns 0 outside [rm, rc].
---------------------------------------------------------------------- */

  inline double dMDF_dr(double r, double rm, double rc)
  {
    if (r <= rm || r >= rc) return 0.0;
    double dr = r - rm;
    double cr = rc - r;
    double range = rc - rm;
    double range5 = range * range * range * range * range;
    return -30.0 * dr * dr * cr * cr / range5;
  }

  /* ----------------------------------------------------------------------
   Apply MDF taper to a coupling value V and its forces fH, fX, fY.
   V is scaled by MDF(rHY); forces receive the product-rule correction.
   Pass fX = nullptr when there is no X-force contribution (e.g. Grimme).
   On entry V should contain V_raw; on exit V = V_raw * MDF.
---------------------------------------------------------------------- */

  inline void apply_taper(double *xH, double *xY, Domain *domain, double rm, double rc, double &V,
                          double *fH, double *fX, double *fY)
  {
    double dHY[3];
    dHY[0] = xH[0] - xY[0];
    dHY[1] = xH[1] - xY[1];
    dHY[2] = xH[2] - xY[2];
    domain->minimum_image(FLERR, dHY[0], dHY[1], dHY[2]);

    double rHY = std::sqrt(dHY[0] * dHY[0] + dHY[1] * dHY[1] + dHY[2] * dHY[2]);
    double mdf_val = MDF(rHY, rm, rc);
    double dmdf = dMDF_dr(rHY, rm, rc);
    double V_raw = V;
    V = V_raw * mdf_val;

    if (rHY > 0.0) {
      double inv_rHY = 1.0 / rHY;
      double taper_pf = dmdf * V_raw;
      for (int d = 0; d < 3; d++) {
        double uHY = dHY[d] * inv_rHY;
        fH[d] = fH[d] * mdf_val - taper_pf * uHY;
        fY[d] = fY[d] * mdf_val + taper_pf * uHY;
        if (fX) fX[d] *= mdf_val;
      }
    } else {
      for (int d = 0; d < 3; d++) {
        fH[d] *= mdf_val;
        fY[d] *= mdf_val;
        if (fX) fX[d] *= mdf_val;
      }
    }
  }

}    // namespace CouplingMDF
}    // namespace LAMMPS_NS

#endif
