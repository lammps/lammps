/* -*- c++ -*- ----------------------------------------------------------
   Grimme2015 coupling for MSEVB: V = a * exp(-b * dE^2)

   Energy-based coupling: uses energy difference between EVB states.
   Forces on ALL atoms via chain rule through energy/force difference.

   dE = E_state_k - E_state_0
   V = a * exp(-b * dE^2)
   dV/ddE = -2*b*dE*V
   coupling_forces = dV/ddE * (F_state_0 - F_state_k)

   Header-only, no external dependencies. All arrays are raw pointers
   for Kokkos compatibility.

   Reference: Grimme, Brandenburg (2015), Phys. Chem. Chem. Phys. 17, 16992
   Python reference: reference/mustard/src/mustard/coupling.py class Grimme2015
------------------------------------------------------------------------- */

#ifndef LMP_COUPLING_GRIMME2015_H
#define LMP_COUPLING_GRIMME2015_H

#include "coupling_mdf.h"
#include "domain.h"
#include <cmath>

namespace LAMMPS_NS {
namespace Grimme2015 {

  /* ----------------------------------------------------------------------
   Compute Grimme2015 coupling value and fill off-diagonal dH block.

   a, b: coupling parameters
   E_parent: potential energy of parent state
   E_child: potential energy of child state
   np: number of partitions (Hamiltonian row stride)
   parent_state: parent EVB state index
   child_state: child EVB state index
   dH: force-derivative tensor (flat array)
   dH_stride: stride per Hamiltonian element (dH_nmax * 3)
   nlocal: number of local atoms
   V: output coupling value

   This function reads diagonal dH blocks [parent,parent] and [child,child]
   to get per-atom forces for each state, then fills off-diagonal blocks
   [parent,child] and [child,parent] with the coupling forces.
---------------------------------------------------------------------- */

  inline void compute(double a, double b, double E_parent, double E_child, int np, int parent_state,
                      int child_state, double *dH, int dH_stride, int nlocal, double &V)
  {
    double dE = E_child - E_parent;
    V = a * std::exp(-b * dE * dE);

    double dVddE = -2.0 * b * dE * V;

    int off_pp = (parent_state * np + parent_state) * dH_stride;
    int off_cc = (child_state * np + child_state) * dH_stride;
    int off_pc = (parent_state * np + child_state) * dH_stride;
    int off_cp = (child_state * np + parent_state) * dH_stride;

    for (int i = 0; i < nlocal * 3; i++) {
      double dF = dH[off_pp + i] - dH[off_cc + i];
      double coupling_force = -dVddE * dF;
      dH[off_pc + i] = coupling_force;
      dH[off_cp + i] = coupling_force;
    }
  }

  /* ----------------------------------------------------------------------
   Scalar-only variant: compute coupling value V and Grimme weight g
   without touching the dH tensor.  Used by weight-based force mixing.

   g = -2*b*dE*V  is the scalar such that off-diagonal coupling forces
   dH[p][c][a] = -g * (f_parent[a] - f_child[a]).
---------------------------------------------------------------------- */

  inline void compute_scalar(double a, double b, double E_parent, double E_child, double &V,
                             double &g)
  {
    double dE = E_child - E_parent;
    V = a * std::exp(-b * dE * dE);
    g = -2.0 * b * dE * V;
  }

  /* ----------------------------------------------------------------------
   Tapered scalar variant: V_tapered = V_raw * MDF(rHY, rm, rc)
   g_tapered = g_raw * MDF(rHY, rm, rc)
   taper_fH/fY: forces from taper derivative on H and Y atoms.
   rm = taper start radius, rc = cutoff radius.
---------------------------------------------------------------------- */

  inline void compute_scalar_tapered(double a, double b, double E_parent, double E_child,
                                     double *xH, double *xY, Domain *domain, double rm, double rc,
                                     double &V, double &g, double *taper_fH, double *taper_fY)
  {
    double dE = E_child - E_parent;
    double V_raw = a * std::exp(-b * dE * dE);
    double g_raw = -2.0 * b * dE * V_raw;

    double dHY[3];
    dHY[0] = xH[0] - xY[0];
    dHY[1] = xH[1] - xY[1];
    dHY[2] = xH[2] - xY[2];
    domain->minimum_image(FLERR, dHY[0], dHY[1], dHY[2]);

    double rHY = std::sqrt(dHY[0] * dHY[0] + dHY[1] * dHY[1] + dHY[2] * dHY[2]);

    double mdf_val = CouplingMDF::MDF(rHY, rm, rc);
    double dmdf = CouplingMDF::dMDF_dr(rHY, rm, rc);

    V = V_raw * mdf_val;
    g = g_raw * mdf_val;

    if (rHY > 0.0) {
      double inv_rHY = 1.0 / rHY;
      double taper_prefactor = dmdf * V_raw;
      for (int d = 0; d < 3; d++) {
        double uHY = dHY[d] * inv_rHY;
        taper_fH[d] = -taper_prefactor * uHY;
        taper_fY[d] = taper_prefactor * uHY;
      }
    } else {
      for (int d = 0; d < 3; d++) {
        taper_fH[d] = 0.0;
        taper_fY[d] = 0.0;
      }
    }
  }

}    // namespace Grimme2015
}    // namespace LAMMPS_NS

#endif
