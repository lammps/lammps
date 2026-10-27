/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_XL_QEQ_H
#define LMP_XL_QEQ_H

#include "lmptype.h"

namespace LAMMPS_NS {
class LAMMPS;

/* ----------------------------------------------------------------------
   settings and warm-up state of the extended-Lagrangian (XL) propagation
   of the auxiliary charge variables, shared by the charge equilibration
   fixes of the QEQ and REAXFF packages and their accelerated versions.
   this is a plain struct, so it can be copied into KOKKOS kernels.
------------------------------------------------------------------------- */

struct XLQEq {
  static constexpr int MAXHIST = 10;    // maximum depth of the auxiliary history

  int flag = 0;                 // 1 if the XL propagation is enabled
  int ncg = 2;                  // CG iterations per solve in XL mode: 2 keep the energy
                                // conservation close to converged solves, 1 is the count
                                // of the original publication. 1 and 3 were less stable
  int kdis = -1;                // dissipation order (0 = plain Verlet, -1 = not set)
  int kappa_set = 0;            // 1 if the xlkappa keyword was used
  double kappa = 2.0;           // coupling constant kappa = omega^2 * dt^2
  double alpha = 0.0;           // dissipation prefactor alpha
  double c[MAXHIST] = {0.0};    // dissipation coefficients c_0 ... c_kdis
  int nhist = 0;                // depth of the auxiliary history, 0 unless enabled
  int nseed = 0;                // consecutive converged solves seeding the history
  int bypass = 0;               // force a fully converged solve when set
  bigint laststep = -1;         // last step the auxiliary history was advanced
  bigint natoms = -1;           // number of atoms when the history was advanced
  double guard = 1.0e-3;        // relative residual of a truncated solve above which the
                                // solve is completed and the propagation restarted, <= 0: off
  bigint nfallback = 0;         // number of such restarts

  // parse keywords xlcg, xldamp, xlkappa, and xlguard: return the number of consumed arguments or 0
  int parse_keyword(LAMMPS *lmp, const char *style, int narg, char **arg, int iarg);

  // check the settings, apply defaults, set the coefficients and history depth
  void finalize(LAMMPS *lmp, const char *style, int nevery);

  // true if the auxiliary history can seed a truncated solve on the current step
  [[nodiscard]] bool ready(const LAMMPS *lmp) const;

  // true if the history must be (re-)seeded with all slots set to the current solution
  [[nodiscard]] bool restart(const LAMMPS *lmp) const;

  // update the warm-up state after the history was propagated (active) or seeded
  void advance(const LAMMPS *lmp, bool active, bool restart);

  // true if the relative residual of a truncated solve indicates that the
  // propagation became unstable and the solve must be completed to convergence
  [[nodiscard]] bool unstable(double resid) const { return (guard > 0.0) && (resid > guard); }

  // restart the warm-up after a truncated solve was completed to convergence
  void fallback(LAMMPS *lmp, const char *style, double resid, int warn);

  // time-reversible (dissipative) Verlet propagation of one auxiliary variable:
  //   theta(n+1) = 2*theta(n) - theta(n-1) + kappa*(x(n) - theta(n))
  //                + alpha*sum_k c_k theta(n-k)
  // hist holds theta(n), theta(n-1), ... (newest first)

  [[nodiscard]] double propagate(double x, const double *hist) const
  {
    double theta = 2.0 * hist[0] - hist[1] + kappa * (x - hist[0]);
    if (kdis)
      for (int k = 0; k <= kdis; ++k) theta += alpha * c[k] * hist[k];
    return theta;
  }

  // store a new value as the newest entry of the history
  void push(double *hist, double value) const
  {
    for (int k = nhist - 1; k > 0; --k) hist[k] = hist[k - 1];
    hist[0] = value;
  }

  // set all entries of the history to the same value
  void fill(double *hist, double value) const
  {
    for (int k = 0; k < nhist; ++k) hist[k] = value;
  }
};
}    // namespace LAMMPS_NS

#endif
