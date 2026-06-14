// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   AIREBO-BC: bond-centric modification of the P_ij term of AIREBO.
   Reference: J. Hur and S. J. Stuart, "Modified reactive empirical
   bond-order potential for heterogeneous bonding environments",
   J. Chem. Phys. 137, 054102 (2012).  https://doi.org/10.1063/1.4738879
------------------------------------------------------------------------- */

#include "pair_airebo_bc.h"

#include "comm.h"
#include "error.h"
#include "potential_file_reader.h"

#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

static constexpr double TOL = 1.0e-9;

/* ---------------------------------------------------------------------- */

PairAIREBObc::PairAIREBObc(LAMMPS *lmp) : PairAIREBO(lmp)
{
  // all setup (flags, variant=AIREBO, neighbor bookkeeping, sig* cutoffs) is
  // performed by the PairAIREBO constructor; the bond-centric variant adds
  // no new heap allocations, so the inherited destructor suffices.
}

/* ----------------------------------------------------------------------
   bond-order P evaluation hook.

   For a C-C bond, evaluate P_CC at the bond-averaged coordination numbers
   Nbar^t = 1/2 (N_this^t + N_oth^t)  (Hur & Stuart 2012, Eqs. 4-5).  For all
   other bond types the bond-centric change does not apply and we defer to the
   standard atom-centric evaluation (this side's coordination numbers).
------------------------------------------------------------------------- */

double PairAIREBObc::Pij_eval(double thisC, double thisH, double othC, double othH,
                              int typei, int typej, double dN2[2])
{
  if (typei == 0 && typej == 0) {
    const double NbarC = 0.5 * (thisC + othC);
    const double NbarH = 0.5 * (thisH + othH);
    return PijSpline(NbarC, NbarH, typei, typej, dN2);
  }
  return PijSpline(thisC, thisH, typei, typej, dN2);
}

/* ----------------------------------------------------------------------
   P_ij spline.

   The C-C branch evaluates the bond-centric P_CC on the half-integer grid;
   every other case (C-H, and the typei==1 early return) is identical to
   PairAIREBO and is delegated to the base implementation.
------------------------------------------------------------------------- */

double PairAIREBObc::PijSpline(double NijC, double NijH, int typei, int typej,
                               double dN2[2])
{
  // only the carbon-carbon correction differs from stock AIREBO
  if (!(typei == 0 && typej == 0))
    return PairAIREBO::PijSpline(NijC, NijH, typei, typej, dN2);

  int x, y;
  double Pij = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;

  // ------------------------------------------------------------------
  // bond-centric P_CC.
  //
  // NijC, NijH arrive here as the *bond-averaged* coordination numbers
  //   Nbar^t = 1/2 (N_ij^t + N_ji^t)   (averaging done in Pij_eval()).
  //
  // The half-integer P_CC grid is stored in doubled-coordinate index space
  // (see spline_init): cell (mC,mH) spans [mC,mC+1]x[mH,mH+1] with
  // mC = 2*N_C, mH = 2*N_H.  We evaluate at u = 2*NijC, v = 2*NijH.
  //
  // Derivative scaling (matches the Fortran exactly): the force code expects
  // dP/dN_ij.  With bond averaging dP/dN_ij = (dP/dNbar)(1/2), and since Nbar
  // maps to u = 2*Nbar, dP/dNbar = 2 (dP/du), so dP/dN_ij = dP/du.  The factor
  // 2 (from u=2*Nbar) and 1/2 (from averaging) cancel, so Spbicubic's df
  // (= dP/du) is returned UNSCALED.
  // ------------------------------------------------------------------

  // clamp to the physical-N domain
  if (NijC < pCCdom_bc[0][0]) NijC = pCCdom_bc[0][0];
  if (NijC > pCCdom_bc[0][1]) NijC = pCCdom_bc[0][1];
  if (NijH < pCCdom_bc[1][0]) NijH = pCCdom_bc[1][0];
  if (NijH > pCCdom_bc[1][1]) NijH = pCCdom_bc[1][1];

  const double u = 2.0 * NijC;   // doubled coordinate
  const double v = 2.0 * NijH;

  x = (int) floor(u);
  y = (int) floor(v);

  if (fabs(u - floor(u)) < TOL && fabs(v - floor(v)) < TOL) {
    // exactly on a half-integer knot -> table lookup
    Pij    = PCCf_bc[x][y];
    dN2[0] = PCCdfdx_bc[x][y];   // derivatives at knots are zero
    dN2[1] = PCCdfdy_bc[x][y];
  } else {
    // interior of a cell -> evaluate the bicubic patch
    if (u == 2.0 * pCCdom_bc[0][1]) --x;   // upper edge belongs to last cell
    if (v == 2.0 * pCCdom_bc[1][1]) --y;
    Pij = Spbicubic(u, v, pCC_bc[x][y], dN2);   // dN2 = dP/du, dP/dv (unscaled)
  }

  return Pij;
}

/* ----------------------------------------------------------------------
   read the bond-centric P_CC knots that follow the standard AIREBO data.

   Called on rank 0 only (PairAIREBO::read_file), with the reader positioned
   just past the Tij section.  File format: an integer count N, then N triples
   "mC mH value" with mC = 2*N_C, mH = 2*N_H (integer indices 0..6).  The
   triples are read with next_dvector so the three tokens on a line are taken
   together (next_int/next_double are line-oriented).
------------------------------------------------------------------------- */

void PairAIREBObc::read_file_extra(PotentialFileReader &reader)
{
  for (int i = 0; i < 7; i++)
    for (int j = 0; j < 7; j++) PCCf_bc[i][j] = 0.0;

  int nbc = reader.next_int();
  if (nbc < 0 || nbc > 49)
    error->one(FLERR, "AIREBO-BC pCC knot count out of range: {}", nbc);
  std::vector<double> bcvals(3 * nbc);
  reader.next_dvector(bcvals.data(), 3 * nbc);
  for (int n = 0; n < nbc; n++) {
    int mC = (int) (bcvals[3*n]   + 0.5);
    int mH = (int) (bcvals[3*n+1] + 0.5);
    if (mC < 0 || mC > 6 || mH < 0 || mH > 6)
      error->one(FLERR, "AIREBO-BC pCC knot index out of range: {} {}", mC, mH);
    PCCf_bc[mC][mH] = bcvals[3*n+2];
  }
}

/* ----------------------------------------------------------------------
   build the spline coefficients.

   Run the standard AIREBO spline setup first, then broadcast the
   bond-centric P_CC knot values (read on rank 0) to all ranks and build the
   bicubic patches on the half-integer grid.
------------------------------------------------------------------------- */

void PairAIREBObc::spline_init()
{
  PairAIREBO::spline_init();

  // knot values were read on rank 0 in read_file_extra(); share them
  MPI_Bcast(&PCCf_bc[0][0], 49, MPI_DOUBLE, 0, world);

  // knot derivatives are zero; clamping domain is N in [0,3]
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 7; j++) {
      PCCdfdx_bc[i][j] = 0.0;
      PCCdfdy_bc[i][j] = 0.0;
    }
  }
  pCCdom_bc[0][0] = 0.0;  pCCdom_bc[0][1] = 3.0;   // N_C in [0,3]
  pCCdom_bc[1][0] = 0.0;  pCCdom_bc[1][1] = 3.0;   // N_H in [0,3]

  // build bicubic patch coefficients on the half-integer grid.
  // 6x6 cells; cell (mC,mH) spans index coords [mC,mC+1] x [mH,mH+1].
  // derivatives at the knots are set to zero (y1 = y2 = 0).
  for (int mH = 0; mH < 6; mH++) {
    for (int mC = 0; mC < 6; mC++) {
      double y[4] = {0}, y1[4] = {0}, y2[4] = {0};
      y[0] = PCCf_bc[mC][mH];
      y[1] = PCCf_bc[mC][mH+1];
      y[2] = PCCf_bc[mC+1][mH];
      y[3] = PCCf_bc[mC+1][mH+1];
      Spbicubic_patch_coeffs(mC, mC+1, mH, mH+1, y, y1, y2, &pCC_bc[mC][mH][0]);
    }
  }
}
