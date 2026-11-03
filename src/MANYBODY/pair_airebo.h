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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo,PairAIREBO);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_H
#define LMP_PAIR_AIREBO_H

#include "math_const.h"
#include "pair.h"
#include <cmath>

namespace LAMMPS_NS {

class PairAIREBO : public Pair {
 public:
  PairAIREBO(class LAMMPS *);
  ~PairAIREBO() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  enum { AIREBO, REBO_2, AIREBO_M };    // for telling class variants apart in shared code

 protected:
  int variant;
  int ljflag, torflag;    // 0/1 if LJ/Morse,torsion terms included
  int morseflag;          // 1 if Morse instead of LJ for non-bonded
  int bcflag;             // 1 for the bond-centric (AIREBO-BC) P_CC variant

  double cutlj;                     // user-specified LJ cutoff
  double sigcut, sigwid, sigmin;    // corresponding cutoff function
  double cutljrebosq;               // cut for when to compute
                                    // REBO neighs of ghost atoms

  double **cutljsq;                     // LJ cutoffs for C,H types
  double **lj1, **lj2, **lj3, **lj4;    // pre-computed LJ coeffs for C,H types
  double cut3rebo;                      // maximum distance for 3rd REBO neigh

  int maxlocal;             // size of numneigh, firstneigh arrays
  int pgsize;               // size of neighbor page
  int oneatom;              // max # of neighbors for one atom
  MyPage<int> *ipage;       // neighbor list pages
  int *REBO_numneigh;       // # of pair neighbors for each atom
  int **REBO_firstneigh;    // ptr to 1st neighbor of each atom

  double *closestdistsq;    // closest owned atom dist to each ghost
  double *nC, *nH;          // sum of weighting fns with REBO neighs

  double smin, Nmin, Nmax, NCmin, NCmax, thmin, thmax;
  double rcmin[2][2], rcmax[2][2], rcmaxsq[2][2], rcmaxp[2][2];
  double Q[2][2], alpha[2][2], A[2][2], rho[2][2], BIJc[2][2][3], Beta[2][2][3];
  double rcLJmin[2][2], rcLJmax[2][2], rcLJmaxsq[2][2], bLJmin[2][2], bLJmax[2][2];
  double epsilon[2][2], sigma[2][2], epsilonT[2][2];

  // parameters for Morse variant

  double epsilonM[2][2], alphaM[2][2], reqM[2][2];

  // spline coefficients

  double gCdom[5], gC1[4][6], gC2[4][6], gHdom[4], gH[3][6];
  double pCCdom[2][2], pCHdom[2][2], pCC[4][4][16], pCH[4][4][16];
  double piCCdom[3][2], piCHdom[3][2], piHHdom[3][2];
  double piCC[4][4][9][64], piCH[4][4][9][64], piHH[4][4][9][64];
  double Tijdom[3][2], Tijc[4][4][9][64];

  // spline knot values

  double PCCf[5][5], PCCdfdx[5][5], PCCdfdy[5][5], PCHf[5][5];
  double PCHdfdx[5][5], PCHdfdy[5][5];
  double piCCf[5][5][11], piCCdfdx[5][5][11];
  double piCCdfdy[5][5][11], piCCdfdz[5][5][11];
  double piCHf[5][5][11], piCHdfdx[5][5][11];
  double piCHdfdy[5][5][11], piCHdfdz[5][5][11];
  double piHHf[5][5][11], piHHdfdx[5][5][11];
  double piHHdfdy[5][5][11], piHHdfdz[5][5][11];
  double Tf[5][5][10], Tdfdx[5][5][10], Tdfdy[5][5][10], Tdfdz[5][5][10];

  // bond-centric P_CC (AIREBO-BC), used only when bcflag is set.  Stored on a
  // half-integer grid in doubled-coordinate index space (mC = 2*N_C, mH = 2*N_H).
  double pCCdom_bc[2][2];     // clamping domain in physical N
  double pCC_bc[6][6][16];    // bicubic patch coeffs (doubled coordinate)
  double PCCf_bc[7][7];       // knot values at the half-integer grid
  double PCCdfdx_bc[7][7];    // d/d(2*N_C) at knots (zero)
  double PCCdfdy_bc[7][7];    // d/d(2*N_H) at knots (zero)

  void REBO_neigh();
  void FREBO(int);
  void FLJ(int);
  void TORSION(int);

  double bondorder(int, int, double *, double, double, double **);
  double bondorderLJ(int, int, double *, double, double, double *, double, double **);

  double gSpline(double, double, int, double *, double *);
  double PijSpline(double, double, int, int, double *);

  // ----------------------------------------------------------------------
  // Bond-order P term evaluation.
  //
  // bondorder() and bondorderLJ() evaluate the P_ij correction through this
  // wrapper rather than calling PijSpline() directly.  Stock (atom-centric)
  // AIREBO forwards the i-side coordination numbers to PijSpline() (the
  // "other"-side counts are ignored).  The bond-centric variant (AIREBO-BC,
  // Hur & Stuart, J. Chem. Phys. 137, 054102 (2012); selected by bcflag)
  // instead evaluates P_CC at the bond-averaged coordination numbers
  // Nbar^t = 1/2 (N_ij^t + N_ji^t), reflecting both sides of the bond.
  //   thisC/thisH : coordination numbers on the side whose P is evaluated
  //   othC/othH   : coordination numbers on the opposite side of the bond
  // ----------------------------------------------------------------------
  double Pij_eval(double thisC, double thisH, double othC, double othH,
                  int typei, int typej, double dN2[2])
  {
    if (bcflag && typei == 0 && typej == 0)
      return PijSpline(0.5 * (thisC + othC), 0.5 * (thisH + othH), typei, typej, dN2);
    return PijSpline(thisC, thisH, typei, typej, dN2);
  }

  // whether the i-j bond uses a bond-averaged P (bond-centric C-C only); gates
  // the cross-side P force below.  False for stock atom-centric P.
  bool Pij_bond_averaged(int itype, int jtype) { return bcflag && itype == 0 && jtype == 0; }

  // ----------------------------------------------------------------------
  // Bond-centric P cross force.
  //
  // When Pij_eval() forms a bond-averaged P (bcflag), the P term depends on
  // the coordination on BOTH sides of the i-j bond and enters both pij and
  // pji.  The atom-centric P-coordination forces in bondorder() and
  // bondorderLJ() then miss two pieces -- the pji term acting on i's neighbors
  // and the pij term acting on j's neighbors -- which this helper adds.  For
  // stock atom-centric P, Pij_bond_averaged() is false and the helper is a
  // no-op, leaving AIREBO/REBO/AIREBO-M results bit-for-bit unchanged.
  // ----------------------------------------------------------------------
  void bondorder_Pij_cross(int i, int j, int itype, int jtype, double VA, double tmppij,
                           double tmppji, const double dN2PIJ[2], const double dN2PJI[2],
                           double **f);

  double piRCSpline(double, double, double, int, int, double *);
  double TijSpline(double, double, double, double *);

  void read_file(char *);

  // read the bond-centric P_CC knots appended after the standard AIREBO data
  // (bcflag only; a no-op otherwise).  Called from read_file() on rank 0 with
  // the reader positioned just past the Tij section.
  void read_file_extra(class PotentialFileReader &);

  double Spbicubic(double, double, double *, double *);
  double Sptricubic(double, double, double, double *, double *);
  void Sptricubic_patch_adjust(double *, double, double, char);
  void Sptricubic_patch_coeffs(double, double, double, double, double, double, double *, double *,
                               double *, double *, double *);
  void Spbicubic_patch_adjust(double *, double, double, char);
  void Spbicubic_patch_coeffs(double, double, double, double, double *, double *, double *,
                              double *);
  virtual void spline_init();

  void allocate();

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // added to header for inlining
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
   fifth order spline evaluation using Horner's rule
   ------------------------------------------------------------------------- */
  double Sp5th(const double &x, const double coeffs[6], double *df) const
  {
    double f = coeffs[5] * x;
    double d = 5.0 * coeffs[5] * x;
    f += coeffs[4];
    d += 4.0 * coeffs[4];
    f *= x;
    d *= x;
    f += coeffs[3];
    d += 3.0 * coeffs[3];
    f *= x;
    d *= x;
    f += coeffs[2];
    d += 2.0 * coeffs[2];
    f *= x;
    d *= x;
    f += coeffs[1];
    d += coeffs[1];
    f *= x;
    f += coeffs[0];
    *df = d;
    return f;
  }

  /* ----------------------------------------------------------------------
     cutoff function Sprime
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp(double Xij, double Xmin, double Xmax, double &dX) const
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = 0.5 * (1.0 + cos(t * MathConst::MY_PI));
      dX = (-0.5 * MathConst::MY_PI * sin(t * MathConst::MY_PI)) / (Xmax - Xmin);
    }
    return cutoff;
  };

  /* ----------------------------------------------------------------------
     LJ cutoff function Sp2
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp2(double Xij, double Xmin, double Xmax, double &dX) const
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = (1.0 - (t * t * (3.0 - 2.0 * t)));
      dX = 6.0 * (t * t - t) / (Xmax - Xmin);
    }
    return cutoff;
  };

  /* kronecker delta function returning a double */

  [[nodiscard]] double kronecker(const int a, const int b) const { return (a == b) ? 1.0 : 0.0; };
};
}    // namespace LAMMPS_NS

#endif
#endif
