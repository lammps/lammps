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
PairStyle(airebo/bc,PairAIREBObc);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_BC_H
#define LMP_PAIR_AIREBO_BC_H

#include "pair_airebo.h"

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
   AIREBO-BC: bond-centric modification of the AIREBO bond-order potential
   (Hur & Stuart, J. Chem. Phys. 137, 054102 (2012)).

   The only physical change relative to AIREBO is that the P_ij coordination
   correction is made bond-centric: P_CC becomes a function of the
   bond-averaged coordination numbers Nbar^t = 1/2 (N_ij^t + N_ji^t) and is
   stored on a half-integer spline grid (Table III of the paper).  All other
   terms (g, pi^rc, pi^dh, LJ, torsion, neighbor handling, forces) are
   inherited unchanged from PairAIREBO.  This class therefore overrides only:

     - Pij_eval()        : form Nbar for C-C bonds before evaluating P
     - PijSpline()       : evaluate the bond-centric P_CC half-integer grid
     - spline_init()     : broadcast + build the bond-centric P_CC patches
     - read_file_extra() : read the bond-centric P_CC knots from the file
------------------------------------------------------------------------- */

class PairAIREBObc : public PairAIREBO {
 public:
  PairAIREBObc(class LAMMPS *);

 protected:
  // bond-centric P_CC on a half-integer grid.  Indexed by the doubled
  // coordination numbers mC = 2*N_C, mH = 2*N_H, with N_C, N_H in [0,3]
  // (=> m in [0,6]).  Patches cover m in [0,5] (6x6 cells).  This mirrors
  // the Fortran xh(icarb, 2*N_H+1, 2*N_C+1) table on 0-based C++ indices.
  double pCCdom_bc[2][2];     // clamping domain in physical N (not doubled)
  double pCC_bc[6][6][16];    // bicubic patch coeffs (doubled coordinate)
  double PCCf_bc[7][7];       // knot values at the half-integer grid
  double PCCdfdx_bc[7][7];    // d/d(2*N_C) at knots (zero)
  double PCCdfdy_bc[7][7];    // d/d(2*N_H) at knots (zero)

  // bond-centric overrides of the AIREBO bond-order machinery
  double Pij_eval(double thisC, double thisH, double othC, double othH,
                  int typei, int typej, double dN2[2]) override;
  // only C-C uses a bond-averaged P, so only C-C needs the cross force
  bool Pij_bond_averaged(int itype, int jtype) override { return itype == 0 && jtype == 0; }
  double PijSpline(double, double, int, int, double *) override;
  void spline_init() override;
  void read_file_extra(class PotentialFileReader &reader) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
