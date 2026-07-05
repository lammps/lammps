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

#ifndef LMP_EWALD_CONST_H
#define LMP_EWALD_CONST_H

namespace LAMMPS_NS::EwaldConst {

static constexpr double EWALD_F = 1.12837917;
static constexpr double EWALD_P = 0.3275911;
static constexpr double A1 = 0.254829592;
static constexpr double A2 = -0.284496736;
static constexpr double A3 = 1.421413741;
static constexpr double A4 = -1.453152027;
static constexpr double A5 = 1.061405429;

// bits in ewald_order / ewald_off flag which interaction order is treated
// with the long-range (k-space) solver: order 1 = Coulomb (1/r),
// order 3 = dipole, order 6 = LJ dispersion (1/r^6).

enum { EWALD_COUL = 1 << 1, EWALD_DIPOLE = 1 << 3, EWALD_DISP = 1 << 6 };

// highest interaction order that can be requested via ewald_order

static constexpr int EWALD_MAXORDER = 6;

// the dispersion-capable k-space solvers (ewald/disp, pppm/disp) compute up
// to four kinds of contributions ("terms"); their termflag[] array records
// which ones are enabled.  the first three indices are common to both
// solvers; index 3 is solver-specific: long-range point dipoles (ewald/disp)
// or unmixed per-type-pair dispersion (pppm/disp).

enum {
  TERM_COUL = 0,          // Coulomb 1/r
  TERM_DISP_GEOM = 1,     // dispersion 1/r^6, geometric mixing
  TERM_DISP_ARITH = 2,    // dispersion 1/r^6, arithmetic mixing
  TERM_DIPOLE = 3,        // point dipoles (ewald/disp only)
  TERM_DISP_NONE = 3,     // dispersion 1/r^6, no mixing (pppm/disp only)
};
static constexpr int EWALD_NTERMS = 4;

}    // namespace LAMMPS_NS::EwaldConst

#endif
