/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */
#ifdef PAIR_CLASS
// clang-format off
PairStyle(mesocnt/viscous, PairMesoCNTViscous);
// clang-format on
#else

#ifndef LMP_PAIR_MESOCNT_VISCOUS_H
#define LMP_PAIR_MESOCNT_VISCOUS_H

#include "pair_mesocnt.h"

namespace LAMMPS_NS {
class PairMesoCNTViscous : public PairMesoCNT {
 public:
  using PairMesoCNT::PairMesoCNT;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

 protected:
  double fvisc_max, kvisc, vvisc, fvisc_shift;

  inline double weight(const double *, const double *, const double *, const double *);
  inline void weight(const double *, const double *, const double *, const double *, double &,
                     double *, double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
