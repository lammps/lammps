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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#ifdef LDD_INDICATOR_CLASS
// clang-format off
LddIndicatorStyle(lucy,LddIndicatorLucy)
// clang-format on
#else

#ifndef LMP_LDD_INDICATOR_LUCY_H
#define LMP_LDD_INDICATOR_LUCY_H

#include "ldd_indicator.h"

namespace LAMMPS_NS {

class LddIndicatorLucy : public LddIndicator {
 public:
  LddIndicatorLucy(class LAMMPS *);
  ~LddIndicatorLucy() override;

  void init_coeffs(double, double, int) override;
  double w(double) override;
  double wp(double) override;
  double wp2(double) override;

 protected:
  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
