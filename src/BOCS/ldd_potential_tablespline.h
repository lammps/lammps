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
    The Pennsylvania State University
   ------------------------------------------------------ */
#ifdef LDD_POTENTIAL_CLASS
// clang-format off
LddPotentialStyle(table/spline,LddPotentialTableSpline);
// clang-format on
#else

#ifndef LDD_POTENTIAL_TABLESPLINE
#define LDD_POTENTIAL_TABLESPLINE

#include "ldd_potential.h"

namespace LAMMPS_NS {

class LddPotentialTableSpline : public LddPotential {
 public:
  LddPotentialTableSpline(class LAMMPS *);
  ~LddPotentialTableSpline() override;

  void setup_potl(int, int, char **) override;
  double u(double) override;
  double f(double) override;
  double splint(double, double, double, double, double, double, double, double, double, double);

  void spline(double *, double *, int, double, double, double *);

 protected:
  virtual void allocate();
};
}    // namespace LAMMPS_NS

#endif
#endif
