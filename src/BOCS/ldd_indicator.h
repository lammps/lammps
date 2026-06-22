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
#ifndef LMP_LDD_INDICATOR_H
#define LMP_LDD_INDICATOR_H

#include "pointers.h"

namespace LAMMPS_NS {

class LddIndicator : protected Pointers {
 public:
  double r0, rc;    // The min/max length scales (resp) where the indicator is defined [user]
  // extra length scales that may be relevent to the indicator [derived]
  // e.g. the radius of small and big sphere inferred from r0 and rc in ldd_indicator_sphere
  double rs, rb;
  double *coeffs;          // Coefficients involved in the indicator function w(r)
  double norm, invnorm;    // Normalization [w] and its inverse
  int n_coeffs;            // length of coeffs
  int allocated;           // 0 or 1, whether the indicator has been allocated

  LddIndicator(class LAMMPS *);
  ~LddIndicator() override;

  // Specific coeffs involved in w(r)
  virtual void init_coeffs(double, double, int) = 0;
  virtual double w(double) = 0;      // w(r) given the coeffs (overriden)
  virtual double wp(double) = 0;     // w'(r)
  virtual double wp2(double) = 0;    // w''(r)
};

}    // namespace LAMMPS_NS

#endif
