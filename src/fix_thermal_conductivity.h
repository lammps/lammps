/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_THERMAL_CONDUCTIVITY_H
#define FIX_THERMAL_CONDUCTIVITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixThermalConductivity : public Fix {
 public:
  FixThermalConductivity(class LAMMPS *, int, char **);
  ~FixThermalConductivity();
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();

 private:
  int me;
  int edim,nbin,periodicity;
  int nswap;
  double prd,boxlo,boxhi;
  double slablo_lo,slablo_hi,slabhi_lo,slabhi_hi;
  double e_exchange;

  int nlo,nhi;
  int *index_lo,*index_hi;
  double *ke_lo,*ke_hi;
};

}

#endif
