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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pressure/cartesian,ComputePressureCartesian);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_CARTESIAN_H
#define LMP_COMPUTE_PRESSURE_CARTESIAN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressureCartesian : public Compute {
 public:
  ComputePressureCartesian(class LAMMPS *, int, char **);
  ~ComputePressureCartesian();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int nbins1, nbins2, dir1, dir2, dims;
  double bin_width1, bin_width2, invV;

  // Number density, kinetic and configurational contribution to the pressure.
  double *dens, *pkxx, *pkyy, *pkzz, *pcxx, *pcyy, *pczz;
  double *tdens, *tpkxx, *tpkyy, *tpkzz, *tpcxx, *tpcyy, *tpczz;
  class NeighList *list;
  void compute_pressure_1d(double, double, double, double, double, double);
  void compute_pressure_2d(double, double, double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No pair style is defined for compute pressure/cartesian

Self-explanatory.

E: Pair style does not support compute pressure/cartesian

The pair style does not have a single() function, so it can
not be invoked by compute pressure/cartesian.

*/
