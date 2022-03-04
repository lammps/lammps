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

ComputeStyle(pressure/spherical,ComputePressureSpherical)

#else

#ifndef LMP_COMPUTE_PRESSURE_SPHERICAL
#define LMP_COMPUTE_PRESSURE_SPHERICAL

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressureSpherical : public Compute {
 public:
  ComputePressureSpherical(class LAMMPS *, int, char **);
  ~ComputePressureSpherical();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int nbins;
  double bin_width, x0, y0, z0, Rmax;
  
  // Number density, kinetic and configurational contribution to the pressure.
  double *invV, *dens, *pkrr, *pktt, *pkpp, *pcrr, *pctt, *pcpp;
  double *tdens, *tpkrr, *tpktt, *tpkpp, *tpcrr, *tpctt, *tpcpp;
  class NeighList *list;
};

}

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

