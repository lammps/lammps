/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pressure/cylinder,ComputePressureCyl)

#else

#ifndef LMP_COMPUTE_PRESSURE_CYLINDER
#define LMP_COMPUTE_PRESSURE_CYLINDER

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressureCyl : public Compute {
 public:
  ComputePressureCyl(class LAMMPS *, int, char **);
  ~ComputePressureCyl();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int nbins,nphi,nzbins;
  double *Pr_temp,*Pr_all,*Pz_temp,*Pz_all,*Pphi_temp,*Pphi_all;
  double *R,*Rinv,*R2,*PrAinv,*PzAinv,PphiAinv;
  double Rmax,bin_width,nktv2p;
  double *R2kin,*density_temp,*invVbin,*density_all;
  double *tangent,*ephi_x,*ephi_y;
  double *binz;

  double zlo,zhi;

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

E: No pair style is defined for compute pressure/cylinder

Self-explanatory.

E: Pair style does not support compute pressure/cylinder

The pair style does not have a single() function, so it can
not be invoked by compute pressure/cylinder.

*/

