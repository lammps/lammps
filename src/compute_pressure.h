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
ComputeStyle(pressure,ComputePressure);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_H
#define LMP_COMPUTE_PRESSURE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePressure : public Compute {
 public:
  ComputePressure(class LAMMPS *, int, char **);
  ~ComputePressure() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  void reset_extra_compute_fix(const char *) override;

 protected:
  double boltz, nktv2p, inv_volume;
  int nvirial, dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];    // ordering: xx,yy,zz,xy,xz,yz
  int pairhybridflag;
  class Pair *pairhybrid;
  int keflag, pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int fixflag, kspaceflag;

  void virial_compute(int, int);

 private:
  char *pstyle;
  int nsub;
};

}    // namespace LAMMPS_NS

#endif
#endif
