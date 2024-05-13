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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/profile,ComputeTempProfile);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_PROFILE_H
#define LMP_COMPUTE_TEMP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempProfile : public Compute {
 public:
  ComputeTempProfile(class LAMMPS *, int, char **);
  ~ComputeTempProfile() override;
  void init() override;
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;
  void compute_array() override;

  void reset_extra_dof() override;
  void remove_bias(int, double *) override;
  void remove_bias_thr(int, double *, double *) override;
  void remove_bias_all() override;
  void restore_bias(int, double *) override;
  void restore_bias_thr(int, double *, double *) override;
  void restore_bias_all() override;
  double memory_usage() override;

 private:
  int xflag, yflag, zflag, ncount, outflag;
  int nbinx, nbiny, nbinz, nbins;
  int ivx, ivy, ivz;
  double tfactor;
  double nstreaming;

  int box_change, triclinic;
  int *periodicity;
  double *boxlo, *boxhi, *prd;
  double invdelta[3];

  int maxatom;
  int *bin;
  double **vbin, **binave;
  double *tbin, *tbinall;

  void dof_compute();
  void bin_average();
  void bin_setup();
  void bin_assign();
};

}    // namespace LAMMPS_NS

#endif
#endif
