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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald,Ewald);
// clang-format on
#else

#ifndef LMP_EWALD_H
#define LMP_EWALD_H

#include "kspace.h"

namespace LAMMPS_NS {

class Ewald : public KSpace {
 public:
  Ewald(class LAMMPS *);
  ~Ewald() override;
  void init() override;
  void setup() override;
  void settings(int, char **) override;
  void compute(int, int) override;
  double memory_usage() override;

  void compute_group_group(int, int, int) override;

 protected:
  int kxmax, kymax, kzmax;
  int kcount, kmax, kmax3d, kmax_created;
  double gsqmx, volume;
  int nmax;

  double unitk[3];
  int *kxvecs, *kyvecs, *kzvecs;
  int kxmax_orig, kymax_orig, kzmax_orig;
  double *ug;
  double **eg, **vg;
  double **ek;
  double *sfacrl, *sfacim, *sfacrl_all, *sfacim_all;
  double ***cs, ***sn;

  // group-group interactions

  int group_allocate_flag;
  double *sfacrl_A, *sfacim_A, *sfacrl_A_all, *sfacim_A_all;
  double *sfacrl_B, *sfacim_B, *sfacrl_B_all, *sfacim_B_all;

  double rms(int, double, bigint, double);
  virtual void eik_dot_r();
  virtual void coeffs();
  virtual void allocate();
  virtual void deallocate();
  void slabcorr();

  // triclinic

  int triclinic;
  void eik_dot_r_triclinic();
  void coeffs_triclinic();

  // group-group interactions

  void slabcorr_groups(int, int, int);
  void allocate_groups();
  void deallocate_groups();
};

}    // namespace LAMMPS_NS

#endif
#endif
