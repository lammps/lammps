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

/* ----------------------------------------------------------------------
   Contributing author: Agilio Padua (ENS de Lyon & CNRS)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(fep,ComputeFEP);
// clang-format on
#else

#ifndef COMPUTE_FEP_H
#define COMPUTE_FEP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFEP : public Compute {
 public:
  ComputeFEP(class LAMMPS *, int, char **);
  ~ComputeFEP() override;
  void init() override;
  void compute_vector() override;

 private:
  int npert;
  int pairflag;
  int chgflag;
  int tailflag, volumeflag;
  int fepinitflag;
  int eflag, vflag;
  double temp_fep;

  int nmax;
  double *q_orig;
  double **f_orig;
  double eng_vdwl_orig, eng_coul_orig;
  double pvirial_orig[6];
  double *peatom_orig, **pvatom_orig;
  double energy_orig;
  double kvirial_orig[6];
  double *keatom_orig, **kvatom_orig;

  class Fix *fixgpu;

  struct Perturb {
    int which, ivar;
    char *var;
    char *pstyle, *pparam;
    int ilo, ihi, jlo, jhi;
    int pdim;
    double **array, **array_orig;
    int aparam;
  };

  Perturb *perturb;

  double compute_epair();
  void perturb_params();
  void backup_params();
  void restore_params();
  void allocate_storage();
  void deallocate_storage();
  void backup_qfev();
  void restore_qfev();
};

}    // namespace LAMMPS_NS

#endif
#endif
