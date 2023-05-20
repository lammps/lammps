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

#ifdef FIX_CLASS
// clang-format off
FixStyle(adapt/dcci,FixAdaptDCCI);
// clang-format on
#else

#ifndef LMP_FIX_ADAPT_DCCI_H
#define LMP_FIX_ADAPT_DCCI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdaptDCCI : public Fix {
 public:
  FixAdaptDCCI(class LAMMPS *, int, char **);
  ~FixAdaptDCCI() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void post_run() override;
  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  double lambda;
  double compute_scalar() override;

 private:
  int nadapt, fscaleflag;
  int anypair;
  int nlevels_respa;

  struct Adapt {
    int which;
    char *pstyle, *pparam;
    int ilo, ihi, jlo, jhi;
    int pdim;
    double *scalar, scalar_orig;
    double **array, **array_orig;
    class Pair *pair;
  };

  Adapt *adapt;

  void change_settings();
  void restore_settings();
};

}    // namespace LAMMPS_NS

#endif
#endif
