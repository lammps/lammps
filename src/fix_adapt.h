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
FixStyle(adapt,FixAdapt);
// clang-format on
#else

#ifndef LMP_FIX_ADAPT_H
#define LMP_FIX_ADAPT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdapt : public Fix {
 public:
  int diamflag;    // 1 if atom diameters will vary, for AtomVecGranular
  int chgflag;

  FixAdapt(class LAMMPS *, int, char **);
  ~FixAdapt() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void post_run() override;
  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  void set_arrays(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  int nadapt, resetflag, scaleflag, massflag;
  int anypair, anybond, anyangle;
  int nlevels_respa;
  char *id_fix_diam, *id_fix_chg;
  class FixStoreAtom *fix_diam, *fix_chg;
  double previous_diam_scale, previous_chg_scale;
  int discflag;

  struct Adapt {
    int which, ivar;
    char *var;
    char *pstyle, *pparam;
    char *bstyle, *bparam;
    char *astyle, *aparam;
    int ilo, ihi, jlo, jhi;
    int pdim, bdim, adim;
    double *scalar, scalar_orig;
    double *vector, *vector_orig;
    double **array, **array_orig;
    int atomparam;
    class Pair *pair;
    class Bond *bond;
    class Angle *angle;
  };

  Adapt *adapt;
  double *kspace_scale;

  void change_settings();
  void restore_settings();
};

}    // namespace LAMMPS_NS

#endif
#endif
