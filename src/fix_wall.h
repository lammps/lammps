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

#ifndef LMP_FIX_WALL_H
#define LMP_FIX_WALL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWall : public Fix {
 public:
  int nwall;
  int wallwhich[6];
  double coord0[6];
  int xflag;    // 1 if any wall position is a variable
  int xstyle[6];
  int xindex[6];
  char *xstr[6];

  FixWall(class LAMMPS *, int, char **);
  ~FixWall() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void pre_force(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

  virtual void precompute(int) = 0;
  virtual void wall_particle(int, int, double) = 0;

 protected:
  double epsilon[6], sigma[6], alpha[6], cutoff[6];
  double ewall[7], ewall_all[7];
  double xscale, yscale, zscale;
  int estyle[6], sstyle[6], astyle[6], wstyle[6];
  int eindex[6], sindex[6];
  char *estr[6], *sstr[6], *astr[6], *lstr[6], *fstr[6], *kstr[6];
  int varflag;    // 1 if any wall position,epsilon,sigma is a variable
  int eflag;      // per-wall flag for energy summation
  int ilevel_respa;
  int fldflag;
};

}    // namespace LAMMPS_NS

#endif
