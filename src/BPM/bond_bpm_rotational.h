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

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/rotational,BondBPMRotational);
// clang-format on
#else

#ifndef LMP_BOND_BPM_ROTATIONAL_H
#define LMP_BOND_BPM_ROTATIONAL_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMRotational : public BondBPM {
 public:
  BondBPMRotational(class LAMMPS *);
  ~BondBPMRotational() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *Kr, *Ks, *Kt, *Kb, *gnorm, *gslide, *groll, *gtwist;
  double *Fcr, *Fcs, *Tct, *Tcb;
  int smooth_flag;

  double acos_limit(double);

  double elastic_forces(int, int, int, double &, double, double, double, double *, double *,
                        double *, double *, double *, double *);
  void damping_forces(int, int, int, double &, double *, double *, double *, double *, double *);

  void allocate();
  void store_data();
  double store_bond(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
