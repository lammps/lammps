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
BondStyle(hybrid,BondHybrid);
// clang-format on
#else

#ifndef LMP_BOND_HYBRID_H
#define LMP_BOND_HYBRID_H

#include "bond.h"

namespace LAMMPS_NS {

class BondHybrid : public Bond {
  friend class Force;

 public:
  int nstyles;        // # of different bond styles
  Bond **styles;      // class list for each Bond style
  char **keywords;    // keyword for each Bond style

  BondHybrid(class LAMMPS *);
  ~BondHybrid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, double, int, int, double &) override;
  double memory_usage() override;

 private:
  int *map;           // which style each bond type points to
  int has_quartic;    // which style, if any is a quartic bond style
  int *nbondlist;     // # of bonds in sub-style bondlists
  int *maxbond;       // max # of bonds sub-style lists can store
  int ***bondlist;    // bondlist for each sub-style

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
