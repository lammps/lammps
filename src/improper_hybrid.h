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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(hybrid,ImproperHybrid);
// clang-format on
#else

#ifndef LMP_IMPROPER_HYBRID_H
#define LMP_IMPROPER_HYBRID_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperHybrid : public Improper {
 public:
  int nstyles;          // # of different improper styles
  Improper **styles;    // class list for each Improper style
  char **keywords;      // keyword for each improper style

  ImproperHybrid(class LAMMPS *);
  ~ImproperHybrid() override;
  void init_style() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double memory_usage() override;

 private:
  int *map;    // which style each improper type points to

  int *nimproperlist;     // # of impropers in sub-style improperlists
  int *maximproper;       // max # of impropers sub-style lists can store
  int ***improperlist;    // improperlist for each sub-style

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
