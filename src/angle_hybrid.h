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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(hybrid,AngleHybrid);
// clang-format on
#else

#ifndef LMP_ANGLE_HYBRID_H
#define LMP_ANGLE_HYBRID_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleHybrid : public Angle {
 public:
  int nstyles;        // # of different angle styles
  Angle **styles;     // class list for each Angle style
  char **keywords;    // keyword for each Angle style

  AngleHybrid(class LAMMPS *);
  ~AngleHybrid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int) override;
  double memory_usage() override;

 private:
  int *map;    // which style each angle type points to

  int *nanglelist;     // # of angles in sub-style anglelists
  int *maxangle;       // max # of angles sub-style lists can store
  int ***anglelist;    // anglelist for each sub-style

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
