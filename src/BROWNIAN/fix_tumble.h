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
FixStyle(tumble,FixTumble);
// clang-format on
#else

#ifndef LMP_FIX_TUMBLE_H
#define LMP_FIX_TUMBLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTumble : public Fix {
 public:
  FixTumble(class LAMMPS *, int, char **);
  ~FixTumble() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  void post_integrate_respa(int, int) override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  int mode;             // per-atom property that stores the orientation
  double rate;          // tumbling rate (inverse time units)
  int seed;             // random number generator seed
  int planar_flag;      // 1 if new orientations are confined to the xy plane in 3d
  double ptumble;       // probability of a tumble during one timestep
  int nlevels_respa;    // number of rRESPA levels
  class RanMars *rng;

  void tumble_dipole();
};

}    // namespace LAMMPS_NS
#endif
#endif
