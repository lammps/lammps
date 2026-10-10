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
FixStyle(propel/ou,FixPropelOU);
// clang-format on
#else

#ifndef LMP_FIX_PROPEL_OU_H
#define LMP_FIX_PROPEL_OU_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropelOU : public Fix {
 public:
  FixPropelOU(class LAMMPS *, int, char **);
  ~FixPropelOU() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void reset_dt() override;

  void write_restart(FILE *) override;
  void restart(char *) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  double magnitude;      // root mean square value of each active force component
  double tau;            // correlation time of the active force
  int seed;              // random number generator seed
  double decay, kick;    // prefactors of the exact update over one timestep
  int ilevel_respa;      // rRESPA level at which the force is applied
  double **factive;      // per-atom active force
  class RanMars *rng;

  void advance();
  void apply_force(int);
};

}    // namespace LAMMPS_NS
#endif
#endif
