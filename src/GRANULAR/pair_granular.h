/* -*- c++ -*- ----------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(granular,PairGranular);
// clang-format on
#else

#ifndef LMP_PAIR_GRANULAR_H
#define LMP_PAIR_GRANULAR_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class PairGranular : public Pair {
 public:
  PairGranular(class LAMMPS *);
  ~PairGranular() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void reset_dt() override;
  double single(int, int, int, int, double, double, double, double &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  double atom2cut(int) override;
  double radii2cut(double, double) override;

 protected:
  int freeze_group_bit;
  int use_history;

  int neighprev;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;
  double **cut;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, null pointer if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void allocate();
  void transfer_history(double *, double *, int, int) override;
  void prune_models();

 private:
  int size_history;
  int heat_flag;

  // granular models
  int nmodels, maxmodels;
  class Granular_NS::GranularModel** models_list;
  int **types_indices;

  // optional user-specified global cutoff, per-type user-specified cutoffs
  double **cutoff_type;
  double cutoff_global;
};

}    // namespace LAMMPS_NS

#endif
#endif
