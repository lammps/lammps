/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(runner,PairRuNNer);
// clang-format on
#else

#ifndef LMP_PAIR_RUNNER_H
#define LMP_PAIR_RUNNER_H

#include "lmptype.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairRuNNer : public Pair {
 public:
  PairRuNNer(class LAMMPS *);
  ~PairRuNNer() override;

  void compute(int eflag, int vflag) override;
  void settings(int narg, char **arg) override;
  void coeff(int narg, char **arg) override;
  void init_style() override;
  double init_one(int, int) override;
  void allocate();

  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) override;
  void unpack_forward_comm(int n, int first, double *buf) override;
  int pack_reverse_comm(int n, int first, double *buf) override;
  void unpack_reverse_comm(int n, int *list, double *buf) override;

  void pack_structure(int rank, int size, int natoms, int inum, int *ilist, tagint *tag, double **x,
                      int *runner_types, double *xyz_global, int *z_global);

  void pack_atomic_property(int rank, int size, int natoms, int inum, int *ilist, tagint *tag,
                            double *local_property, double *global_property);

  void unpack_local_atomic_properties(int rank, int size, int natoms, int inum, int *ilist,
                                      tagint *tag, int nprop, double *global_properties,
                                      double *local_properties);

  double memory_usage() override;

 private:
  double cflength;    // Length conversion factor.
  double cfenergy;    // Energy conversion factor.
  bool
      luse_prev_q;    // Use charges from previous timestep as initial guess for iterative qeq solvers.
  bool lwrite_f_comm;      // Write committee forces into f_comm array
  bool lwrite_q_comm;      // Write committee charges into q_comm array
  bool lcheck_extrap;      // Flag enabling checks for feature extrapolation
  bigint max_extrap;       // Maximal number of allowed timesteps with feature extrapolations
  bool lshow_ew;           // Flag enabling output of extrapolation warnings to log file
  bigint sum_ew_freq;      // Frequency where extrapolation warning summary is printed to log file
  bigint reset_ew_freq;    // Frequency where extrapolation count is reseted to 0
  bigint local_extrap_sum;    // Sum of recorded extrapolations per process over multiple time steps
  double cutoff;              // Max feature map cutoff.
  double total_charge;        // The total charge of the structure. Must be 0 for periodic systems.
  char *directory;            // directory containing RuNNer potential files
  int *map;                   // Mapping from atom types to elements
  int nmax;                   // Allocated size of per-atom arrays.
  static int instances;       // count pair style instances, since we currently
                              // only support one instance at a time

  // Additional per-atom arrays
  double *atomic_charge, *hirshfeld_volume, *electronegativity, *lagrange_charges, *de_dq,
      *screening_de_dq, *committee_storage;
  bool lhirshfeld_vdw;
  bool ltwo_body;
  int nnp_generation;
  int num_committee_members;    // specified in input.nn
  int commstyle;                // communication flag for forward and reverse communication
};
}    // namespace LAMMPS_NS
#endif
#endif
