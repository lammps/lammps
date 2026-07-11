// clang-format off
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

#ifndef LMP_TUNE_KOKKOS_H
#define LMP_TUNE_KOKKOS_H

#include "pointers.h"

namespace LAMMPS_NS {

class TuneKokkos : protected Pointers {
 public:
  TuneKokkos(class LAMMPS *, int kernel_type, int nevery, int nparams=2,
    const char* name=nullptr, int nsamples=5, int mode=0, double rel_tol=0.2);
  ~TuneKokkos() override;
  void allocate(int);
  void tuning_kernel_params();
  int get_current_combination_idx() { return combination_idx; }
  void set_relative_tolerance(double tol) { relative_tolerance = tol; }

  // override the default parameter value sets
  void set_team_size_values(const std::vector<int>& tsizes);
  void set_vector_size_values(const std::vector<int>& vsizes);

  enum { PAIR, BOND, NBIN, FIX, COMPUTE, GENERIC };

  int kernel_type;  // type of kernel being tuned: PAIR, BOND, NBIN, FIX, COMPUTE
  int interval;     // # of timesteps to run the simulation with a given parameter set
                    //   should be sufficiently large (~100) for estimating the overall performance in TPS
                    //   should be small enough to avoid significant overhead during auto-tuning

  std::vector<int> team_sizes;   // parameter values for the team size (typically, thread block size)
  std::vector<int> vector_sizes; // parameter values for the vector size (the 2nd dimension of thread block)

  int num_params;            // number of parameters to tune: 1 (team size only) or 2 (team size and threads per atom)
  double** performance;      // array to store the performance data for each parameter set
  int ncombinations;         // total number of parameter combinations
  int combination_idx;       // current combination index during scanning
  int scanning_completed;    // 0 if still scanning, 1 if scanning completed
  int allocated;             // 1 if the performance array is allocated and param values set up
  int nsamples;              // number of samples collected for each parameter combination
  int sample_idx;            // current sample index for the current parameter combination
  int mode;                  // how to determine the optimal performance among multiple samples for each parameter combination: MAX_VALUE, AVERAGE_VALUE, or MEDIAN_VALUE
  double opt_perf;           // stored the optimal performance
  int opt_combination_idx;   // stored the index of the optimal parameter combination
  double relative_tolerance; // acceptable threshold for performance degradation wrt opt_perf
  int nperf_degraded;        // number of consecutive times the performance is degraded beyond the relative tolerance

  bigint last_step;          // last timestep when timing info was collected
  double last_cpu;           // last CPU time when timing info was collected
  int firststep;             // 1 if first timestep for timing info collection

  FILE *tuning_logfile;      // logfile
  std::string my_name;       // name of the tuner instance

  double get_timing_info();                   // get the elapsed time from the last call
  int get_current_team_size();                // get the team size for the current combination index
  int get_current_vector_size();              // get the vector size for the current combination index
  void get_current_params(int, auto&, auto&); // get the team size and vector size for a given combination index
  void set_param_values(int);                 // set the KOKKOS kernel parameters based on the combination index
  int get_optimal_combination_idx();          // find the combination index with the best performance
  void regular_performance_check();           // monitor the performance during normal simulation and re-trigger scanning if needed
};

}

#endif

