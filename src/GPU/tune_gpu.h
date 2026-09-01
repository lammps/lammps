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

#ifndef LMP_TUNE_GPU_H
#define LMP_TUNE_GPU_H

#include "pointers.h"

#include <vector>

namespace LAMMPS_NS {

class TuneGPU : protected Pointers {
 public:
  TuneGPU(class LAMMPS *, int nevery, int nsamples = 5, int mode = 0,
          double rel_tol = 0.2, int max_nthreads = 0);
  ~TuneGPU() override;

  void tuning_kernel_params();

  int interval;                 // # of timesteps a parameter set is timed for
                                //   large enough (~100) for a meaningful rate,
                                //   small enough to keep the scan affordable

  std::vector<int> tpa_values;      // candidate values for threads per atom
  std::vector<int> block_values;    // candidate values for the pair block size
  std::vector<int> nthreads_values; // candidate values for host OpenMP threads

  double **performance;      // performance data for each parameter set
  int ncombinations;         // total number of parameter combinations
  int combination_idx;       // current combination index during scanning
  int scanning_completed;    // 0 while still scanning, 1 when done
  int nsamples;              // number of samples per parameter combination
  int sample_idx;            // current sample index
  int mode;                  // how to reduce the samples of a combination:
                             //   MAX_VALUE, AVERAGE_VALUE or MEDIAN_VALUE
  double opt_perf;           // best performance found so far
  int opt_combination_idx;   // combination index of the best performance
  double relative_tolerance; // accepted performance drop wrt opt_perf
  int nperf_degraded;        // # of consecutive degraded measurements

  int checked_simd;          // 1 once the final SIMD width has been checked
  bigint last_call_step;     // timestep of the last state machine update
  bigint last_beginstep;     // first timestep of the run seen last
  int window_open;           // 1 while a timing window is being collected
  bigint window_step;        // timestep at which the timing window started
  double window_cpu;         // elapsed CPU time when the window started

  FILE *tuning_logfile;      // logfile for the tuning history

 private:
  void allocate();
  void check_simd_size();                     // trim tpa values to SIMD width
  double close_window();                      // steps per second of a window
  void get_params(int, int &, int &, int &);  // parameters of a combination
  void set_param_values(int);                 // apply them to the GPU library
  int get_optimal_combination_idx();          // best performing combination
  void regular_performance_check();           // watch for degradation
};

}    // namespace LAMMPS_NS

#endif
