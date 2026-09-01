// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "tune_gpu.h"

#include "comm.h"
#include "error.h"
#include "gpu_extra.h"
#include "lammps_gpu.h"
#include "memory.h"
#include "neighbor.h"
#include "timer.h"
#include "update.h"

#include <algorithm>

#if (LAL_USE_OMP == 1)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;

// how the samples of a parameter combination are reduced to one value

enum { MAX_VALUE=0, AVERAGE_VALUE=1, MEDIAN_VALUE=2 };

/* ---------------------------------------------------------------------- */

TuneGPU::TuneGPU(LAMMPS *lmp, int nevery, int _nsamples, int _mode,
                 double _rel_tol, int max_nthreads) : Pointers(lmp),
  interval(nevery), performance(nullptr), tuning_logfile(nullptr)
{
  nsamples = _nsamples;
  mode = _mode;
  relative_tolerance = _rel_tol;

  ncombinations = 0;
  combination_idx = 0;
  sample_idx = 0;
  scanning_completed = 0;
  opt_perf = 0.0;
  opt_combination_idx = 0;
  nperf_degraded = 0;

  checked_simd = 0;
  last_call_step = -1;
  last_beginstep = -1;
  window_open = 0;
  window_step = 0;
  window_cpu = 0.0;

  if (interval <= 0)
    error->all(FLERR,"Interval for GPU auto-tuning must be positive");
  if (nsamples <= 0)
    error->all(FLERR,"Number of samples for GPU auto-tuning must be positive");
  if (relative_tolerance <= 0.0)
    error->all(FLERR,"Relative tolerance for GPU auto-tuning must be positive");

  // candidate values for threads per atom: powers of two up to the SIMD width

  const int max_tpa = lmp_gpu_simd_size();
  for (int tpa = 1; tpa <= max_tpa; tpa *= 2) tpa_values.push_back(tpa);
  if (tpa_values.empty()) tpa_values.push_back(1);

  // candidate values for the host threads used by the GPU library:
  // powers of two up to the thread count requested with the "omp" keyword,
  // or up to what OpenMP offers this MPI rank when "omp" was not used

  #if (LAL_USE_OMP == 1)
  if (max_nthreads <= 0) max_nthreads = omp_get_max_threads();
  #else
  max_nthreads = 1;
  #endif
  if (max_nthreads < 1) max_nthreads = 1;
  for (int nt = 1; nt <= max_nthreads; nt *= 2) nthreads_values.push_back(nt);
  if (nthreads_values.back() != max_nthreads)
    nthreads_values.push_back(max_nthreads);

  if (comm->me == 0) {
    tuning_logfile = fopen("tuning-gpu.log","w");
    if (tuning_logfile == nullptr)
      error->one(FLERR,"Cannot open tuning logfile tuning-gpu.log: {}",
                 utils::getsyserror());
  }

  allocate();
  set_param_values(combination_idx);
}

/* ---------------------------------------------------------------------- */

TuneGPU::~TuneGPU()
{
  memory->destroy(performance);

  if (tuning_logfile) {
    fclose(tuning_logfile);
    tuning_logfile = nullptr;
  }
}

/* ----------------------------------------------------------------------
   allocate the performance array for the current parameter value sets
------------------------------------------------------------------------- */

void TuneGPU::allocate()
{
  ncombinations = tpa_values.size() * nthreads_values.size();

  // rows = parameter combinations, ordered with threads per atom varying
  //   fastest, cols = samples collected for that combination

  memory->destroy(performance);
  memory->create(performance, ncombinations, nsamples, "tune_gpu:performance");

  for (int i = 0; i < ncombinations; i++)
    for (int s = 0; s < nsamples; s++) performance[i][s] = 0.0;

  scanning_completed = 0;
  combination_idx = 0;
  sample_idx = 0;
  window_open = 0;
}

/* ----------------------------------------------------------------------
   the SIMD width reported by the device may be lowered when the pair kernels
   are compiled, so drop any threads per atom value that became invalid
   once all pair styles have been initialized
------------------------------------------------------------------------- */

void TuneGPU::check_simd_size()
{
  checked_simd = 1;

  const int max_tpa = lmp_gpu_simd_size();
  if (tpa_values.back() <= max_tpa) return;

  std::vector<int> values;
  for (auto tpa : tpa_values)
    if (tpa <= max_tpa) values.push_back(tpa);
  if (values.empty()) values.push_back(1);

  tpa_values = values;
  allocate();
  set_param_values(combination_idx);
}

/* ----------------------------------------------------------------------
   scan the parameter combinations and settle on the fastest one

   called by fix GPU on every timestep when auto-tuning is enabled

   a change of the threads per atom setting only takes effect at the next
   neighbor list rebuild, because the layout of the packed neighbor list on
   the device depends on it.  timing windows therefore start and end on a
   rebuild step, so that a window only covers steps that ran with the
   parameter set being timed and always spans whole rebuild cycles.
------------------------------------------------------------------------- */

void TuneGPU::tuning_kernel_params()
{
  // each run restarts the timer, so a partially collected timing window
  // must not be carried over into the next run

  if (update->beginstep != last_beginstep) {
    last_beginstep = update->beginstep;
    last_call_step = -1;
    window_open = 0;
    if (tuning_logfile) {
      utils::print(tuning_logfile,"A new run starts...\n");
      fflush(tuning_logfile);
    }
  }

  // with run style respa this is called once per level and loop iteration,
  // but the state machine must only advance once per timestep

  if (update->ntimestep == last_call_step) return;
  last_call_step = update->ntimestep;

  if (!checked_simd) check_simd_size();

  // no timing information is available yet on the first step of a run

  if (update->ntimestep == update->beginstep) return;

  if (!scanning_completed) {

    // wait for the neighbor list rebuild that puts the pending parameters
    // in effect before opening the timing window

    if (!window_open) {
      if (neighbor->ago == 0) {
        window_step = update->ntimestep;
        window_cpu = timer->elapsed(Timer::TOTAL);
        window_open = 1;
      }
      return;
    }

    if (update->ntimestep - window_step < interval) return;
    if (neighbor->ago != 0) return;

    double perf = close_window();

    // discard a window that produced no usable timing and time it again

    if (perf <= 0.0) return;

    performance[combination_idx][sample_idx] = perf;

    if (tuning_logfile) {
      int tpa, nthreads;
      get_params(combination_idx, tpa, nthreads);
      utils::print(tuning_logfile,"t = {}: combination_idx {} sample {}: "
                   "tpa = {} omp = {} perf = {:.1f} TPS\n", update->ntimestep,
                   combination_idx, sample_idx, tpa, nthreads, perf);
      fflush(tuning_logfile);
    }

    // move on to the next parameter set, and to the next sample once all
    // combinations have been timed

    combination_idx++;
    if (combination_idx == ncombinations) {
      combination_idx = 0;
      sample_idx++;
      if (sample_idx == nsamples) scanning_completed = 1;
    }

    if (!scanning_completed) {
      set_param_values(combination_idx);
      return;
    }

    // scanning is done: switch to the best parameter set found

    opt_combination_idx = get_optimal_combination_idx();
    set_param_values(opt_combination_idx);

    if (tuning_logfile) {
      int tpa, nthreads;
      get_params(opt_combination_idx, tpa, nthreads);
      utils::print(tuning_logfile,"Finished tuning at t = {}. Found the "
                   "optimal params: tpa = {} omp = {} perf = {:.1f} TPS\n",
                   update->ntimestep, tpa, nthreads, opt_perf);
      fflush(tuning_logfile);
    }

    // be ready for another scan if the performance degrades later on

    combination_idx = 0;
    sample_idx = 0;
    nperf_degraded = 0;
    return;
  }

  regular_performance_check();
}

/* ----------------------------------------------------------------------
   close the current timing window and return its performance in timesteps
   per second, reduced over all MPI ranks so that every rank selects the
   same parameter combination
------------------------------------------------------------------------- */

double TuneGPU::close_window()
{
  double cpu_diff = timer->elapsed(Timer::TOTAL) - window_cpu;
  bigint step_diff = update->ntimestep - window_step;

  window_open = 0;

  // the slowest rank sets the pace of the simulation

  double max_cpu_diff;
  MPI_Allreduce(&cpu_diff,&max_cpu_diff,1,MPI_DOUBLE,MPI_MAX,world);

  if (max_cpu_diff <= 0.0) return 0.0;
  return static_cast<double>(step_diff) / max_cpu_diff;
}

/* ----------------------------------------------------------------------
   parameter values of a combination index
------------------------------------------------------------------------- */

void TuneGPU::get_params(int cidx, int &tpa, int &nthreads)
{
  const int ntpa = tpa_values.size();
  tpa = tpa_values[cidx % ntpa];
  nthreads = nthreads_values[cidx / ntpa];
}

/* ----------------------------------------------------------------------
   apply the parameter values of a combination index

   the number of host threads takes effect right away, the threads per atom
   value is picked up by the pair styles at the next neighbor list rebuild
------------------------------------------------------------------------- */

void TuneGPU::set_param_values(int cidx)
{
  int tpa, nthreads;
  get_params(cidx, tpa, nthreads);

  lmp_gpu_set_threads_per_atom(tpa);

  #if (LAL_USE_OMP == 1)
  omp_set_num_threads(nthreads);
  #else
  (void) nthreads;
  #endif
}

/* ----------------------------------------------------------------------
   find the best performing parameter combination
------------------------------------------------------------------------- */

int TuneGPU::get_optimal_combination_idx()
{
  if (performance == nullptr || ncombinations == 0 || nsamples == 0)
    error->all(FLERR,"No performance data available for GPU kernel tuning");

  opt_perf = 0.0;
  int opt_idx = 0;

  for (int i = 0; i < ncombinations; i++) {

    double p_i = 0.0;
    if (mode == MAX_VALUE) {
      p_i = performance[i][0];
      for (int s = 1; s < nsamples; s++)
        if (performance[i][s] > p_i) p_i = performance[i][s];

    } else if (mode == AVERAGE_VALUE) {
      double ave = 0.0;
      for (int s = 0; s < nsamples; s++) ave += performance[i][s];
      p_i = ave / nsamples;

    } else if (mode == MEDIAN_VALUE) {
      std::sort(performance[i], performance[i] + nsamples);
      if (nsamples % 2 != 0) p_i = performance[i][nsamples / 2];
      else p_i = 0.5 * (performance[i][nsamples / 2 - 1] +
                        performance[i][nsamples / 2]);
    }

    if (p_i > opt_perf) {
      opt_perf = p_i;
      opt_idx = i;
    }
  }

  return opt_idx;
}

/* ----------------------------------------------------------------------
   watch the performance after scanning has completed and start a new scan
   when it drops below the accepted threshold for long enough
------------------------------------------------------------------------- */

void TuneGPU::regular_performance_check()
{
  if (!window_open) {
    if (neighbor->ago == 0) {
      window_step = update->ntimestep;
      window_cpu = timer->elapsed(Timer::TOTAL);
      window_open = 1;
    }
    return;
  }

  if (update->ntimestep - window_step < interval) return;
  if (neighbor->ago != 0) return;

  double perf = close_window();
  if (perf <= 0.0) return;

  if (tuning_logfile) {
    int tpa, nthreads;
    get_params(opt_combination_idx, tpa, nthreads);
    utils::print(tuning_logfile,"Using the optimal params at timestep {}: "
                 "tpa = {} omp = {} current perf = {:.1f} TPS\n",
                 update->ntimestep, tpa, nthreads, perf);
    fflush(tuning_logfile);
  }

  double diff = 0.0;
  if (opt_perf > 0.0) diff = (opt_perf - perf) / opt_perf;
  if (diff <= relative_tolerance) return;

  // tolerate a few degraded windows before paying for a new scan

  if (nperf_degraded < nsamples) {
    nperf_degraded++;
    if (tuning_logfile) {
      utils::print(tuning_logfile,"t = {}: Performance degraded by {:.2f} "
                   "percent. opt perf = {:.1f} current perf = {:.1f}. "
                   "Continue collecting samples for the current parameter "
                   "set.\n", update->ntimestep, diff * 100.0, opt_perf, perf);
      fflush(tuning_logfile);
    }
    return;
  }

  scanning_completed = 0;
  combination_idx = 0;
  sample_idx = 0;
  nperf_degraded = 0;
  set_param_values(combination_idx);

  if (tuning_logfile) {
    utils::print(tuning_logfile,"t = {}: Performance degraded by {:.2f} "
                 "percent. opt perf = {:.1f} current perf = {:.1f}\n"
                 "Triggering a re-scan (disabled by setting the relative "
                 "tolerance to 1.0)..\n", update->ntimestep, diff * 100.0,
                 opt_perf, perf);
    fflush(tuning_logfile);
  }
}
