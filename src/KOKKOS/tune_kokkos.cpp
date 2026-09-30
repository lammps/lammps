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

#include "tune_kokkos.h"

#include "comm.h"
#include "error.h"
#include "lammps.h"
#include "kokkos.h"
#include "memory.h"
#include "timer.h"
#include "update.h"

#include <algorithm>

using namespace LAMMPS_NS;

// determine the representative value based on the collected samples for each parameter combination
enum { MAX_VALUE=0, AVERAGE_VALUE=1, MEDIAN_VALUE=2 };

//#define TUNE_DEBUG

/* ---------------------------------------------------------------------- */

TuneKokkos::TuneKokkos(LAMMPS *lmp, int _kernel_type, int nevery,
  int _nparams, const char* _name, int _nsamples, int _mode, double _rel_tol) : Pointers(lmp),
  kernel_type(_kernel_type), interval(nevery), num_params(_nparams),
  performance(nullptr), tuning_logfile(nullptr)
{
  ncombinations = 0;
  nsamples = _nsamples;
  mode = _mode;

  allocated = 0;
  firststep = 1;
  opt_perf = 0.0;
  scanning_completed = 0;

  // trigger rescaning if performance drop below relative_tolerance
  //   compared to the optimal performance obtained from the last scanning

  relative_tolerance = _rel_tol;
  nperf_degraded = 0;

  // tuner name used for logging

  if (_name) my_name = _name;
  else my_name = "default";

  // ensure that relevant kokkos parameters are allowed to be specified

  if (kernel_type == PAIR) {
    if (!lmp->kokkos->neigh_thread_set) {
      lmp->kokkos->neigh_thread_set = 1;
      lmp->kokkos->neigh_thread = 1;
    }
    if (!lmp->kokkos->threads_per_atom_set) lmp->kokkos->threads_per_atom_set = 1;
    if (!lmp->kokkos->pair_team_size_set) lmp->kokkos->pair_team_size_set = 1;

  } else if (kernel_type == BOND) {
    if (!lmp->kokkos->bond_chunk_size_set) lmp->kokkos->bond_chunk_size_set = 1;

  } else if (kernel_type == GENERIC) {
    // no specific kokkos parameters to set for generic kernels
    // leaving it for the pair/bond/fix/compute style to decide

  } else
    error->all(FLERR,"KOKKOS tuning_kernel_params: kernel type not yet supported");

  // error check the auto-tuning parameters

  if (interval <= 0)
    error->all(FLERR,"Interval for KOKKOS auto-tuning must be positive");

  if (nsamples <= 0)
    error->all(FLERR,"Number of samples for KOKKOS auto-tuning must be positive");

  if (relative_tolerance <= 0.0)
    error->all(FLERR,"Relative tolerance for KOKKOS auto-tuning must be positive");

  if (comm->me == 0) {
    std::string filename = fmt::format("tuning-{}.log", my_name);
    tuning_logfile = fopen(filename.c_str(),"w");
    if (tuning_logfile == nullptr)
      error->all(FLERR,"Cannot open tuning logfile {}: {}",
        filename.c_str(), utils::getsyserror());
  } else {
    tuning_logfile = nullptr;
  }

  allocate(num_params);
}

/* ---------------------------------------------------------------------- */

TuneKokkos::~TuneKokkos()
{
  memory->destroy(performance);

  if (tuning_logfile) {
    fclose(tuning_logfile);
    tuning_logfile = nullptr;
  }
}

/* ----------------------------------------------------------------------
   allocate arrays to store the parameter combinations and performance data
   num_params = 1: only team size (bond compute and nbin kernels)
                2: team size and threads per atom (typical pair compute kernels)
------------------------------------------------------------------------- */

void TuneKokkos::allocate(int num_params)
{
  if (num_params < 1 || num_params > 2)
    error->all(FLERR,"Illegal number of parameters for Kokkos kernel tuning");

  if (allocated) return;

  // initialize the possible values of the parameters to scan

  team_sizes.clear();
  vector_sizes.clear();

  if (num_params == 2) {

    // team sizes are used for pair team size
    // maximum team size is often limited by the GPU architecture, 512 is a safe choice

    int max_team_size = 512;
    for (int ts = 64; ts <= max_team_size; ts += 32)
      team_sizes.push_back(ts);

    // vector sizes are used for threads per atom in pair compute kernels

    #if defined(KOKKOS_ENABLE_HIP)
    int max_vectorsize = 64;
    #else
    int max_vectorsize = 32;
    #endif

    for (int vs = 1; vs <= max_vectorsize; vs *= 2)
      vector_sizes.push_back(vs);

  } else {

    // for num_params = 1 (bond compute), team sizes are used for ChunkSize

    int max_team_size = 16;
    for (int ts = 1; ts <= max_team_size; ts *= 2)
      team_sizes.push_back(ts);

     if (lmp->kokkos->threads_per_atom_set)
      vector_sizes.push_back(lmp->kokkos->threads_per_atom);
    else
      vector_sizes.push_back(1);
  }

  // compute the total number of parameter combinations

  int num_team_sizes = team_sizes.size();
  int num_vector_sizes = vector_sizes.size();

  ncombinations = num_team_sizes * num_vector_sizes;

  // allocate the 2-d performance as a 1-d array
  //   cols = team sizes (pair/team/size or bond/chunk/size = 64, 96, ..., 512)
  //   rows = vector sizes (threads/per/atom = 1, 2, 4, 8, .., max_vectorsize)

  memory->destroy(performance);
  memory->create(performance, ncombinations, nsamples, "tune_kokkos:performance");

  scanning_completed = 0;
  combination_idx = 0;
  sample_idx = 0;
  allocated = 1;
}

/* ----------------------------------------------------------------------
   tuning the pair compute kernel parameters
   this function is called by the /kk pair style at every timestep
   if auto-tuning is enabled.

   NOTE: For pair hybrid, each kk pair style creates its own TuneKokkos object,
   which sets lmp->kokkos->pair_team_size and lmp->kokkos->threads_per_atom.
   This is possible because these parameters are effective right before
   the pair's own pair compute kernel is launched in compute().
------------------------------------------------------------------------- */

void TuneKokkos::tuning_kernel_params()
{
  // retrieve the current parameter set from combination_idx

  int current_team_size, current_vector_size;
  get_current_params(combination_idx, current_team_size, current_vector_size);

  bigint elapsed_steps = update->ntimestep - update->beginstep;
  if (elapsed_steps == 0) {
    if (tuning_logfile) {
      utils::print(tuning_logfile, "A new run starts...\n");
      fflush(tuning_logfile);
    }
  }

  if (!scanning_completed) {

    // set the KOKKOS kernel parameters

    set_param_values(combination_idx);

    // wait for interval timesteps to collect timing info

    if (elapsed_steps % interval == 0) {

      double tps = get_timing_info();

      // skip for the first call when no timing info is available

      if (tps == 0.0) return;

      // store the performance of the current parameter set

      performance[combination_idx][sample_idx] = 1.0 / tps;

      if (tuning_logfile) {
        std::string mesg = fmt::format("t = {}: combination_idx {} sample {}: team size = {} ",
                            update->ntimestep, combination_idx, sample_idx, current_team_size);
        mesg += fmt::format("vector size = {} ", current_vector_size);
        mesg += fmt::format("perf = {:.1f} TPS\n", performance[combination_idx][sample_idx]);
        utils::print(tuning_logfile, "{}", mesg.c_str());
        fflush(tuning_logfile);

        #ifdef TUNE_DEBUG
        utils::logmesg(lmp, mesg);
        #endif
      }

      // move to the next parameter set

      combination_idx++;

      if (combination_idx == ncombinations) {

        // move to the next sample
        // reset combination_idx to zero to start the next scan
        //   if there are more samples to collect

        sample_idx++;
        combination_idx = 0;

        // if all samples collected, mark scanning as completed
        //   to trigger the selection of optimal parameter set
        if (sample_idx == nsamples)
          scanning_completed = 1;
      }

    }
  }

  // if scanning has just completed, find the parameter set with the optimal performance

  if (scanning_completed && sample_idx == nsamples) {
    opt_combination_idx = get_optimal_combination_idx();
    set_param_values(opt_combination_idx);

    if (tuning_logfile) {
      int opt_ts = 0, opt_vs = 0;
      get_current_params(opt_combination_idx, opt_ts, opt_vs);
      std::string mesg = fmt::format("Finished tuning at t = {}. Found the optimal params: ", update->ntimestep);
      mesg += fmt::format("team size = {} vector size = {} ",
                      opt_ts, opt_vs);
      mesg += fmt::format(" perf = {:.1f} TPS\n", opt_perf);
      utils::print(tuning_logfile, "{}", mesg.c_str());
      fflush(tuning_logfile);
    }

    // reset combination_idx and sample_idx to zero to be ready for another scan if needed
    // and to avoid repeating this block

    combination_idx = 0;
    sample_idx = 0;
    nperf_degraded = 0;
  }

  // check if the performance is within acceptable range of the optimal performance
  // if not, re-trigger the scanning process

  regular_performance_check();
}

/* ----------------------------------------------------------------------
   set the possible values of team sizes for tuning
------------------------------------------------------------------------- */

void TuneKokkos::set_team_size_values(const std::vector<int>& tsizes)
{
  team_sizes.clear();
  team_sizes = tsizes;

  int num_team_sizes = team_sizes.size();
  int num_vector_sizes = vector_sizes.size();

  ncombinations = num_team_sizes * num_vector_sizes;

  // allocate the 2-d performance as a 1-d array
  //   cols = team sizes (pair/team/size or bond/block/size = 64, 96, ..., 512)
  //   rows = vector sizes (threads/per/atom = 1, 2, 4, 8, .., max_vectorsize)

  memory->destroy(performance);
  memory->create(performance, ncombinations, nsamples, "tune_kokkos:performance");

  scanning_completed = 0;
  combination_idx = 0;
  sample_idx = 0;
}

/* ----------------------------------------------------------------------
   set the possible values of vector sizes for tuning
------------------------------------------------------------------------- */

void TuneKokkos::set_vector_size_values(const std::vector<int>& vsizes)
{
  vector_sizes.clear();
  vector_sizes = vsizes;

  int num_team_sizes = team_sizes.size();
  int num_vector_sizes = vector_sizes.size();

  ncombinations = num_team_sizes * num_vector_sizes;

  // allocate the 2-d performance as a 1-d array
  //   cols = team sizes (pair/team/size or bond/block/size = 64, 96, ..., 512)
  //   rows = vector sizes (threads/per/atom = 1, 2, 4, 8, .., max_vectorsize)

  memory->destroy(performance);
  memory->create(performance, ncombinations, nsamples, "tune_kokkos:performance");

  scanning_completed = 0;
  combination_idx = 0;
  sample_idx = 0;
}

/* ----------------------------------------------------------------------
   figure out CPU time per timestep since last time checked
------------------------------------------------------------------------- */

double TuneKokkos::get_timing_info()
{
  double dvalue;
  double new_cpu;
  bigint new_step = update->ntimestep;

  if (firststep == 1) {
    new_cpu = 0.0;
    dvalue = 0.0;
    firststep = 0;
  } else {
    new_cpu = timer->elapsed(Timer::TOTAL);
    double cpu_diff = new_cpu - last_cpu;
    bigint step_diff = new_step - last_step;
    if (step_diff > 0.0) dvalue = cpu_diff / step_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_cpu = new_cpu;

  return dvalue;
}

/* ----------------------------------------------------------------------
   get the current team size based on the combination index
------------------------------------------------------------------------- */

int TuneKokkos::get_current_team_size()
{
  int num_team_sizes = team_sizes.size();
  int current_team_size = team_sizes[combination_idx % num_team_sizes];
  return current_team_size;
}

/* ----------------------------------------------------------------------
   get the current vector size based on the combination index
------------------------------------------------------------------------- */

int TuneKokkos::get_current_vector_size()
{
  int num_team_sizes = team_sizes.size();
  int current_vector_size = vector_sizes[combination_idx / num_team_sizes];
  return current_vector_size;
}


/* ----------------------------------------------------------------------
   get the current params based on the combination index
   NOTE: using auto& with an eye toward supporting different types of params
   in the future, not restricted to int's
------------------------------------------------------------------------- */

void TuneKokkos::get_current_params(int cidx, auto& team_size, auto& vector_size)
{
  int num_team_sizes = team_sizes.size();
  team_size = team_sizes[cidx % num_team_sizes];
  vector_size = vector_sizes[cidx / num_team_sizes];
}

/* ----------------------------------------------------------------------
   set the kernel parameters based on the combination index and kernel type
------------------------------------------------------------------------- */

void TuneKokkos::set_param_values(int cidx)
{
  int num_team_sizes = team_sizes.size();
  int current_team_size = team_sizes[cidx % num_team_sizes];
  int current_vector_size = vector_sizes[cidx / num_team_sizes];

  if (kernel_type == PAIR) {
    lmp->kokkos->pair_team_size = current_team_size;
    lmp->kokkos->threads_per_atom = current_vector_size;

  } else if (kernel_type == BOND) {
    lmp->kokkos->bond_chunk_size = current_team_size;

  } else if (kernel_type == NBIN) {

    // leave it to the caller to decide how to use the team size and/or vector size

  } else if (kernel_type == FIX) {

  } else if (kernel_type == COMPUTE) {

    // leave it to the caller to decide how to use the team size and/or vector size

  } else if (kernel_type == GENERIC) {

    // leave it to the caller to decide how to use the team size and/or vector size
  }

}

/* ----------------------------------------------------------------------
   find the optimal performance from the stored performance data
   return the index of the optimal parameter set
------------------------------------------------------------------------- */

int TuneKokkos::get_optimal_combination_idx()
{
  if (performance == nullptr || ncombinations == 0 || nsamples == 0)
    error->all(FLERR,"No performance data available for Kokkos kernel tuning");

  opt_perf = performance[0][0];
  int opt_idx = 0;

  // iterate through the combinations to find the one with the best performance (highest TPS)

  for (int i = 0; i < ncombinations; i++) {

    double p_i = 0;
    if (mode == MAX_VALUE) {
      double p = performance[i][0];
      for (int s = 1; s < nsamples; s++) {
        if (performance[i][s] > p)
          p = performance[i][s];
      }
      p_i = p;

    } else if (mode == AVERAGE_VALUE) {
      double ave = 0.0;
      for (int s = 0; s < nsamples; s++)
        ave += performance[i][s];
      p_i = ave / nsamples;

    } else if (mode == MEDIAN_VALUE) {
      // sort the array performance[i]
      std::sort(performance[i], performance[i] + nsamples);
      if (nsamples % 2 != 0)
        p_i = performance[i][nsamples / 2];
      else
        p_i = (performance[i][nsamples / 2 - 1] + performance[i][nsamples / 2]) / 2.0;
    }

    if (p_i > opt_perf) {
      opt_perf = p_i;
      opt_idx = i;
    }
  }

  return opt_idx;
}

/* ----------------------------------------------------------------------
   regularly monitor the current performance after scanning is completed
   if performance degraded beyond acceptable threshold,
     re-trigger the scanning process
   otherwise, continue using the optimal parameter set,
     including with a new run
------------------------------------------------------------------------- */

void TuneKokkos::regular_performance_check()
{
  bigint elapsed_steps = update->ntimestep - update->beginstep;
  if (elapsed_steps == 0) return;

  if (!scanning_completed || elapsed_steps % interval != 0) return;

  double tps = get_timing_info();
  if (tps <= 0.0) return;
  double perf = 1.0 / tps;

  if (tuning_logfile) {
    int opt_ts = 0, opt_vs = 0;
    int opt_idx = opt_combination_idx;
    get_current_params(opt_idx, opt_ts, opt_vs);
    std::string mesg = fmt::format("Using the optimal params at timestep {}: ",
                                    update->ntimestep);
    mesg += fmt::format("team size = {} vector size = {} ",
                        opt_ts, opt_vs);
    mesg += fmt::format(" current perf = {:.1f} TPS\n", perf);
    utils::print(tuning_logfile, "{}", mesg.c_str());
    fflush(tuning_logfile);
    #ifdef TUNE_DEBUG
    utils::logmesg(lmp, mesg);
    #endif
  }

  // compute the relative performance difference wrt the optimal performance

  double diff;
  if (opt_perf > 0.0)
    diff = (opt_perf - perf) / opt_perf;
  else
    diff = 0.0;

  // if performance degraded beyond acceptable threshold after a number of samples

  if (diff > relative_tolerance) {
    if (nperf_degraded < nsamples) {
      if (comm->me == 0) {
        std::string mesg = fmt::format("t = {}: Performance degraded by {:.2f} after {} steps, "
                                       " but still collecting samples for the current parameter set. ",
                                         update->ntimestep, diff * 100.0, elapsed_steps);
        mesg += fmt::format("opt perf = {:.1f} current perf = {:.1f} ",
                          opt_perf, perf);
        mesg += fmt::format("Continue collecting samples for the current parameter set.\n");
        utils::print(tuning_logfile, "{}", mesg.c_str());
        fflush(tuning_logfile);
        #ifdef TUNE_DEBUG
        utils::logmesg(lmp, mesg);
        #endif
      }
      nperf_degraded++;
      return;
    }

    scanning_completed = 0;
    combination_idx = 0;
    sample_idx = 0;
    nperf_degraded = 0;
    firststep = 1;

    if (comm->me == 0) {
      std::string mesg = fmt::format("t = {}: Performance degraded by {:.2f} after {} steps.",
                        update->ntimestep, diff * 100.0, elapsed_steps);
      mesg += fmt::format(" opt perf = {:.1f} current perf = {:.1f} ", opt_perf, perf);
      mesg += fmt::format("\nTriggering a re-scan (disabled by setting relative tolerance to be 1.0)..\n");
      utils::print(tuning_logfile, "{}", mesg.c_str());
      fflush(tuning_logfile);
      #ifdef TUNE_DEBUG
      utils::logmesg(lmp, mesg);
      #endif
    }

  }
}


