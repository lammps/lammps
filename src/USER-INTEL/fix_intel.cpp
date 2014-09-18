/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_hybrid_overlay.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "fix_intel.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "suffix.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#ifdef __INTEL_OFFLOAD
#ifndef _LMP_INTEL_OFFLOAD
#warning "Not building Intel package with Xeon Phi offload support."
#endif
#endif

enum{NSQ,BIN,MULTI};

/* ---------------------------------------------------------------------- */

FixIntel::FixIntel(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal package intel command");

  int ncops = force->inumeric(FLERR,arg[3]);

  _precision_mode = PREC_MODE_MIXED;
  _offload_balance = 1.0;
  _overflow_flag[LMP_OVERFLOW] = 0;
  _off_overflow_flag[LMP_OVERFLOW] = 0;

  _offload_affinity_balanced = 0;
  _offload_threads = 1;
  _offload_tpc = 4;

  #ifdef _LMP_INTEL_OFFLOAD
  if (ncops < 1) error->all(FLERR,"Illegal package intel command");
  _offload_affinity_set = 0;
  _off_force_array_s = 0;
  _off_force_array_m = 0;
  _off_force_array_d = 0;
  _off_ev_array_s = 0;
  _off_ev_array_d = 0;
  _balance_fixed = 0.0;

  _cop = 0;

  int max_offload_threads, offload_cores;
  #pragma offload target(mic:_cop) mandatory \
    out(max_offload_threads,offload_cores)
  {
    offload_cores = omp_get_num_procs();
    omp_set_num_threads(offload_cores);
    max_offload_threads = omp_get_max_threads();
  }
  _max_offload_threads = max_offload_threads;
  _offload_cores = offload_cores;
  _offload_threads = offload_cores;
  #endif

  // optional keywords

  int nomp = 0;
  _allow_separate_buffers = 1;
  _offload_ghost = -1;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"omp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      nomp = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      if (strcmp(arg[iarg+1],"single") == 0) 
        _precision_mode = PREC_MODE_SINGLE;
      else if (strcmp(arg[iarg+1],"mixed") == 0) 
        _precision_mode = PREC_MODE_MIXED;
      else if (strcmp(arg[iarg+1],"double") == 0) 
        _precision_mode = PREC_MODE_DOUBLE;
      else error->all(FLERR,"Illegal package intel command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"balance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _offload_balance = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ghost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      if (strcmp(arg[iarg+1],"yes") == 0) _offload_ghost = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) _offload_ghost = 0;
      else error->all(FLERR,"Illegal package intel command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "tpc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _offload_tpc = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tptask") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _offload_threads = atoi(arg[iarg+1]);
      iarg += 2;
    }

    // undocumented options

    else if (strcmp(arg[iarg],"offload_affinity_balanced") == 0) {
      _offload_affinity_balanced = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"buffers") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _allow_separate_buffers = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal package intel command");
  }

  // error check

  if (_offload_balance > 1.0 || _offload_threads <= 0 ||
      _offload_tpc <= 0 || _offload_tpc > 4)
    error->all(FLERR,"Illegal package intel command");

  // set OpenMP threads


  // NOTE: nomp is user setting, default = 0



  // set offload params

  #ifdef _LMP_INTEL_OFFLOAD
  _ncops = ncops;
  if (_offload_balance < 0.0) {
    _balance_neighbor = 0.9;
    _balance_pair = 0.9;
  } else {
    _balance_neighbor = _offload_balance;
    _balance_pair = _offload_balance;
  }

  _tscreen = screen;
  zero_timers();
  _setup_time_cleared = false;
  _timers_allocated = false;
  #else
  _offload_balance = 0.0;
  #endif

  // set precision

  if (_precision_mode == PREC_MODE_SINGLE)
    _single_buffers = new IntelBuffers<float,float>(lmp);
  else if (_precision_mode == PREC_MODE_MIXED)
    _mixed_buffers = new IntelBuffers<float,double>(lmp);
  else
    _double_buffers = new IntelBuffers<double,double>(lmp);
}

/* ---------------------------------------------------------------------- */

FixIntel::~FixIntel()
{
  #ifdef _LMP_INTEL_OFFLOAD
  output_timing_data();
  if (_timers_allocated) {
    double *time1 = off_watch_pair();
    double *time2 = off_watch_neighbor();
    int *overflow = get_off_overflow_flag();
    if (time1 != NULL && time2 != NULL && overflow != NULL) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(time1,time2,overflow:alloc_if(0) free_if(1))
    }
  }
  #endif

  if (_precision_mode == PREC_MODE_SINGLE)
    delete _single_buffers;
  else if (_precision_mode == PREC_MODE_MIXED)
    delete _mixed_buffers;
  else
    delete _double_buffers;
}

/* ---------------------------------------------------------------------- */

int FixIntel::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIntel::init()
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_offload_balance != 0.0) atom->sortfreq = 1;
  
  if (force->newton_pair == 0)
    _offload_noghost = 0;
  else if (_offload_ghost == 0)
    _offload_noghost = 1;

  set_offload_affinity();

  output_timing_data();
  if (!_timers_allocated) {
    double *time1 = off_watch_pair();
    double *time2 = off_watch_neighbor();
    int *overflow = get_off_overflow_flag();
    if (time1 != NULL && time2 != NULL && overflow != NULL) {
      #pragma offload_transfer target(mic:_cop)  \
        nocopy(time1,time2:length(1) alloc_if(1) free_if(0)) \
        in(overflow:length(5) alloc_if(1) free_if(0))
    }
    _timers_allocated = true;
  }

  char kmode[80];
  if (_precision_mode == PREC_MODE_SINGLE)
    strcpy(kmode, "single");
  else if (_precision_mode == PREC_MODE_MIXED)
    strcpy(kmode, "mixed");
  else
    strcpy(kmode, "double");

  // print summary of settings
  if (comm->me == 0) {
    if (screen) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (_offload_balance != 0.0) {
        fprintf(screen,"using offload with %d threads per core, ",_offload_tpc);
        fprintf(screen,"%d threads per task\n",_offload_threads);
      }
      #endif
    }
  }
  if (update->whichflag == 2 && _offload_balance != 0.0) {
    if (_offload_balance == 1.0 && _offload_noghost == 0)
      _sync_at_pair = 1;
    else
      _sync_at_pair = 2;
  } else {
    _sync_at_pair = 0;
    if (strstr(update->integrate_style,"intel") == 0)
      error->all(FLERR,
		 "Specified run_style does not support the Intel package.");
  }
  #endif

  if (neighbor->style != BIN)
    error->all(FLERR,
	    "Currently, neighbor style BIN must be used with Intel package.");
  if (neighbor->exclude_setting() != 0)
    error->all(FLERR,
	    "Currently, cannot use neigh_modify exclude with Intel package.");
  int nstyles = 0;
  if (force->pair_match("hybrid", 1) != NULL) {
    PairHybrid *hybrid = (PairHybrid *) force->pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strstr(hybrid->keywords[i], "/intel") == NULL)
        nstyles++;
  } else if (force->pair_match("hybrid/overlay", 1) != NULL) {
    PairHybridOverlay *hybrid = (PairHybridOverlay *) force->pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strstr(hybrid->keywords[i], "/intel") == NULL)
        nstyles++;
      else
	force->pair->no_virial_fdotr_compute = 1;
  }
  if (nstyles > 1)
    error->all(FLERR,
	       "Currently, cannot use more than one intel style with hybrid.");

  neighbor->fix_intel = (void *)this;
  _nthreads = comm->nthreads;

  check_neighbor_intel();
  if (_precision_mode == PREC_MODE_SINGLE)
    _single_buffers->zero_ev();
  else if (_precision_mode == PREC_MODE_MIXED)
    _mixed_buffers->zero_ev();
  else
    _double_buffers->zero_ev();
}

/* ---------------------------------------------------------------------- */

void FixIntel::check_neighbor_intel()
{
  #ifdef _LMP_INTEL_OFFLOAD
  _full_host_list = 0;
  #endif
  const int nrequest = neighbor->nrequest;

  for (int i = 0; i < nrequest; ++i) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_offload_balance != 0.0 && neighbor->requests[i]->intel == 0) {
      _full_host_list = 1;
      _offload_noghost = 0;
    }	
    #endif
    if (neighbor->requests[i]->skip)
      error->all(FLERR, "Cannot yet use hybrid styles with Intel package.");
  }
}

/* ---------------------------------------------------------------------- */

void FixIntel::sync_coprocessor()
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_offload_balance != 0.0) {
    if (_off_force_array_m != 0) {
      add_off_results(_off_force_array_m, _off_ev_array_d);
      _off_force_array_m = 0;
    } else if (_off_force_array_d != 0) {
      add_off_results(_off_force_array_d, _off_ev_array_d);
      _off_force_array_d = 0;
    } else if (_off_force_array_s != 0) {
      add_off_results(_off_force_array_s, _off_ev_array_s);
      _off_force_array_s = 0;
    }
  }
  #endif
}

/* ---------------------------------------------------------------------- */

double FixIntel::memory_usage()
{
  double bytes;
  if (_precision_mode == PREC_MODE_SINGLE)
    bytes = _single_buffers->memory_usage(_nthreads);
  else if (_precision_mode == PREC_MODE_MIXED)
    bytes = _mixed_buffers->memory_usage(_nthreads);
  else
    bytes = _double_buffers->memory_usage(_nthreads);

  return bytes;
}

/* ---------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD

void FixIntel::output_timing_data() {
  if (_im_real_space_task == 0 || _offload_affinity_set == 0) return;

  double timer_total = 0.0;
  int size, rank;
  double timers[NUM_ITIMERS];
  MPI_Comm_size(_real_space_comm, &size);
  MPI_Comm_rank(_real_space_comm, &rank);
  MPI_Allreduce(&_timers, &timers, NUM_ITIMERS, MPI_DOUBLE, MPI_SUM,
                _real_space_comm);
  for (int i=0; i < NUM_ITIMERS; i++) {
    timers[i] /= size;
    timer_total += timers[i];
  }
  #ifdef TIME_BALANCE
  double timers_min[NUM_ITIMERS], timers_max[NUM_ITIMERS];
  MPI_Allreduce(&_timers, &timers_max, NUM_ITIMERS, MPI_DOUBLE, MPI_MAX,
                _real_space_comm);
  MPI_Allreduce(&_timers, &timers_min, NUM_ITIMERS, MPI_DOUBLE, MPI_MIN,
                _real_space_comm);
  #endif

  if (timer_total > 0.0) {
    double balance_out[2], balance_in[2];
    balance_out[0] = _balance_pair;
    balance_out[1] = _balance_neighbor;
    MPI_Reduce(balance_out, balance_in, 2, MPI_DOUBLE, MPI_SUM,
	       0, _real_space_comm);
    balance_in[0] /= size;
    balance_in[1] /= size;

    if (rank == 0 && _tscreen) {
      fprintf(_tscreen, "\n------------------------------------------------\n");
      fprintf(_tscreen, "               Offload Timing Data\n");
      fprintf(_tscreen, "------------------------------------------------\n");
      fprintf(_tscreen, "  Data Pack/Cast Seconds    %f\n",
              timers[TIME_PACK]);
      if (_offload_balance != 0.0) {
        fprintf(_tscreen, "  Host Neighbor Seconds     %f\n",
                timers[TIME_HOST_NEIGHBOR]);
        fprintf(_tscreen, "  Host Pair Seconds         %f\n",
                timers[TIME_HOST_PAIR]);
        fprintf(_tscreen, "  Offload Neighbor Seconds  %f\n",
                timers[TIME_OFFLOAD_NEIGHBOR]);
        fprintf(_tscreen, "  Offload Pair Seconds      %f\n",
                timers[TIME_OFFLOAD_PAIR]);
        fprintf(_tscreen, "  Offload Wait Seconds      %f\n",
                timers[TIME_OFFLOAD_WAIT]);
        fprintf(_tscreen, "  Offload Latency Seconds   %f\n",
                timers[TIME_OFFLOAD_LATENCY]);
        fprintf(_tscreen, "  Offload Neighbor Balance  %f\n",
                balance_in[1]);
        fprintf(_tscreen, "  Offload Pair Balance      %f\n",
                balance_in[0]);
	fprintf(_tscreen, "  Offload Ghost Atoms       ");
	if (_offload_noghost) fprintf(_tscreen,"No\n");
	else fprintf(_tscreen,"Yes\n");
        #ifdef TIME_BALANCE
        fprintf(_tscreen, "  Offload Imbalance Seconds %f\n",
                timers[TIME_IMBALANCE]);
	fprintf(_tscreen, "  Offload Min/Max Seconds   ");
	for (int i = 0; i < NUM_ITIMERS; i++)
	  fprintf(_tscreen, "[%f, %f] ",timers_min[i],timers_max[i]);
	fprintf(_tscreen, "\n");
        #endif
      }
      fprintf(_tscreen, "------------------------------------------------\n");
    }
    zero_timers();
    _setup_time_cleared = false;
  }
}

/* ---------------------------------------------------------------------- */

int FixIntel::get_ppn(int &node_rank) {
  int nprocs;
  int rank;
  MPI_Comm_size(_real_space_comm, &nprocs);
  MPI_Comm_rank(_real_space_comm, &rank);

  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(node_name,&name_length);
  node_name[name_length] = '\0';
  char *node_names = new char[MPI_MAX_PROCESSOR_NAME*nprocs];
  MPI_Allgather(node_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, node_names,
		MPI_MAX_PROCESSOR_NAME, MPI_CHAR, _real_space_comm);
  int ppn = 0;
  node_rank = 0;
  for (int i = 0; i < nprocs; i++) {
    if (strcmp(node_name, node_names + i * MPI_MAX_PROCESSOR_NAME) == 0) {
      ppn++;
      if (i < rank)
	node_rank++;
    }
  }

  return ppn;
}

/* ---------------------------------------------------------------------- */

void FixIntel::set_offload_affinity()
{
  _separate_buffers = 0;
  if (_allow_separate_buffers)
    if (_offload_balance != 0.0 && _offload_balance < 1.0)
      _separate_buffers = 1;

  _im_real_space_task = 1;
  if (strncmp(update->integrate_style,"verlet/split",12) == 0) {
    _real_space_comm = world;
    if (universe->iworld != 0) {
      _im_real_space_task = 0;
      return;
    }
  } else
    _real_space_comm = universe->uworld;

  if (_offload_balance == 0.0) _cop = -1;
  if (_offload_balance == 0.0 || _offload_affinity_set == 1)
    return;

  _offload_affinity_set = 1;
  int node_rank;
  int ppn = get_ppn(node_rank);

  if (ppn % _ncops != 0)
    error->all(FLERR, "MPI tasks per node must be multiple of offload_cards");
  ppn = ppn / _ncops;
  _cop = node_rank / ppn;
  node_rank = node_rank % ppn;

  int max_threads_per_task = _offload_cores / 4 * _offload_tpc / ppn;
  if (_offload_threads > max_threads_per_task)
    _offload_threads = max_threads_per_task;
  if (_offload_threads > _max_offload_threads)
    _offload_threads = _max_offload_threads;

  int offload_threads = _offload_threads;
  int offload_tpc = _offload_tpc;
  int offload_affinity_balanced = _offload_affinity_balanced;
  #pragma offload target(mic:_cop) mandatory \
    in(node_rank,offload_threads,offload_tpc,offload_affinity_balanced)
  {
    omp_set_num_threads(offload_threads);
    #pragma omp parallel
    {
      int tnum = omp_get_thread_num();
      kmp_affinity_mask_t mask;
      kmp_create_affinity_mask(&mask);
      int proc;
      if (offload_affinity_balanced) {
	proc = offload_threads * node_rank + tnum;
	proc = proc * 4 - (proc / 60) * 240 + proc / 60 + 1;
      } else {
	proc = offload_threads * node_rank + tnum;
	proc += (proc / 4) * (4 - offload_tpc) + 1;
      }
      kmp_set_affinity_mask_proc(proc, &mask);
      if (kmp_set_affinity(&mask) != 0)
	printf("Could not set affinity on rank %d thread %d to %d\n",
	       node_rank, tnum, proc);
    }
  }
  if (_precision_mode == PREC_MODE_SINGLE)
    _single_buffers->set_off_params(offload_threads, _cop, _separate_buffers);
  else if (_precision_mode == PREC_MODE_MIXED)
    _mixed_buffers->set_off_params(offload_threads, _cop, _separate_buffers);
  else
    _double_buffers->set_off_params(offload_threads, _cop, _separate_buffers);
}

#endif
