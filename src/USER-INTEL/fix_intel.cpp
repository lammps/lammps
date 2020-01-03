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
                        Anupama Kurpad (Intel) - Host Affinitization
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

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef _LMP_INTEL_OFFLOAD
#ifndef INTEL_OFFLOAD_NOAFFINITY
#include <unistd.h>
#endif
#endif

#include "suffix.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#ifdef __INTEL_OFFLOAD
#ifndef _LMP_INTEL_OFFLOAD
#warning "Not building Intel package with Xeon Phi offload support."
#endif
#endif

/* ---------------------------------------------------------------------- */

FixIntel::FixIntel(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal package intel command");

  int ncops = force->inumeric(FLERR,arg[3]);

  _nbor_pack_width = 1;
  _three_body_neighbor = 0;
  _pair_intel_count = 0;
  _hybrid_nonpair = 0;
  _print_pkg_info = 1;

  _precision_mode = PREC_MODE_MIXED;
  _offload_balance = -1.0;
  _overflow_flag[LMP_OVERFLOW] = 0;
  _off_overflow_flag[LMP_OVERFLOW] = 0;

  _offload_affinity_balanced = 0;
  _offload_threads = 0;
  _offload_tpc = 4;

  _force_array_s = 0;
  _force_array_m = 0;
  _force_array_d = 0;
  _ev_array_s = 0;
  _ev_array_d = 0;

  #ifdef _LMP_INTEL_OFFLOAD
  if (ncops < 0) error->all(FLERR,"Illegal package intel command");
  _offload_affinity_set = 0;
  _off_force_array_s = 0;
  _off_force_array_m = 0;
  _off_force_array_d = 0;
  _off_ev_array_s = 0;
  _off_ev_array_d = 0;
  _balance_fixed = 0.0;
  _cop = 0;
  #endif

  // optional keywords

  int nomp = 0, no_affinity = 0;
  _allow_separate_buffers = 1;
  _offload_ghost = -1;
  _lrt = 0;

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
    } else if (strcmp(arg[iarg],"no_affinity") == 0) {
      no_affinity = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "lrt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      if (strcmp(arg[iarg+1],"yes") == 0) _lrt = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) _lrt = 0;
      else error->all(FLERR,"Illegal package intel command");
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

  // if ncops is zero, just run on the cpu
  if (ncops < 1) {
    ncops = -1;
    _offload_balance = 0.0;
  }

  // if using LRT mode, create the integrate style
  if (_lrt) {
    char *str;
    str = (char *) "verlet/lrt/intel";
    update->create_integrate(1,&str,0);
  }

  // error check

  if (_offload_balance > 1.0 || _offload_threads < 0 ||
      _offload_tpc <= 0 || _offload_tpc > 4 || nomp < 0)
    error->all(FLERR,"Illegal package intel command");

  #ifdef _LMP_INTEL_OFFLOAD
  _ncops = ncops;
  if (_offload_balance != 0.0) {
    _real_space_comm = MPI_COMM_WORLD;
    if (no_affinity == 0)
      if (set_host_affinity(nomp) != 0)
        error->all(FLERR,"Could not set host affinity for offload tasks");
  }

  int max_offload_threads = 0, offload_cores = 0;
  if (_offload_balance != 0.0) {
    #pragma offload target(mic:_cop) mandatory \
      out(max_offload_threads,offload_cores)
    {
      offload_cores = omp_get_num_procs();
      omp_set_num_threads(offload_cores);
      max_offload_threads = omp_get_max_threads();
      #ifdef __AVX512F__
      if ( (offload_cores / 4) % 2 == 1) {
        offload_cores += 4;
        max_offload_threads += 4;
      }
      #endif
    }
    _max_offload_threads = max_offload_threads;
    _offload_cores = offload_cores;
    if (_offload_threads == 0) _offload_threads = offload_cores;
    if (_offload_cores > 244 && _offload_tpc > 2)
      _offload_tpc = 2;
  }
  #endif

  // set OpenMP threads
  // nomp is user setting, default = 0

  #if defined(_OPENMP)
  #if defined(__INTEL_COMPILER)
  kmp_set_blocktime(0);
  #endif
  if (nomp != 0) {
    omp_set_num_threads(nomp);
    comm->nthreads = nomp;
  } else {
    int nthreads;
    #pragma omp parallel default(none) shared(nthreads)
    nthreads = omp_get_num_threads();
    comm->nthreads = nthreads;
  }
  #endif

  // set offload params

  #ifdef _LMP_INTEL_OFFLOAD
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
    if (_offload_balance != 0.0) {
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
  mask |= PRE_REVERSE;
  mask |= MIN_PRE_REVERSE;
  #ifdef _LMP_INTEL_OFFLOAD
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  #endif
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIntel::init()
{
  #ifdef _LMP_INTEL_OFFLOAD
  output_timing_data();
  _sync_mode = 0;
  if (offload_balance() != 0.0) {
    if (offload_noghost() || force->newton_pair == 0)
      _sync_mode = 2;
    else
      _sync_mode = 1;
    if (update->whichflag == 2) _sync_mode = 1;
  }
  #endif

  const int nstyles = _pair_intel_count;
  if (force->pair_match("hybrid", 1) != NULL) {
    _pair_hybrid_flag = 1;
    if (force->newton_pair != 0 && force->pair->no_virial_fdotr_compute)
      error->all(FLERR,
                 "Intel package requires fdotr virial with newton on.");
  } else if (force->pair_match("hybrid/overlay", 1) != NULL) {
    _pair_hybrid_flag = 1;
    if (force->newton_pair != 0 && force->pair->no_virial_fdotr_compute)
      error->all(FLERR,
                 "Intel package requires fdotr virial with newton on.");
  } else
    _pair_hybrid_flag = 0;

  if (nstyles > 1 && _pair_hybrid_flag) _pair_hybrid_flag = 2;
  else if (force->newton_pair == 0) _pair_hybrid_flag = 0;

  _pair_hybrid_zero = 0;
  _zero_master = 0;

  if (_pair_hybrid_flag && _hybrid_nonpair)
    if (_pair_hybrid_flag > 1 || force->newton_pair == 0)
      _pair_hybrid_zero = 1;
  _hybrid_nonpair = 0;

  _pair_intel_count = 0;

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_balance() != 0.0) {
    _pair_hybrid_zero = 0;
    if (force->newton_pair == 0) _pair_hybrid_flag = 0;
    if (nstyles > 1)
      error->all(FLERR,
        "Currently, cannot offload more than one intel style with hybrid.");
  }
  #endif

  check_neighbor_intel();

  int off_mode = 0;
  if (_offload_balance != 0.0) off_mode = 1;
  if (_precision_mode == PREC_MODE_SINGLE) {
    _single_buffers->zero_ev();
    _single_buffers->grow_ncache(off_mode, comm->nthreads);
    _single_buffers->free_list_ptrs();
  } else if (_precision_mode == PREC_MODE_MIXED) {
    _mixed_buffers->zero_ev();
    _mixed_buffers->grow_ncache(off_mode, comm->nthreads);
    _mixed_buffers->free_list_ptrs();
  } else {
    _double_buffers->zero_ev();
    _double_buffers->grow_ncache(off_mode, comm->nthreads);
    _double_buffers->free_list_ptrs();
  }

  _need_reduce = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::setup(int vflag)
{
  if (neighbor->style != Neighbor::BIN)
    error->all(FLERR,
            "Currently, neighbor style BIN must be used with Intel package.");
  if (vflag > 3)
   error->all(FLERR,
               "Cannot currently get per-atom virials with Intel package.");
  #ifdef _LMP_INTEL_OFFLOAD
  if (neighbor->exclude_setting() != 0)
    error->all(FLERR,
     "Currently, cannot use neigh_modify exclude with Intel package offload.");
  post_force(vflag);
  #endif
}

/* ---------------------------------------------------------------------- */

void FixIntel::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ---------------------------------------------------------------------- */

void FixIntel::pair_init_check(const bool cdmessage)
{
  #ifdef INTEL_VMASK
  atom->sortfreq = 1;
  #endif

  _nbor_pack_width = 1;

  #ifdef _LMP_INTEL_OFFLOAD
  if (_offload_balance != 0.0) atom->sortfreq = 1;

  _offload_noghost = 0;
  if (force->newton_pair && _offload_ghost == 0)
    _offload_noghost = 1;

  set_offload_affinity();

  if (!_timers_allocated) {
    double *time1 = off_watch_pair();
    double *time2 = off_watch_neighbor();
    int *overflow = get_off_overflow_flag();
    if (_offload_balance !=0.0) {
      #pragma offload_transfer target(mic:_cop)  \
        nocopy(time1,time2:length(1) alloc_if(1) free_if(0)) \
        in(overflow:length(5) alloc_if(1) free_if(0))
    }
    _timers_allocated = true;
  }
  #endif
  _nthreads = comm->nthreads;

  if (_offload_balance != 0.0 && comm->me == 0) {
    #ifndef __INTEL_COMPILER_BUILD_DATE
    error->warning(FLERR, "Unknown Intel Compiler Version\n");
    #else
    if (__INTEL_COMPILER_BUILD_DATE != 20131008 &&
        __INTEL_COMPILER_BUILD_DATE < 20141023)
      error->warning(FLERR, "Unsupported Intel Compiler.");
    #endif
    #if !defined(__INTEL_COMPILER)
    error->warning(FLERR, "Unsupported Intel Compiler.");
    #endif
  }

  #ifndef LMP_INTEL_NBOR_COMPAT
  if (force->pair->manybody_flag && atom->molecular) {
    int flag = 0;
    if (atom->nbonds > 0 && force->special_lj[1] == 0.0 &&
        force->special_coul[1] == 0.0) flag = 1;
    if (atom->nangles > 0 && force->special_lj[2] == 0.0 &&
        force->special_coul[2] == 0.0) flag = 1;
    if (atom->ndihedrals > 0 && force->special_lj[3] == 0.0 &&
        force->special_coul[3] == 0.0) flag = 1;
    if (flag)
      error->all(FLERR,"Add -DLMP_INTEL_NBOR_COMPAT to build for special_bond"
                 " exclusions with Intel");
  }
  #endif

  int need_tag = 0;
  if (atom->molecular) need_tag = 1;

  // Clear buffers used for pair style
  char kmode[80];
  if (_precision_mode == PREC_MODE_SINGLE) {
    strcpy(kmode, "single");
    get_single_buffers()->need_tag(need_tag);
  } else if (_precision_mode == PREC_MODE_MIXED) {
    strcpy(kmode, "mixed");
    get_mixed_buffers()->need_tag(need_tag);
  } else {
    strcpy(kmode, "double");
    get_double_buffers()->need_tag(need_tag);
  }

  _pair_intel_count++;

  #ifdef _LMP_INTEL_OFFLOAD
  set_offload_affinity();
  #endif

  if (_print_pkg_info && comm->me == 0) {
    if (screen) {
      fprintf(screen,
              "----------------------------------------------------------\n");
      if (_offload_balance != 0.0) {
        fprintf(screen,"Using Intel Coprocessor with %d threads per core, ",
                _offload_tpc);
        fprintf(screen,"%d threads per task\n",_offload_threads);
      } else {
        fprintf(screen,"Using Intel Package without Coprocessor.\n");
      }
      fprintf(screen,"Precision: %s\n",kmode);
      if (cdmessage) {
        #ifdef LMP_USE_AVXCD
        fprintf(screen,"AVX512 CD Optimizations: Enabled\n");
        #else
        fprintf(screen,"AVX512 CD Optimizations: Disabled\n");
        #endif
      }
      fprintf(screen,
              "----------------------------------------------------------\n");
    }
  }
  _print_pkg_info = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::bond_init_check()
{
  if (_offload_balance != 0.0 && atom->molecular &&
      force->newton_pair != force->newton_bond)
    error->all(FLERR,
      "USER-INTEL package requires same setting for newton bond and non-bond.");

  int intel_pair = 0;
  if (force->pair_match("/intel", 0) != NULL)
    intel_pair = 1;
  else if (force->pair_match("hybrid", 1) != NULL) {
    _hybrid_nonpair = 1;
    if (_pair_intel_count) intel_pair = 1;
  } else if (force->pair_match("hybrid/overlay", 1) != NULL) {
    _hybrid_nonpair = 1;
    if (_pair_intel_count) intel_pair = 1;
  }

  if (intel_pair == 0)
    error->all(FLERR, "Intel styles for bond/angle/dihedral/improper "
      "require intel pair style.");
}

/* ---------------------------------------------------------------------- */

void FixIntel::kspace_init_check()
{
  int intel_pair = 0;
  if (force->pair_match("/intel", 0) != NULL)
    intel_pair = 1;
  else if (force->pair_match("hybrid", 1) != NULL) {
    _hybrid_nonpair = 1;
    if (_pair_intel_count) intel_pair = 1;
  } else if (force->pair_match("hybrid/overlay", 1) != NULL) {
    _hybrid_nonpair = 1;
    if (_pair_intel_count) intel_pair = 1;
  }

  if (intel_pair == 0)
    error->all(FLERR, "Intel styles for kspace require intel pair style.");
}

/* ---------------------------------------------------------------------- */

void FixIntel::check_neighbor_intel()
{
  #ifdef _LMP_INTEL_OFFLOAD
  _full_host_list = 0;

  const int nrequest = neighbor->nrequest;
  for (int i = 0; i < nrequest; ++i) {
    if (_offload_balance != 0.0 && neighbor->requests[i]->intel == 0) {
      _full_host_list = 1;
      _offload_noghost = 0;
    }
    if (neighbor->requests[i]->skip && _offload_balance != 0.0)
      error->all(FLERR, "Cannot yet use hybrid styles with Intel offload.");

    // avoid flagging a neighbor list as both USER-INTEL and USER-OMP
    if (neighbor->requests[i]->intel)
      neighbor->requests[i]->omp = 0;
  }
  #else
  // avoid flagging a neighbor list as both USER-INTEL and USER-OMP
  const int nrequest = neighbor->nrequest;
  for (int i = 0; i < nrequest; ++i)
    if (neighbor->requests[i]->intel)
      neighbor->requests[i]->omp = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

void FixIntel::_sync_main_arrays(const int prereverse)
{
  if (!prereverse) _zero_master = 1;
  int done_this_step = prereverse;
  if (_pair_hybrid_zero == 0) done_this_step = 1;
  if (_force_array_m != 0) {
    if (_need_reduce) {
      reduce_results(&_force_array_m[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_m, _ev_array_d, _results_eatom, _results_vatom,0);
    if (done_this_step) _force_array_m = 0;
    else _ev_array_d = 0;
  } else if (_force_array_d != 0) {
    if (_need_reduce) {
      reduce_results(&_force_array_d[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_d, _ev_array_d, _results_eatom, _results_vatom,0);
    if (done_this_step) _force_array_d = 0;
    else _ev_array_d = 0;
  } else if (_force_array_s != 0) {
    if (_need_reduce) {
      reduce_results(&_force_array_s[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_s, _ev_array_s, _results_eatom, _results_vatom,0);
    if (done_this_step) _force_array_s = 0;
    else _ev_array_s = 0;
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_sync_mode == 1) sync_coprocessor();
  #endif
}

/* ---------------------------------------------------------------------- */

void FixIntel::pre_reverse(int /*eflag*/, int /*vflag*/)
{
  _sync_main_arrays(1);
}

/* ---------------------------------------------------------------------- */

template <class acc_t>
void FixIntel::reduce_results(acc_t * _noalias const f_scalar)
{
  int o_range, f_stride;
  if (force->newton_pair)
    o_range = atom->nlocal + atom->nghost;
  else
    o_range = atom->nlocal;
  IP_PRE_get_stride(f_stride, o_range, (sizeof(acc_t)*4), lmp->atom->torque);

  o_range *= 4;
  const int f_stride4 = f_stride * 4;

  if (_nthreads <= INTEL_HTHREADS) {
    acc_t *f_scalar2 = f_scalar + f_stride4;
    if (_nthreads == 4) {
      acc_t *f_scalar3 = f_scalar2 + f_stride4;
      acc_t *f_scalar4 = f_scalar3 + f_stride4;
      _use_simd_pragma("vector aligned")
      _use_simd_pragma("simd")
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n] + f_scalar3[n] + f_scalar4[n];
    } else if (_nthreads == 2) {
      _use_simd_pragma("vector aligned")
      _use_simd_pragma("simd")
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n];
    } else {
      acc_t *f_scalar3 = f_scalar2 + f_stride4;
      _use_simd_pragma("vector aligned")
      _use_simd_pragma("simd")
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n] + f_scalar3[n];
    }
  } else {
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id_align(iifrom, iito, tid, o_range, _nthreads,
                                sizeof(acc_t));

      acc_t *f_scalar2 = f_scalar + f_stride4;
      for (int t = 1; t < _nthreads; t++) {
        _use_simd_pragma("vector aligned")
        _use_simd_pragma("simd")
        for (int n = iifrom; n < iito; n++)
          f_scalar[n] += f_scalar2[n];
        f_scalar2 += f_stride4;
      }
    }
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

template <class ft, class acc_t>
void FixIntel::add_results(const ft * _noalias const f_in,
                           const acc_t * _noalias const ev_global,
                           const int eatom, const int vatom,
                           const int offload) {
  start_watch(TIME_PACK);
  int f_length;
  #ifdef _LMP_INTEL_OFFLOAD
  if (_separate_buffers) {
    if (offload) {
      if (force->newton_pair) {
        add_oresults(f_in, ev_global, eatom, vatom, 0, _offload_nlocal);
        const acc_t * _noalias const enull = 0;
        int offset = _offload_nlocal;
        if (atom->torque) offset *= 2;
        add_oresults(f_in + offset, enull, eatom, vatom,
                     _offload_min_ghost, _offload_nghost);
      } else
        add_oresults(f_in, ev_global, eatom, vatom, 0, offload_end_pair());
    } else {
      if (force->newton_pair) {
        add_oresults(f_in, ev_global, eatom, vatom,
                     _host_min_local, _host_used_local);
        const acc_t * _noalias const enull = 0;
        int offset = _host_used_local;
        if (atom->torque) offset *= 2;
        add_oresults(f_in + offset, enull, eatom,
                     vatom, _host_min_ghost, _host_used_ghost);
      } else {
        int start = host_start_pair();
        add_oresults(f_in, ev_global, eatom, vatom, start, atom->nlocal-start);
      }
    }
    stop_watch(TIME_PACK);
    return;
  }
  int start;
  if (offload) {
    start = 0;
    if (force->newton_pair) {
      if (_offload_noghost == 0)
        f_length = atom->nlocal + atom->nghost;
      else
        f_length = atom->nlocal;
    } else
      f_length = offload_end_pair();
  } else {
    if (force->newton_pair) {
      start = 0;
      f_length = atom->nlocal + atom->nghost;
    } else {
      start = host_start_pair();
      f_length = atom->nlocal - start;
    }
  }
  add_oresults(f_in, ev_global, eatom, vatom, start, f_length);
  #else
  if (force->newton_pair)
    f_length = atom->nlocal + atom->nghost;
  else
    f_length = atom->nlocal;
  add_oresults(f_in, ev_global, eatom, vatom, 0, f_length);
  #endif
  stop_watch(TIME_PACK);
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_oresults(const ft * _noalias const f_in,
                            const acc_t * _noalias const ev_global,
                            const int eatom, const int /*vatom*/,
                            const int out_offset, const int nall) {
  lmp_ft * _noalias const f = (lmp_ft *) lmp->atom->f[0] + out_offset;
  if (atom->torque) {
    if (f_in[1].w)
    {
      if (f_in[1].w == 1)
        error->all(FLERR,"Bad matrix inversion in mldivide3");
      else
        error->all(FLERR,
                   "Sphere particles not yet supported for gayberne/intel");
    }
  }

  int packthreads;
  if (_nthreads > INTEL_HTHREADS) packthreads = _nthreads;
  else packthreads = 1;
  #if defined(_OPENMP)
  #pragma omp parallel if(packthreads > 1)
  #endif
  {
    #if defined(_OPENMP)
    const int tid = omp_get_thread_num();
    #else
    const int tid = 0;
    #endif
    int ifrom, ito;
    IP_PRE_omp_range_align(ifrom, ito, tid, nall, packthreads, sizeof(acc_t));
    if (atom->torque) {
      int ii = ifrom * 2;
      lmp_ft * _noalias const tor = (lmp_ft *) lmp->atom->torque[0] +
        out_offset;
      if (eatom) {
        double * _noalias const lmp_eatom = force->pair->eatom + out_offset;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[ii].x;
          f[i].y += f_in[ii].y;
          f[i].z += f_in[ii].z;
          lmp_eatom[i] += f_in[ii].w;
          tor[i].x += f_in[ii+1].x;
          tor[i].y += f_in[ii+1].y;
          tor[i].z += f_in[ii+1].z;
          ii += 2;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[ii].x;
          f[i].y += f_in[ii].y;
          f[i].z += f_in[ii].z;
          tor[i].x += f_in[ii+1].x;
          tor[i].y += f_in[ii+1].y;
          tor[i].z += f_in[ii+1].z;
          ii += 2;
        }
      }
    } else {
      if (eatom) {
        double * _noalias const lmp_eatom = force->pair->eatom + out_offset;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[i].x;
          f[i].y += f_in[i].y;
          f[i].z += f_in[i].z;
          lmp_eatom[i] += f_in[i].w;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[i].x;
          f[i].y += f_in[i].y;
          f[i].z += f_in[i].z;
        }
      }
    }
  }

  if (ev_global != NULL) {
    force->pair->eng_vdwl += ev_global[0];
    force->pair->eng_coul += ev_global[1];
    force->pair->virial[0] += ev_global[2];
    force->pair->virial[1] += ev_global[3];
    force->pair->virial[2] += ev_global[4];
    force->pair->virial[3] += ev_global[5];
    force->pair->virial[4] += ev_global[6];
    force->pair->virial[5] += ev_global[7];
  }
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

/* ---------------------------------------------------------------------- */

void FixIntel::post_force(int vflag)
{
  if (_sync_mode == 2) sync_coprocessor();
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_off_results(const ft * _noalias const f_in,
                               const acc_t * _noalias const ev_global) {
  if (_offload_balance < 0.0)
    _balance_other_time = MPI_Wtime() - _balance_other_time;

  start_watch(TIME_OFFLOAD_WAIT);
  #ifdef _LMP_INTEL_OFFLOAD
  if (neighbor->ago == 0) {
    #pragma offload_wait target(mic:_cop) wait(atom->tag,f_in)
  } else {
    #pragma offload_wait target(mic:_cop) wait(f_in)
  }
  #endif
  double wait_time = stop_watch(TIME_OFFLOAD_WAIT);

  int nlocal = atom->nlocal;
  if (neighbor->ago == 0) {
    if (_off_overflow_flag[LMP_OVERFLOW])
      error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    _offload_nlocal = _off_overflow_flag[LMP_LOCAL_MAX] + 1;
    _offload_min_ghost = _off_overflow_flag[LMP_GHOST_MIN];
    _offload_nghost = _off_overflow_flag[LMP_GHOST_MAX] + 1 -
      _offload_min_ghost;
    if (_offload_nghost < 0) _offload_nghost = 0;
    _offload_nall = _offload_nlocal + _offload_nghost;
      _offload_nlocal;
  }

  if (atom->torque)
    if (f_in[1].w < 0.0)
      error->all(FLERR, "Bad matrix inversion in mldivide3");
  add_results(f_in, ev_global, _off_results_eatom, _off_results_vatom, 1);

  // Load balance?
  if (_offload_balance < 0.0) {
    if (neighbor->ago == 0)
      _balance_pair = _balance_neighbor;
    double mic_time;
    mic_time = *_stopwatch_offload_pair;
    if (_balance_pair_time + _balance_other_time < mic_time) {
      double ft = _balance_pair_time + _balance_other_time + wait_time -
          mic_time;
      _balance_fixed = (1.0 - INTEL_LB_MEAN_WEIGHT) * _balance_fixed +
          INTEL_LB_MEAN_WEIGHT * ft;
    }

    double ctps = _balance_pair_time / (1.0-_balance_pair);
    double otps = mic_time / _balance_pair;
    double new_balance = (ctps + _balance_other_time - _balance_fixed) /
        (otps + ctps);
    _balance_neighbor = (1.0 - INTEL_LB_MEAN_WEIGHT) *_balance_neighbor +
        INTEL_LB_MEAN_WEIGHT * new_balance;
  }

  #ifdef TIME_BALANCE
  start_watch(TIME_IMBALANCE);
  MPI_Barrier(_real_space_comm);
  stop_watch(TIME_IMBALANCE);
  #endif
  acc_timers();
}

/* ---------------------------------------------------------------------- */

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
        double ht = timers[TIME_HOST_NEIGHBOR] + timers[TIME_HOST_PAIR] +
          timers[TIME_OFFLOAD_WAIT];
        double ct = timers[TIME_OFFLOAD_NEIGHBOR] +
          timers[TIME_OFFLOAD_PAIR];
        double tt = MAX(ht,ct);
        if (timers[TIME_OFFLOAD_LATENCY] / tt > 0.07 && _separate_coi == 0)
          error->warning(FLERR,
                 "Leaving a core free can improve performance for offload");
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
  int offload_cores = _offload_cores;
  #pragma offload target(mic:_cop) mandatory \
    in(node_rank,offload_threads,offload_tpc,offload_affinity_balanced, \
       offload_cores)
  {
    omp_set_num_threads(offload_threads);
    #pragma omp parallel
    {
      int tnum = omp_get_thread_num();
      kmp_affinity_mask_t mask;
      kmp_create_affinity_mask(&mask);
      int proc = offload_threads * node_rank + tnum;
      #ifdef __AVX512F__
      proc = (proc / offload_tpc) + (proc % offload_tpc) *
             ((offload_cores) / 4);
      proc += 68;
      #else
      if (offload_affinity_balanced)
        proc = proc * 4 - (proc / 60) * 240 + proc / 60 + 1;
      else
        proc += (proc / 4) * (4 - offload_tpc) + 1;
      #endif
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

/* ---------------------------------------------------------------------- */

int FixIntel::set_host_affinity(const int nomp)
{
  #ifndef INTEL_OFFLOAD_NOAFFINITY
  _separate_coi = 1;
  int rank = comm->me;
  int node_rank;
  int ppn = get_ppn(node_rank);
  int cop = node_rank / (ppn / _ncops);

  // Get a sorted list of logical cores
  int proc_list[INTEL_MAX_HOST_CORE_COUNT];
  int ncores;
  FILE *p;
  char cmd[512];
  char readbuf[INTEL_MAX_HOST_CORE_COUNT*5];
  sprintf(cmd, "lscpu -p | grep -v '#' |"
          "sort -t, -k 3,3n -k 2,2n | awk -F, '{print $1}'");
  p = popen(cmd, "r");
  if (p == NULL) return -1;
  ncores = 0;
  while(fgets(readbuf, 512, p)) {
    proc_list[ncores] = atoi(readbuf);
    ncores++;
  }
  pclose(p);

  // Sanity checks for core list
  if (ncores < 2) return -1;
  int nzero = 0;
  for (int i = 0; i < ncores; i++) {
    if (proc_list[i] == 0) nzero++;
    if (proc_list[i] < 0 || proc_list[i] >= ncores) return -1;
  }
  if (nzero > 1) return -1;

  // Determine the OpenMP/MPI configuration
  char *estring;
  int nthreads = nomp;
  if (nthreads == 0) {
    estring = getenv("OMP_NUM_THREADS");
    if (estring != NULL) {
      nthreads = atoi(estring);
      if (nthreads < 2) nthreads = 1;
    } else
      nthreads = 1;
  }

  // Determine how many logical cores for COI and MPI tasks
  int coi_cores = 0, mpi_cores;
  int subscription = nthreads * ppn;
  if (subscription > ncores) {
    if (rank == 0)
      error->warning(FLERR,
                     "More MPI tasks/OpenMP threads than available cores");
    return 0;
  }
  if (subscription == ncores)
    _separate_coi = 0;

  if (subscription > ncores / 2) {
    coi_cores = ncores - subscription;
    if (coi_cores > INTEL_MAX_COI_CORES) coi_cores = INTEL_MAX_COI_CORES;
  }
  mpi_cores = (ncores - coi_cores) / ppn;

  // Get ids of all LWPs that COI spawned and affinitize
  int lwp = 0, plwp = 0, nlwp = 0, mlwp = 0, fail = 0;
  cpu_set_t cpuset;
  pid_t pid = getpid();
  if (coi_cores) {
    sprintf(cmd, "ps -Lp %d -o lwp | awk ' (NR > 2) {print}'", pid);
    p = popen(cmd, "r");
    if (p == NULL) return -1;

    while(fgets(readbuf, 512, p)) {
      lwp = atoi(readbuf);
      int first = coi_cores + node_rank * mpi_cores;
      CPU_ZERO(&cpuset);
      for (int i = first; i < first + mpi_cores; i++)
        CPU_SET(proc_list[i], &cpuset);
      if (sched_setaffinity(lwp, sizeof(cpu_set_t), &cpuset)) {
        fail = 1;
        break;
      }
      plwp++;
    }
    pclose(p);

    // Do async offload to create COI threads
    int sig1, sig2;
    float *buf1;
    int pragma_size = 1024;
    buf1 = (float*) malloc(sizeof(float)*pragma_size);

    #pragma offload target (mic:0) mandatory \
      in(buf1:length(pragma_size) alloc_if(1) free_if(0))       \
      signal(&sig1)
    { buf1[0] = 0.0; }
    #pragma offload_wait target(mic:0) wait(&sig1)

    #pragma offload target (mic:0) mandatory \
      out(buf1:length(pragma_size) alloc_if(0) free_if(1))      \
      signal(&sig2)
    { buf1[0] = 1.0; }
    #pragma offload_wait target(mic:0) wait(&sig2)
    free(buf1);

    p = popen(cmd, "r");
    if (p == NULL) return -1;

    while(fgets(readbuf, 512, p)) {
      lwp = atoi(readbuf);
      nlwp++;
      if (nlwp <= plwp) continue;

      CPU_ZERO(&cpuset);
      for(int i=0; i<coi_cores; i++)
        CPU_SET(proc_list[i], &cpuset);

      if (sched_setaffinity(lwp, sizeof(cpu_set_t), &cpuset)) {
        fail = 1;
        break;
      }
    }
    pclose(p);
    nlwp -= plwp;

    // Get stats on the number of LWPs per process
    MPI_Reduce(&nlwp, &mlwp, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  if (screen && rank == 0) {
    if (coi_cores)
      fprintf(screen,"Intel Package: Affinitizing %d Offload Threads to %d Cores\n",
              mlwp, coi_cores);
    fprintf(screen,"Intel Package: Affinitizing MPI Tasks to %d Cores Each\n",mpi_cores);
  }
  if (fail) return -1;

  // Affinitize MPI Ranks
  CPU_ZERO(&cpuset);
  int first = coi_cores + node_rank * mpi_cores;
  for (int i = first; i < first+mpi_cores; i++)
    CPU_SET(proc_list[i], &cpuset);
  if (sched_setaffinity(pid, sizeof(cpu_set_t), &cpuset))
    return -1;

  #endif
  return 0;
}

#endif
