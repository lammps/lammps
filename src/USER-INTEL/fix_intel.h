/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(INTEL,FixIntel)

#else

#ifndef LMP_FIX_INTEL_H
#define LMP_FIX_INTEL_H

#include "fix.h"
#include "intel_buffers.h"
#include "force.h"
#include "pair.h"
#include "error.h"
#include "update.h"

namespace LAMMPS_NS {

class IntelData;
template <class flt_t, class acc_t> class IntelBuffers;
 
class FixIntel : public Fix {
 public:
  FixIntel(class LAMMPS *, int, char **);
  virtual ~FixIntel();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  void pair_init_check();

  // Get all forces, calculation results from coprocesser
  void sync_coprocessor();

  double memory_usage();

  typedef struct { double x,y,z; } lmp_ft;

  enum {PREC_MODE_SINGLE, PREC_MODE_MIXED, PREC_MODE_DOUBLE};
  
  inline int precision() { return _precision_mode; }
  inline IntelBuffers<float,float> * get_single_buffers() 
    { return _single_buffers; }
  inline IntelBuffers<float,double> * get_mixed_buffers() 
    { return _mixed_buffers; }
  inline IntelBuffers<double,double> * get_double_buffers() 
    { return _double_buffers; }

 protected:
  IntelBuffers<float,float> *_single_buffers;
  IntelBuffers<float,double> *_mixed_buffers;
  IntelBuffers<double,double> *_double_buffers;

  int _precision_mode, _nthreads;

 public:
  inline int* get_overflow_flag() { return _overflow_flag; }
  inline int* get_off_overflow_flag() { return _off_overflow_flag; }
  inline void add_result_array(IntelBuffers<double,double>::vec3_acc_t *f_in,
                               double *ev_in, const int offload,
                               const int eatom = 0, const int vatom = 0);
  inline void add_result_array(IntelBuffers<float,double>::vec3_acc_t *f_in,
                               double *ev_in, const int offload,
                               const int eatom = 0, const int vatom = 0);
  inline void add_result_array(IntelBuffers<float,float>::vec3_acc_t *f_in,
                               float *ev_in, const int offload,
                               const int eatom = 0, const int vatom = 0);
  inline void get_buffern(const int offload, int &nlocal, int &nall, 
			  int &minlocal);

  #ifdef _LMP_INTEL_OFFLOAD
  inline int coprocessor_number() { return _cop; }
  inline int full_host_list() { return _full_host_list; }
  void set_offload_affinity();
  inline double offload_balance() { return _offload_balance; }
  inline int offload_end_neighbor() { return _balance_neighbor * atom->nlocal; }
  inline int offload_end_pair();
  inline int host_start_neighbor()
    { if (_offload_noghost) return 0; else return offload_end_neighbor(); }
  inline int host_start_pair()
    { if (_offload_noghost) return 0; else return offload_end_pair(); }
  inline int offload_nlocal() { return _offload_nlocal; }
  inline int offload_nall() { return _offload_nall; }
  inline int offload_min_ghost() { return _offload_min_ghost; }
  inline int host_min_local() { return _host_min_local; }
  inline int host_min_ghost() { return _host_min_ghost; }
  inline int host_used_local() { return _host_used_local; }
  inline int host_used_ghost() { return _host_used_ghost; }
  inline int host_nall() { return _host_nall; }
  inline int separate_buffers() { return _separate_buffers; }
  inline int offload_noghost() { return _offload_noghost; }
  inline void set_offload_noghost(const int v) 
    { if (_offload_ghost < 0) _offload_noghost = v; }
  inline void set_neighbor_host_sizes();

  inline void zero_timers()
    { memset(_timers, 0, sizeof(double) * NUM_ITIMERS); }
  inline void start_watch(const int which) { _stopwatch[which] = MPI_Wtime(); }
  inline double stop_watch(const int which);
  inline double * off_watch_pair() { return _stopwatch_offload_pair; }
  inline double * off_watch_neighbor() { return _stopwatch_offload_neighbor; }
  inline void balance_stamp();
  inline void acc_timers();
  #else
  inline int offload_end_neighbor() { return 0; }
  inline int offload_end_pair() { return 0; }
  inline int host_start_neighbor() { return 0; }
  inline int host_start_pair() { return 0; }
  inline void zero_timers() {}
  inline void start_watch(const int which) {}
  inline double stop_watch(const int which) { return 0.0; }
  double * off_watch_pair() { return NULL; }
  double * off_watch_neighbor() { return NULL; }
  inline void balance_stamp() {}
  inline void acc_timers() {}
  inline int separate_buffers() { return 0; }
  #endif

 protected:
  int _overflow_flag[5];
  _alignvar(int _off_overflow_flag[5],64);
  int _allow_separate_buffers, _offload_ghost;
  #ifdef _LMP_INTEL_OFFLOAD
  double _balance_pair_time, _balance_other_time;
  int _offload_nlocal, _offload_nall, _offload_min_ghost, _offload_nghost;
  int _host_min_local, _host_min_ghost, _host_nall;
  int _host_used_local, _host_used_ghost;
  int _separate_buffers, _offload_noghost, _sync_at_pair, _separate_coi;
  bool _setup_time_cleared, _timers_allocated;
  void output_timing_data();
  FILE *_tscreen;

  IntelBuffers<float,float>::vec3_acc_t *_off_force_array_s;
  IntelBuffers<float,double>::vec3_acc_t *_off_force_array_m;
  IntelBuffers<double,double>::vec3_acc_t *_off_force_array_d;
  float *_off_ev_array_s;
  double *_off_ev_array_d;
  int _off_results_eatom, _off_results_vatom;
  int _full_host_list, _cop, _ncops;

  int get_ppn(int &);
  int set_host_affinity(const int);
  #endif
  void check_neighbor_intel();

  double _offload_balance, _balance_neighbor, _balance_pair, _balance_fixed;
  double _timers[NUM_ITIMERS];
  double _stopwatch[NUM_ITIMERS];
  _alignvar(double _stopwatch_offload_neighbor[1],64);
  _alignvar(double _stopwatch_offload_pair[1],64);

  template <class ft, class acc_t>
  inline void add_results(const ft * _noalias const f_in,
                          const acc_t * _noalias const ev_global,
                          const int eatom, const int vatom,
			  const int offload);

  template <class ft, class acc_t>
  inline void add_oresults(const ft * _noalias const f_in,
			   const acc_t * _noalias const ev_global,
			   const int eatom, const int vatom,
			   const int out_offset, const int nall);

  int _offload_affinity_balanced, _offload_threads, _offload_tpc;
  #ifdef _LMP_INTEL_OFFLOAD
  int _max_offload_threads, _offload_cores, _offload_affinity_set;
  int _im_real_space_task;
  MPI_Comm _real_space_comm;
  template <class ft, class acc_t>
  inline void add_off_results(const ft * _noalias const f_in,
                              const acc_t * _noalias const ev_global);
  #endif
};

/* ---------------------------------------------------------------------- */

void FixIntel::get_buffern(const int offload, int &nlocal, int &nall,
			   int &minlocal) {
  #ifdef _LMP_INTEL_OFFLOAD
  if (_separate_buffers) {
    if (offload) {
      if (neighbor->ago != 0) {
	nlocal = _offload_nlocal;
	nall = _offload_nall;
      } else {
	nlocal = atom->nlocal;
	nall = nlocal + atom->nghost;
      }
      minlocal = 0;
    } else {
      nlocal = atom->nlocal;
      nall = _host_nall;
      minlocal = _host_min_local;
    }
    return;
  }
  if (_offload_noghost && offload)
    nall = atom->nlocal;
  else
  #endif
    nall = atom->nlocal + atom->nghost;
  nlocal = atom->nlocal;
  minlocal = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<double,double>::vec3_acc_t *f_in,
                                double *ev_in, const int offload,
                                const int eatom, const int vatom) {
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_d = f_in;
    _off_ev_array_d = ev_in;
    if (_sync_at_pair == 1) sync_coprocessor();
    return;
  }
  #endif
  add_results(f_in, ev_in, eatom, vatom, 0);
  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  #ifdef _LMP_INTEL_OFFLOAD
  if (_sync_at_pair) sync_coprocessor();
  #endif
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float,double>::vec3_acc_t *f_in,
                                double *ev_in, const int offload,
                                const int eatom, const int vatom) {
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_m = f_in;
    _off_ev_array_d = ev_in;
    if (_sync_at_pair == 1) sync_coprocessor();
    return;
  }
  #endif
  add_results(f_in, ev_in, eatom, vatom, 0);
  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  #ifdef _LMP_INTEL_OFFLOAD
  if (_sync_at_pair) sync_coprocessor();
  #endif
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float,float>::vec3_acc_t *f_in,
                                float *ev_in, const int offload,
                                const int eatom, const int vatom) {
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_s = f_in;
    _off_ev_array_s = ev_in;
    if (_sync_at_pair == 1) sync_coprocessor();
    return;
  }
  #endif
  add_results(f_in, ev_in, eatom, vatom, 0);
  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  #ifdef _LMP_INTEL_OFFLOAD
  if (_sync_at_pair) sync_coprocessor();
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
      add_oresults(f_in, ev_global, eatom, vatom, 0, _offload_nlocal);
      if (force->newton_pair) {
	const acc_t * _noalias const enull = 0;
	int offset = _offload_nlocal;
	if (atom->torque) offset *= 2;
	add_oresults(f_in + offset, enull, eatom, vatom, 
		     _offload_min_ghost, _offload_nghost);
      }
    } else {
      add_oresults(f_in, ev_global, eatom, vatom,
		   _host_min_local, _host_used_local);
      if (force->newton_pair) {
	const acc_t * _noalias const enull = 0;
	int offset = _host_used_local;
	if (atom->torque) offset *= 2;
	add_oresults(f_in + offset, enull, eatom, 
		     vatom, _host_min_ghost, _host_used_ghost);
      }
    }
    stop_watch(TIME_PACK);
    return;
  }
  if (force->newton_pair && (_offload_noghost == 0 || offload == 0))
    f_length = atom->nlocal + atom->nghost;
  else
    f_length = atom->nlocal;
  #else
  if (force->newton_pair)
    f_length = atom->nlocal + atom->nghost;
  else
    f_length = atom->nlocal;
  #endif

  add_oresults(f_in, ev_global, eatom, vatom, 0, f_length);
  stop_watch(TIME_PACK);
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_oresults(const ft * _noalias const f_in,
			    const acc_t * _noalias const ev_global,
			    const int eatom, const int vatom,
			    const int out_offset, const int nall) {
  lmp_ft * _noalias const f = (lmp_ft *) lmp->atom->f[0] + out_offset;
  if (atom->torque) {
    if (f_in[1].w)
      if (f_in[1].w == 1)
        error->all(FLERR,"Bad matrix inversion in mldivide3");
      else
        error->all(FLERR,
                   "Sphere particles not yet supported for gayberne/intel");
  }

  #if defined(_OPENMP)
  #pragma omp parallel default(none)
  #endif
  {
    #if defined(_OPENMP)
    const int tid = omp_get_thread_num();
    #else
    const int tid = 0;
    #endif
    int ifrom, ito;
    IP_PRE_omp_range_align(ifrom, ito, tid, nall, _nthreads, sizeof(acc_t));
    if (atom->torque) {
      int ii = ifrom * 2;
      lmp_ft * _noalias const tor = (lmp_ft *) lmp->atom->torque[0] +
	out_offset;
      if (eatom) {
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[ii].x;
          f[i].y += f_in[ii].y;
          f[i].z += f_in[ii].z;
          force->pair->eatom[i] += f_in[ii].w;
          tor[i].x += f_in[ii+1].x;
          tor[i].y += f_in[ii+1].y;
          tor[i].z += f_in[ii+1].z;
          ii += 2;
        }
      } else {
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
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[i].x;
          f[i].y += f_in[i].y;
          f[i].z += f_in[i].z;
          force->pair->eatom[i] += f_in[i].w;
        }
      } else {
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

#ifdef _LMP_INTEL_OFFLOAD

/* ---------------------------------------------------------------------- */

int FixIntel::offload_end_pair() {
  if (neighbor->ago == 0) return _balance_neighbor * atom->nlocal;
  else return _balance_pair * atom->nlocal;
}

/* ---------------------------------------------------------------------- */

double FixIntel::stop_watch(const int which) {
  double elapsed = MPI_Wtime() - _stopwatch[which];
  _timers[which] += elapsed;
  return elapsed;
}

/* ---------------------------------------------------------------------- */

void FixIntel::balance_stamp() {
  if (_offload_balance < 0.0) {
    double ct = MPI_Wtime();
    _balance_other_time = ct;
    _balance_pair_time = ct - _stopwatch[TIME_HOST_PAIR];
  }
}

/* ---------------------------------------------------------------------- */

void FixIntel::acc_timers() {
  _timers[TIME_OFFLOAD_PAIR] += *_stopwatch_offload_pair;
  if (neighbor->ago == 0) {
    _timers[TIME_OFFLOAD_NEIGHBOR] += *_stopwatch_offload_neighbor;
    if (_setup_time_cleared == false) {
      zero_timers();
      _setup_time_cleared = true;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixIntel::set_neighbor_host_sizes() {
  _host_min_local = _overflow_flag[LMP_LOCAL_MIN];
  _host_min_ghost = _overflow_flag[LMP_GHOST_MIN];
  _host_used_local = atom->nlocal - _host_min_local;
  _host_used_ghost = _overflow_flag[LMP_GHOST_MAX] + 1 - _host_min_ghost;
  if (_host_used_ghost < 0) _host_used_ghost = 0;
  _host_nall = atom->nlocal + _host_used_ghost;
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_off_results(const ft * _noalias const f_in,
                               const acc_t * _noalias const ev_global) {
  if (_offload_balance < 0.0)
    _balance_other_time = MPI_Wtime() - _balance_other_time;

  start_watch(TIME_OFFLOAD_WAIT);
  #ifdef _LMP_INTEL_OFFLOAD
  #pragma offload_wait target(mic:_cop) wait(f_in)
  #endif
  double wait_time = stop_watch(TIME_OFFLOAD_WAIT);

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
  
  int nlocal = atom->nlocal;
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
    if (new_balance < 0.01) new_balance = 0.01;
    else if (new_balance > 0.99) new_balance = 0.99;
    _balance_neighbor = (1.0 - INTEL_LB_MEAN_WEIGHT) *_balance_neighbor +
        INTEL_LB_MEAN_WEIGHT * new_balance;
  }

  #ifdef TIME_BALANCE
  start_watch(TIME_IMBALANCE);
  MPI_Barrier(_real_space_comm);
  stop_watch(TIME_IMBALANCE);
  #endif
  acc_timers();
  if (atom->torque)
    if (f_in[1].w < 0.0)
      error->all(FLERR, "Bad matrix inversion in mldivide3");
  add_results(f_in, ev_global, _off_results_eatom, _off_results_vatom, 1);
}

#endif

}

#endif
#endif

/* ERROR/WARNING messages:

E: The 'package intel' command is required for /intel styles

Self-explanatory.

W: Could not set host affinity for offload tasks

When using offload to a coprocessor, the application will try to set affinity
for host MPI tasks and OpenMP threads and will generate a warning if unable
to do so successfully. In the unsuccessful case, you might wish to set
affinity outside of the application and performance might suffer if
hyperthreading is disable on the CPU.

E: Neighbor list overflow, boost neigh_modify one

Increase the value for neigh_modify one to allow for larger allocations for
neighbor list builds. The value required can be different for the Intel
package in order to support offload to a coprocessor.

E: Bad matrix inversion in mldivide3

This error should not occur unless the matrix is badly formed.

E: Illegal package intel command

The format for the package intel command is incorrect. Please see the
documentation.

E: fix intel has to operate on group 'all'

Self explanatory.

E: Illegal package intel mode requested

The format for the package intel command is incorrect. Please see the
documentation.

E: Specified run_style does not support the Intel package.

When using offload to a coprocessor, the Intel package requires a run style
with the intel suffix.

E: Currently, neighbor style BIN must be used with Intel package.

This is the only neighbor style that has been implemented for the Intel
package.

E: Currently, cannot use neigh_modify exclude with Intel package.

This is a current restriction of the Intel package.

W: Unknown Intel Compiler Version

The compiler version used to build LAMMPS has not been tested with
offload to a coprocessor.

W: Unsupported Intel Compiler

The compiler version used to build LAMMPS is not supported when using
offload to a coprocessor. There could be performance or correctness
issues. Please use 14.0.1.106 or 15.1.133 or later.

E: Currently, cannot use more than one intel style with hybrid.

Currently, hybrid pair styles can only use the intel suffix for one of the
pair styles.

E: Cannot yet use hybrid styles with Intel package.

The hybrid pair style configuration is not yet supported by the Intel 
package. Support is limited to hybrid/overlay or a hybrid style that does 
not require a skip list.

W: Leaving a core/node free can improve performance for offload

When each CPU is fully subscribed with MPI tasks and OpenMP threads,
context switching with threads used for offload can sometimes decrease
performance. If you see this warning, try using fewer MPI tasks/OpenMP threads
per node to leave a physical CPU core free on each node.

E: MPI tasks per node must be multiple of offload_cards

For offload to multiple coprocessors on a single node, the Intel package
requires that each coprocessor is used by the same number of MPI tasks.

W: More MPI tasks/OpenMP threads than available cores

Using more MPI tasks/OpenMP threads than available cores will typically
decrease performance.

*/
