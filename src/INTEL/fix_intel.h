// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#ifdef FIX_CLASS
// clang-format off
FixStyle(INTEL,FixIntel);
// clang-format on
#else

#ifndef LMP_FIX_INTEL_H
#define LMP_FIX_INTEL_H

#include "error.h"
#include "fix.h"
#include "force.h"
#include "intel_buffers.h"
#include "pair.h"
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
  inline void min_setup(int in) { setup(in); }
  void setup_pre_reverse(int eflag = 0, int vflag = 0);

  bool pair_hybrid_check();
  void pair_init_check(const bool cdmessage = false);
  void bond_init_check();
  void kspace_init_check();

  void pre_reverse(int eflag = 0, int vflag = 0);
  inline void min_pre_reverse(int eflag = 0, int vflag = 0) { pre_reverse(eflag, vflag); }

  void post_run() { _print_pkg_info = 1; }

  // Get all forces, calculation results from coprocesser
  void sync_coprocessor();

  double memory_usage();

  typedef struct {
    double x, y, z;
  } lmp_ft;

  enum { PREC_MODE_SINGLE, PREC_MODE_MIXED, PREC_MODE_DOUBLE };

  inline int precision() { return _precision_mode; }
  inline IntelBuffers<float, float> *get_single_buffers() { return _single_buffers; }
  inline IntelBuffers<float, double> *get_mixed_buffers() { return _mixed_buffers; }
  inline IntelBuffers<double, double> *get_double_buffers() { return _double_buffers; }

  inline int nbor_pack_width() const { return _nbor_pack_width; }
  inline void nbor_pack_width(const int w) { _nbor_pack_width = w; }
  inline int three_body_neighbor() { return _three_body_neighbor; }
  inline void three_body_neighbor(const int i) { _three_body_neighbor = i; }

  inline int need_zero(const int tid)
  {
    if (_need_reduce == 0 && tid > 0)
      return 1;
    else if (_zero_master && tid == 0) {
      _zero_master = 0;
      return 1;
    } else
      return 0;
  }
  inline void set_reduce_flag()
  {
    if (_nthreads > 1) _need_reduce = 1;
  }
  inline int lrt()
  {
    if (force->kspace_match("^pppm/.*intel$", 0) && update->whichflag == 1)
      return _lrt;
    else
      return 0;
  }
  inline int pppm_table()
  {
    if (force->kspace_match("^pppm/.*intel$", 0))
      return INTEL_P3M_TABLE;
    else
      return 0;
  }

 protected:
  IntelBuffers<float, float> *_single_buffers;
  IntelBuffers<float, double> *_mixed_buffers;
  IntelBuffers<double, double> *_double_buffers;

  int _precision_mode, _nthreads, _nbor_pack_width, _three_body_neighbor;
  int _pair_intel_count, _pair_hybrid_flag, _print_pkg_info;
  // These should be removed in subsequent update w/ simpler hybrid arch
  int _pair_hybrid_zero, _hybrid_nonpair, _zero_master;

 public:
  inline int *get_overflow_flag() { return _overflow_flag; }
  inline int *get_off_overflow_flag() { return _off_overflow_flag; }
  inline void add_result_array(IntelBuffers<double, double>::vec3_acc_t *f_in, double *ev_in,
                               const int offload, const int eatom = 0, const int vatom = 0,
                               const int rflag = 0);
  inline void add_result_array(IntelBuffers<float, double>::vec3_acc_t *f_in, double *ev_in,
                               const int offload, const int eatom = 0, const int vatom = 0,
                               const int rflag = 0);
  inline void add_result_array(IntelBuffers<float, float>::vec3_acc_t *f_in, float *ev_in,
                               const int offload, const int eatom = 0, const int vatom = 0,
                               const int rflag = 0);
  inline void get_buffern(const int offload, int &nlocal, int &nall, int &minlocal);

#ifdef _LMP_INTEL_OFFLOAD
  void post_force(int vflag);
  inline int coprocessor_number() { return _cop; }
  inline int full_host_list() { return _full_host_list; }
  void set_offload_affinity();
  inline double offload_balance() { return _offload_balance; }
  inline int offload_end_neighbor();
  inline int offload_end_pair();
  inline int host_start_neighbor()
  {
    if (_offload_noghost)
      return 0;
    else
      return offload_end_neighbor();
  }
  inline int host_start_pair()
  {
    if (_offload_noghost)
      return 0;
    else
      return offload_end_pair();
  }
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
  {
    if (_offload_ghost < 0) _offload_noghost = v;
  }
  inline void set_neighbor_host_sizes();

  inline void zero_timers() { memset(_timers, 0, sizeof(double) * NUM_ITIMERS); }
  inline void start_watch(const int which) { _stopwatch[which] = MPI_Wtime(); }
  inline double stop_watch(const int which);
  inline double *off_watch_pair() { return _stopwatch_offload_pair; }
  inline double *off_watch_neighbor() { return _stopwatch_offload_neighbor; }
  inline void balance_stamp();
  inline void acc_timers();
#else
  inline int offload_end_neighbor() { return 0; }
  inline int offload_end_pair() { return 0; }
  inline int host_start_neighbor() { return 0; }
  inline int host_start_pair() { return 0; }
  inline void zero_timers() {}
  inline void start_watch(const int /*which*/) {}
  inline double stop_watch(const int /*which*/) { return 0.0; }
  double *off_watch_pair() { return nullptr; }
  double *off_watch_neighbor() { return nullptr; }
  inline void balance_stamp() {}
  inline void acc_timers() {}
  inline int separate_buffers() { return 0; }
#endif

 protected:
  int _overflow_flag[5];
  _alignvar(int _off_overflow_flag[5], 64);
  int _allow_separate_buffers, _offload_ghost, _lrt;

  IntelBuffers<float, float>::vec3_acc_t *_force_array_s;
  IntelBuffers<float, double>::vec3_acc_t *_force_array_m;
  IntelBuffers<double, double>::vec3_acc_t *_force_array_d;
  float *_ev_array_s;
  double *_ev_array_d;
  int _results_eatom, _results_vatom;
  int _need_reduce;

#ifdef _LMP_INTEL_OFFLOAD
  double _balance_pair_time, _balance_other_time;
  int _offload_nlocal, _offload_nall, _offload_min_ghost, _offload_nghost;
  int _host_min_local, _host_min_ghost, _host_nall;
  int _host_used_local, _host_used_ghost, _sync_mode;
  int _separate_buffers, _offload_noghost, _separate_coi;
  bool _setup_time_cleared, _timers_allocated;
  void output_timing_data();
  FILE *_tscreen;

  IntelBuffers<float, float>::vec3_acc_t *_off_force_array_s;
  IntelBuffers<float, double>::vec3_acc_t *_off_force_array_m;
  IntelBuffers<double, double>::vec3_acc_t *_off_force_array_d;
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
  _alignvar(double _stopwatch_offload_neighbor[1], 64);
  _alignvar(double _stopwatch_offload_pair[1], 64);

  void _sync_main_arrays(const int prereverse);

  template <class ft> void reduce_results(ft *_noalias const f_in);

  template <class ft, class acc_t>
  inline void add_results(const ft *_noalias const f_in, const acc_t *_noalias const ev_global,
                          const int eatom, const int vatom, const int offload);

  template <class ft, class acc_t>
  inline void add_oresults(const ft *_noalias const f_in, const acc_t *_noalias const ev_global,
                           const int eatom, const int vatom, const int out_offset, const int nall);

  int _offload_affinity_balanced, _offload_threads, _offload_tpc;
#ifdef _LMP_INTEL_OFFLOAD
  int _max_offload_threads, _offload_cores, _offload_affinity_set;
  int _im_real_space_task;
  MPI_Comm _real_space_comm;
  template <class ft, class acc_t>
  inline void add_off_results(const ft *_noalias const f_in, const acc_t *_noalias const ev_global);
#endif
};

/* ---------------------------------------------------------------------- */

void FixIntel::get_buffern(const int offload, int &nlocal, int &nall, int &minlocal)
{
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
      if (force->newton)
        minlocal = _host_min_local;
      else
        minlocal = host_start_pair();
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

void FixIntel::add_result_array(IntelBuffers<double, double>::vec3_acc_t *f_in, double *ev_in,
                                const int offload, const int eatom, const int vatom,
                                const int rflag)
{
#ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_d = f_in;
    _off_ev_array_d = ev_in;
    if (_pair_hybrid_flag && force->pair->fdotr_is_set()) _sync_main_arrays(1);
    return;
  }
#endif

  _force_array_d = f_in;
  _ev_array_d = ev_in;
  _results_eatom = eatom;
  _results_vatom = vatom;
#ifndef _LMP_INTEL_OFFLOAD
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;
#endif

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float, double>::vec3_acc_t *f_in, double *ev_in,
                                const int offload, const int eatom, const int vatom,
                                const int rflag)
{
#ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_m = f_in;
    _off_ev_array_d = ev_in;
    if (_pair_hybrid_flag && force->pair->fdotr_is_set()) _sync_main_arrays(1);
    return;
  }
#endif

  _force_array_m = f_in;
  _ev_array_d = ev_in;
  _results_eatom = eatom;
  _results_vatom = vatom;
#ifndef _LMP_INTEL_OFFLOAD
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;
#endif

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float, float>::vec3_acc_t *f_in, float *ev_in,
                                const int offload, const int eatom, const int vatom,
                                const int rflag)
{
#ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _off_results_eatom = eatom;
    _off_results_vatom = vatom;
    _off_force_array_s = f_in;
    _off_ev_array_s = ev_in;
    if (_pair_hybrid_flag && force->pair->fdotr_is_set()) _sync_main_arrays(1);
    return;
  }
#endif

  _force_array_s = f_in;
  _ev_array_s = ev_in;
  _results_eatom = eatom;
  _results_vatom = vatom;
#ifndef _LMP_INTEL_OFFLOAD
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;
#endif

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

/* ---------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD

/* ---------------------------------------------------------------------- */

int FixIntel::offload_end_neighbor()
{
  if (_offload_balance < 0.0) {
    if (atom->nlocal < 2) error->one(FLERR, "Too few atoms for load balancing offload");
    double granularity = 1.0 / atom->nlocal;
    if (_balance_neighbor < granularity)
      _balance_neighbor = granularity + 1e-10;
    else if (_balance_neighbor > 1.0 - granularity)
      _balance_neighbor = 1.0 - granularity + 1e-10;
  }
  return _balance_neighbor * atom->nlocal;
}

int FixIntel::offload_end_pair()
{
  if (neighbor->ago == 0)
    return _balance_neighbor * atom->nlocal;
  else
    return _balance_pair * atom->nlocal;
}

/* ---------------------------------------------------------------------- */

double FixIntel::stop_watch(const int which)
{
  double elapsed = MPI_Wtime() - _stopwatch[which];
  _timers[which] += elapsed;
  return elapsed;
}

/* ---------------------------------------------------------------------- */

void FixIntel::balance_stamp()
{
  if (_offload_balance < 0.0) {
    double ct = MPI_Wtime();
    _balance_other_time = ct;
    _balance_pair_time = ct - _stopwatch[TIME_HOST_PAIR];
  }
}

/* ---------------------------------------------------------------------- */

void FixIntel::acc_timers()
{
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

void FixIntel::set_neighbor_host_sizes()
{
  _host_min_local = _overflow_flag[LMP_LOCAL_MIN];
  _host_min_ghost = _overflow_flag[LMP_GHOST_MIN];
  _host_used_local = atom->nlocal - _host_min_local;
  _host_used_ghost = _overflow_flag[LMP_GHOST_MAX] + 1 - _host_min_ghost;
  if (_host_used_ghost < 0) _host_used_ghost = 0;
  _host_nall = atom->nlocal + _host_used_ghost;
}

/* ---------------------------------------------------------------------- */

#endif

}    // namespace LAMMPS_NS

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

E: Currently, neighbor style BIN must be used with Intel package.

This is the only neighbor style that has been implemented for the Intel
package.

E: Currently, cannot use neigh_modify exclude with Intel package offload.

This is a current restriction of the Intel package when built for offload.

W: Unknown Intel Compiler Version

The compiler version used to build LAMMPS has not been tested with
offload to a coprocessor.

W: Unsupported Intel Compiler

The compiler version used to build LAMMPS is not supported when using
offload to a coprocessor. There could be performance or correctness
issues. Please use 14.0.1.106 or 15.1.133 or later.

E: Currently, cannot offload more than one intel style with hybrid.

Currently, when using offload, hybrid pair styles can only use the intel
suffix for one of the pair styles.

E: Cannot yet use hybrid styles with Intel offload.

The hybrid pair style configuration is not yet supported when using offload
within the Intel package. Support is limited to hybrid/overlay or a hybrid
style that does not require a skip list.

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

E: INTEL package requires same setting for newton bond and non-bond.

The newton setting must be the same for both pairwise and bonded forces.

E: Intel styles for bond/angle/dihedral/improper require intel pair style."

You cannot use the INTEL package for bond calculations without a
INTEL supported pair style.

E: Intel styles for kspace require intel pair style.

You cannot use the INTEL package for kspace calculations without a
INTEL supported pair style.

E: Cannot currently get per-atom virials with intel package.

The Intel package does not yet support per-atom virial calculation.

E: Too few atoms for load balancing offload.

When using offload to a coprocessor, each MPI task must have at least 2
atoms throughout the simulation.

E: Intel package requires fdotr virial with newton on.

This error can occur with a hybrid pair style that mixes styles that are
incompatible with the newton pair setting turned on. Try turning the
newton pair setting off.

E: Add -DLMP_INTEL_NBOR_COMPAT to build for special_bond exclusions with Intel

When using a manybody pair style, bonds/angles/dihedrals, and special_bond
exclusions, LAMMPS should be built with the above compile flag for compatible
results.

*/
