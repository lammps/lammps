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
  ~FixIntel() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  inline void min_setup(int in) override { setup(in); }
  void setup_pre_reverse(int eflag = 0, int vflag = 0) override;

  bool pair_hybrid_check();
  void pair_init_check(const bool cdmessage = false);
  void bond_init_check();
  void kspace_init_check();

  void pre_reverse(int eflag = 0, int vflag = 0) override;
  inline void min_pre_reverse(int eflag = 0, int vflag = 0) override { pre_reverse(eflag, vflag); }

  void post_force(int vflag) override;
  void post_run() override { _print_pkg_info = 1; }

  double memory_usage() override;

  typedef struct {
    double x, y, z;
  } lmp_ft;

  enum { PREC_MODE_SINGLE, PREC_MODE_MIXED, PREC_MODE_DOUBLE };

  inline int precision() { return _precision_mode; }
  inline IntelBuffers<float, float> *get_single_buffers() { return _single_buffers; }
  inline IntelBuffers<float, double> *get_mixed_buffers() { return _mixed_buffers; }
  inline IntelBuffers<double, double> *get_double_buffers() { return _double_buffers; }

  [[nodiscard]] int nbor_pack_width() const { return _nbor_pack_width; }
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
      return _p3m_table;
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
  int _pair_hybrid_zero, _hybrid_nonpair, _zero_master, _torque_flag;

 public:
  inline int *get_overflow_flag() { return _overflow_flag; }
  inline void add_result_array(IntelBuffers<double, double>::vec3_acc_t *f_in, double *ev_in,
                               const int eatom = 0, const int rflag = 0);
  inline void add_result_array(IntelBuffers<float, double>::vec3_acc_t *f_in, double *ev_in,
                               const int eatom = 0, const int rflag = 0);
  inline void add_result_array(IntelBuffers<float, float>::vec3_acc_t *f_in, float *ev_in,
                               const int eatom = 0, const int rflag = 0);
  inline void get_buffern(int &nlocal, int &nall, int &minlocal);

  inline void start_watch(const int /*which*/) {}
  inline double stop_watch(const int /*which*/) { return 0.0; }
  inline void balance_stamp() {}
  inline int separate_buffers() { return 0; }

 protected:
  int _overflow_flag[5];
  int _lrt, _p3m_table;

  IntelBuffers<float, float>::vec3_acc_t *_force_array_s;
  IntelBuffers<float, double>::vec3_acc_t *_force_array_m;
  IntelBuffers<double, double>::vec3_acc_t *_force_array_d;
  float *_ev_array_s;
  double *_ev_array_d;
  int _results_eatom;
  int _need_reduce;

  void _sync_main_arrays(const int prereverse);

  template <class ft> void reduce_results(ft *_noalias const f_in);

  template <class ft, class acc_t>
  inline void add_results(const ft *_noalias const f_in, const acc_t *_noalias const ev_global,
                          const int eatom);

  template <class ft, class acc_t>
  inline void add_oresults(const ft *_noalias const f_in, const acc_t *_noalias const ev_global,
                           const int eatom, const int out_offset, const int nall);
};

/* ---------------------------------------------------------------------- */

void FixIntel::get_buffern(int &nlocal, int &nall, int &minlocal)
{
  nall = atom->nlocal + atom->nghost;
  nlocal = atom->nlocal;
  minlocal = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<double, double>::vec3_acc_t *f_in, double *ev_in,
                                const int eatom, const int rflag)
{
  _force_array_d = f_in;
  _ev_array_d = ev_in;
  _results_eatom = eatom;
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float, double>::vec3_acc_t *f_in, double *ev_in,
                                const int eatom, const int rflag)
{
  _force_array_m = f_in;
  _ev_array_d = ev_in;
  _results_eatom = eatom;
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

/* ---------------------------------------------------------------------- */

void FixIntel::add_result_array(IntelBuffers<float, float>::vec3_acc_t *f_in, float *ev_in,
                                const int eatom, const int rflag)
{
  _force_array_s = f_in;
  _ev_array_s = ev_in;
  _results_eatom = eatom;
  if (rflag != 2 && _nthreads > 1 && force->newton) _need_reduce = 1;

  if (_overflow_flag[LMP_OVERFLOW])
    error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

  if (_pair_hybrid_flag > 1 || (_pair_hybrid_flag && force->pair->fdotr_is_set()))
    _sync_main_arrays(0);
}

}    // namespace LAMMPS_NS

#endif
#endif
