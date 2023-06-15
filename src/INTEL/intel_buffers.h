// clang-format off
/* -*- c++ -*- -------------------------------------------------------------
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

#ifndef LMP_INTEL_BUFFERS_H
#define LMP_INTEL_BUFFERS_H

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "intel_preprocess.h"

namespace LAMMPS_NS {

#define ATOM_T typename IntelBuffers<flt_t,acc_t>::atom_t
#define QUAT_T typename IntelBuffers<flt_t,acc_t>::quat_t
#define FORCE_T typename IntelBuffers<flt_t,acc_t>::vec3_acc_t

struct IntelNeighListPtrs {
  void *list_ptr;
  int *cnumneigh;
  int *numneighhalf;
  int size;
};

// May not need a separate force array for mixed/double
template <class flt_t, class acc_t>
class IntelBuffers {
 public:
  typedef struct { flt_t x,y,z; int w; } atom_t;
  typedef struct { flt_t w,i,j,k; } quat_t;
  typedef struct { flt_t x,y; } vec2_t;
  typedef struct { flt_t x,y,z,w; } vec3_t;
  typedef struct { flt_t x,y,z,w; } vec4_t;
  typedef struct { acc_t x,y,z,w; } vec3_acc_t;

  IntelBuffers(class LAMMPS *lmp_in);
  ~IntelBuffers();

  inline int get_stride(int nall) {
    int stride;
    IP_PRE_get_stride(stride, nall, sizeof(vec3_acc_t),
                         lmp->atom->torque);
    return stride;
  }

  template <class stype>
  inline int get_scalar_stride(const int n) {
    int stride;
    IP_PRE_get_stride(stride, n, sizeof(stype), 0);
    return stride;
  }

  void free_buffers();
  void free_nmax();
  inline void set_bininfo(int *atombin, int *binpacked) {
    _atombin = atombin;
    _binpacked = binpacked;
    _neigh_list_ptrs[0].numneighhalf = atombin;
  }

  inline void grow(const int nall, const int nlocal, const int nthreads,
                   const int offload_end) {
    if (nall >= _buf_size || nlocal >= _buf_local_size)
      _grow(nall, nlocal, nthreads, offload_end);
    #ifdef _LMP_INTEL_OFFLOAD
    if (lmp->atom->nmax > _host_nmax)
      _grow_nmax(offload_end);
    #endif
  }

  inline void free_all_nbor_buffers() {
    free_nbor_list();
    free_nmax();
    free_list_local();
    free_ncache();
    free_list_ptrs();
  }

  inline void grow_list(NeighList *list, const int nlocal, const int nthreads,
                        const int three_body, const int offload_end,
                        const int pack_width=1) {
    grow_list_local(list, three_body, offload_end);
    grow_nbor_list(list, nlocal, nthreads, offload_end, pack_width);
  }

  void free_list_local();
  inline void grow_list_local(NeighList *list, const int three_body,
                              const int offload_end) {
    _neigh_list_ptrs[0].list_ptr = (void *)list;
    if (list->get_maxlocal() > _off_map_listlocal)
      _grow_list_local(list, three_body, offload_end);
  }

  void free_ccache();
  void grow_ccache(const int off_flag, const int nthreads, const int width=1);
  inline int ccache_stride() { return _ccache_stride; }
  inline flt_t * get_ccachex() { return _ccachex; }
  inline flt_t * get_ccachey() { return _ccachey; }
  inline flt_t * get_ccachez() { return _ccachez; }
  inline flt_t * get_ccachew() { return _ccachew; }
  inline int * get_ccachei() { return _ccachei; }
  inline int * get_ccachej() { return _ccachej; }
  #ifdef LMP_USE_AVXCD
  inline int ccache_stride3() { return _ccache_stride3; }
  inline acc_t * get_ccachef() { return _ccachef; }
  #endif

  void free_ncache();
  void grow_ncache(const int off_flag, const int nthreads);
  void grow_ncachetag(const int off_flag, const int nthreads);
  inline int ncache_stride() { return _ncache_stride; }
  inline flt_t * get_ncachex() { return _ncachex; }
  inline flt_t * get_ncachey() { return _ncachey; }
  inline flt_t * get_ncachez() { return _ncachez; }
  inline int * get_ncachej() { return _ncachej; }
  inline int * get_ncachejtype() { return _ncachejtype; }
  inline tagint * get_ncachetag() { return _ncachetag; }

  inline int get_max_nbors() {
    int mn = lmp->neighbor->oneatom * sizeof(int) /
        (INTEL_ONEATOM_FACTOR * INTEL_DATA_ALIGN);
    return mn * INTEL_DATA_ALIGN / sizeof(int);
  }

  void free_nbor_list();

  inline void grow_nbor_list(NeighList *list, const int nlocal,
                             const int nthreads, const int offload_end,
                             const int pack_width) {
    if (nlocal > _list_alloc_atoms)
      _grow_nbor_list(list, nlocal, nthreads, offload_end, pack_width);
  }

  void set_ntypes(const int ntypes, const int use_ghost_cut = 1);

  inline int * intel_list(const NeighList * /*list*/) { return _list_alloc; }
  inline int * get_atombin() { return _atombin; }
  inline int * get_binpacked() { return _binpacked; }
  inline int * cnumneigh() { return _neigh_list_ptrs[0].cnumneigh; }
  inline void get_list_data3(const NeighList *list, int *&numneighhalf,
                             int *&cnumneigh) {
    for (int i = 0; i < _n_list_ptrs; i++)
      if ((void *)list == _neigh_list_ptrs[i].list_ptr) {
        numneighhalf = _neigh_list_ptrs[i].numneighhalf;
        cnumneigh = _neigh_list_ptrs[i].cnumneigh;
      }
  }
  void grow_data3(NeighList *list, int *&numneighhalf, int *&cnumneigh);
  void free_list_ptrs();

  inline atom_t * get_x(const int offload = 1) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_separate_buffers && offload == 0) return _host_x;
    #endif
    return _x;
  }
  inline flt_t * get_q(const int offload = 1) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_separate_buffers && offload == 0) return _host_q;
    #endif
    return _q;
  }
  inline quat_t * get_quat(const int offload = 1) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_separate_buffers && offload == 0) return _host_quat;
    #endif
    return _quat;
  }
  inline vec3_acc_t * get_f() { return _f; }
  inline acc_t * get_ev_global() { return _ev_global; }
  inline acc_t * get_ev_global_host() { return _ev_global_host; }
  inline void zero_ev()
    { for (int i = 0; i < 8; i++) _ev_global[i] = _ev_global_host[i] = 0.0; }
  inline flt_t ** get_cutneighsq() { return _cutneighsq; }
  inline flt_t ** get_cutneighghostsq() { return _cutneighghostsq; }
  inline int get_off_threads() { return _off_threads; }
  #ifdef _LMP_INTEL_OFFLOAD
  inline void set_off_params(const int n, const int cop,
                             const int separate_buffers)
    { _off_threads = n; _cop = cop; _separate_buffers = separate_buffers; }
  inline vec3_acc_t * get_off_f() { return _off_f; }
  #endif

  inline void thr_pack(const int ifrom, const int ito, const int ago) {
    if (ago == 0) {
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = lmp->atom->x[i][0];
        _x[i].y = lmp->atom->x[i][1];
        _x[i].z = lmp->atom->x[i][2];
        _x[i].w = lmp->atom->type[i];
      }
      if (lmp->atom->q != nullptr)
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++)
          _q[i] = lmp->atom->q[i];
    } else {
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = lmp->atom->x[i][0];
        _x[i].y = lmp->atom->x[i][1];
        _x[i].z = lmp->atom->x[i][2];
      }
    }
  }

  inline void thr_pack_q(const int ifrom, const int ito) {
    if (lmp->atom->q != nullptr)
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int i = ifrom; i < ito; i++)
        _q[i] = lmp->atom->q[i];
  }

  #ifndef _LMP_INTEL_OFFLOAD
  void fdotr_reduce_l5(const int lf, const int lt, const int nthreads,
                       const int f_stride, acc_t &ov0, acc_t &ov1,
                       acc_t &ov2, acc_t &ov3, acc_t &ov4, acc_t &ov5);
  void fdotr_reduce(const int nall, const int nthreads, const int f_stride,
                    acc_t &ov0, acc_t &ov1, acc_t &ov2, acc_t &ov3,
                    acc_t &ov4, acc_t &ov5);
  #endif

  #ifdef _LMP_INTEL_OFFLOAD
  inline void thr_pack_cop(const int ifrom, const int ito,
                           const int offset, const bool dotype = false) {
    double ** x = lmp->atom->x + offset;
    if (dotype == false) {
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = x[i][0];
        _x[i].y = x[i][1];
        _x[i].z = x[i][2];
      }
    } else {
      int *type = lmp->atom->type + offset;
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = x[i][0];
        _x[i].y = x[i][1];
        _x[i].z = x[i][2];
        _x[i].w = type[i];
      }
    }
  }

  inline void thr_pack_host(const int ifrom, const int ito,
                            const int offset) {
    double ** x = lmp->atom->x + offset;
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma ivdep
    #endif
    for (int i = ifrom; i < ito; i++) {
      _host_x[i].x = x[i][0];
      _host_x[i].y = x[i][1];
      _host_x[i].z = x[i][2];
    }
  }

  inline void pack_sep_from_single(const int host_min_local,
                                   const int used_local,
                                   const int host_min_ghost,
                                   const int used_ghost) {
    memcpy(_host_x + host_min_local, _x + host_min_local,
           used_local * sizeof(atom_t));
    memcpy(_host_x + host_min_local + used_local, _x + host_min_ghost,
           used_ghost * sizeof(atom_t));
    int nall = used_local + used_ghost + host_min_local;
    _host_x[nall].x = INTEL_BIGP;
    _host_x[nall].y = INTEL_BIGP;
    _host_x[nall].z = INTEL_BIGP;
    _host_x[nall].w = 1;
    if (lmp->atom->q != nullptr) {
      memcpy(_host_q + host_min_local, _q + host_min_local,
             used_local * sizeof(flt_t));
      memcpy(_host_q + host_min_local + used_local, _q + host_min_ghost,
             used_ghost * sizeof(flt_t));
    }
  }
  #endif

  inline int need_tag() { return _need_tag; }
  inline void need_tag(const int nt) { _need_tag = nt; }

  double memory_usage(const int nthreads);

  tagint _special_holder;
  int _nspecial_holder;

 protected:
  LAMMPS *lmp;
  atom_t *_x;
  flt_t *_q;
  quat_t *_quat;
  vec3_acc_t * _f;
  int _off_threads, _off_map_listlocal;

  int _list_alloc_atoms;
  int *_list_alloc, *_cnumneigh, *_atombin, *_binpacked;

  IntelNeighListPtrs *_neigh_list_ptrs;
  int _n_list_ptrs, _max_list_ptrs;

  flt_t **_cutneighsq, **_cutneighghostsq;
  int _ntypes;

  int _ccache_stride;
  flt_t *_ccachex, *_ccachey, *_ccachez, *_ccachew;
  int *_ccachei, *_ccachej;

  int _ncache_stride, _ncache_alloc;
  flt_t *_ncachex, *_ncachey, *_ncachez;
  int *_ncachej, *_ncachejtype;
  tagint *_ncachetag;

  int _need_tag, _host_nmax;

  #ifdef LMP_USE_AVXCD
  int _ccache_stride3;
  acc_t * _ccachef;
  #endif

  #ifdef _LMP_INTEL_OFFLOAD
  int _separate_buffers;
  atom_t *_host_x;
  flt_t *_host_q;
  quat_t *_host_quat;
  vec3_acc_t *_off_f;
  int _off_map_nmax, _cop, _off_ccache, _off_ncache;
  int *_off_map_ilist, *_off_map_nspecial;
  tagint *_off_map_tag, *_off_map_special;
  int **_off_map_firstneigh, *_off_map_numneigh;
  bool _off_list_alloc;
  #endif

  int _buf_size, _buf_local_size;
  _alignvar(acc_t _ev_global[8],64);
  _alignvar(acc_t _ev_global_host[8],64);

  void _grow(const int nall, const int nlocal, const int nthreads,
             const int offload_end);
  void _grow_nmax(const int offload_end);
  void _grow_list_local(NeighList *list, const int three_body,
                        const int offload_end);
  void _grow_nbor_list(NeighList *list, const int nlocal, const int nthreads,
                       const int offload_end, const int pack_width);
};

}

#endif
