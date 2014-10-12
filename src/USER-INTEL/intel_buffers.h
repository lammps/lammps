/* -*- c++ -*- -------------------------------------------------------------
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

#ifndef LMP_INTEL_BUFFERS_H
#define LMP_INTEL_BUFFERS_H

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "intel_preprocess.h"
#include <cstring>

namespace LAMMPS_NS {

#define ATOM_T typename IntelBuffers<flt_t,acc_t>::atom_t
#define QUAT_T typename IntelBuffers<flt_t,acc_t>::quat_t
#define FORCE_T typename IntelBuffers<flt_t,acc_t>::vec3_acc_t

// May not need a separate force array for mixed/double
template <class flt_t, class acc_t>
class IntelBuffers {
 public:
  typedef struct { flt_t x,y,z; int w; } atom_t;
  typedef struct { flt_t w,i,j,k; } quat_t;
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

  void free_buffers();

  inline void grow(const int nall, const int nlocal, const int nthreads,
                   const int offload_end) {
    if (nall >= _buf_size || nlocal >= _buf_local_size)
      _grow(nall, nlocal, nthreads, offload_end);
  }

  inline void free_all_nbor_buffers() {
    free_nbor_list();
    free_nmax();
    free_binhead();
    free_local();
  }

  inline void grow_nbor(NeighList *list, const int nlocal,
                        const int offload_end) {
    grow_local(list, offload_end);
    if (offload_end) {
      grow_nmax();
      grow_binhead();
    }
    grow_nbor_list(list, nlocal, offload_end);
  }

  void free_nmax();

  inline void grow_nmax() {
    #ifdef _LMP_INTEL_OFFLOAD
    if (lmp->atom->nmax > _off_map_nmax)
      _grow_nmax();
    #endif
  }

  void free_local();

  inline void grow_local(NeighList *list, const int offload_end) {
    if (list->get_maxlocal() > _off_map_maxlocal)
      _grow_local(list, offload_end);
  }

  void free_binhead();
  
  inline void grow_binhead() {
    #ifdef _LMP_INTEL_OFFLOAD
    if (lmp->neighbor->maxhead > _off_map_maxhead)
      _grow_binhead();
    #endif
  }

  inline int get_max_nbors() {
    int mn = lmp->neighbor->oneatom * sizeof(int) /
        (INTEL_ONEATOM_FACTOR * INTEL_DATA_ALIGN);
    return mn * INTEL_DATA_ALIGN / sizeof(int);
  }
  
  void free_nbor_list();

  inline void grow_nbor_list(NeighList *list, const int nlocal,
                             const int offload_end) {
    if (nlocal > _list_alloc_atoms)
      _grow_nbor_list(list, nlocal, offload_end);
    #ifdef _LMP_INTEL_OFFLOAD
    else if (offload_end > 0 && _off_map_stencil != list->stencil)
      _grow_stencil(list);
    #endif
  }

  void set_ntypes(const int ntypes);

  inline int * firstneigh(const NeighList *list) { return _list_alloc; }
  inline int * cnumneigh(const NeighList *list) { return _cnumneigh; }

  inline int * get_atombin() { return _atombin; }
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
  inline int get_off_threads() { return _off_threads; }
  #ifdef _LMP_INTEL_OFFLOAD
  inline void set_off_params(const int n, const int cop, 
			     const int separate_buffers) 
    { _off_threads = n; _cop = cop; _separate_buffers = separate_buffers; } 
  inline vec3_acc_t * get_off_f() { return _off_f; }
  #endif

  inline void thr_pack(const int ifrom, const int ito, const int ago) {
    if (ago == 0) {
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = lmp->atom->x[i][0];
        _x[i].y = lmp->atom->x[i][1];
        _x[i].z = lmp->atom->x[i][2];
        _x[i].w = lmp->atom->type[i];
      }
      if (lmp->atom->q != NULL)
        for (int i = ifrom; i < ito; i++)
          _q[i] = lmp->atom->q[i];
    } else {
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = lmp->atom->x[i][0];
        _x[i].y = lmp->atom->x[i][1];
        _x[i].z = lmp->atom->x[i][2];
      }
    }
  }

  #ifdef _LMP_INTEL_OFFLOAD
  inline void thr_pack_cop(const int ifrom, const int ito, 
			   const int offset, const bool dotype = false) {
    double ** x = lmp->atom->x + offset;
    if (dotype == false) {
      #pragma vector nontemporal
      for (int i = ifrom; i < ito; i++) {
        _x[i].x = x[i][0];
        _x[i].y = x[i][1];
        _x[i].z = x[i][2];
      }
    } else {
      int *type = lmp->atom->type + offset;
      #pragma vector nontemporal
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
    if (lmp->atom->q != NULL) {
      memcpy(_host_q + host_min_local, _q + host_min_local,
	     used_local * sizeof(flt_t));
      memcpy(_host_q + host_min_local + used_local, _q + host_min_ghost,
	     used_ghost * sizeof(flt_t));
    }
  }
  #endif

  double memory_usage(const int nthreads);

  tagint _special_holder;
  int _nspecial_holder;

 protected:
  LAMMPS *lmp;
  atom_t *_x;
  flt_t *_q;
  quat_t *_quat;
  vec3_acc_t * _f;
  int _off_threads, _off_map_maxlocal;

  int _list_alloc_atoms;
  int * _list_alloc;
  int * _cnumneigh;
  int * _atombin;

  flt_t **_cutneighsq;
  int _ntypes;

  #ifdef _LMP_INTEL_OFFLOAD
  int _separate_buffers;
  atom_t *_host_x;
  flt_t *_host_q;
  quat_t *_host_quat;
  vec3_acc_t *_off_f;
  int _off_map_nmax, _off_map_maxhead, _cop;
  int *_off_map_ilist;
  int *_off_map_stencil, *_off_map_special, *_off_map_nspecial, *_off_map_tag;
  int *_off_map_binhead, *_off_map_bins, *_off_map_numneigh;
  bool _off_list_alloc;
  #endif
  
  int _buf_size, _buf_local_size;
  _alignvar(acc_t _ev_global[8],64);
  _alignvar(acc_t _ev_global_host[8],64);

  void _grow(const int nall, const int nlocal, const int nthreads,
	     const int offload_end);
  void _grow_nmax();
  void _grow_local(NeighList *list, const int offload_end);
  void _grow_binhead();
  void _grow_nbor_list(NeighList *list, const int nlocal,
                       const int offload_end);
  void _grow_stencil(NeighList *list);
};

}

#endif
