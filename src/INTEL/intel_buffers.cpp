// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "intel_buffers.h"

#include "force.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
IntelBuffers<flt_t, acc_t>::IntelBuffers(class LAMMPS *lmp_in) :
    lmp(lmp_in), _x(nullptr), _q(nullptr), _quat(nullptr), _f(nullptr), _off_threads(0),
    _n_list_ptrs(1), _max_list_ptrs(4), _buf_size(0), _buf_local_size(0) {
  _neigh_list_ptrs = new IntelNeighListPtrs[_max_list_ptrs];
  _neigh_list_ptrs[0].cnumneigh = nullptr;
  _list_alloc_atoms = 0;
  _ntypes = 0;
  _off_map_listlocal = 0;
  _ccachex = nullptr;
  _ncache_alloc = 0;
  _ncachetag = nullptr;
  _cutneighsq = nullptr;
  _cutneighghostsq = nullptr;
  _need_tag = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  _separate_buffers = 0;
  _off_f = 0;
  _off_map_ilist = 0;
  _off_map_nmax = 0;
  _off_list_alloc = false;
  _off_threads = 0;
  _off_ccache = 0;
  _off_ncache = 0;
  _host_nmax = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
IntelBuffers<flt_t, acc_t>::~IntelBuffers()
{
  free_buffers();
  free_all_nbor_buffers();
  free_ccache();
  set_ntypes(0);
  delete []_neigh_list_ptrs;
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_buffers()
{
  if (_buf_size > 0) {
    atom_t * x = get_x();
    flt_t * q = get_q();
    quat_t * quat = get_quat();

    #ifdef _LMP_INTEL_OFFLOAD
    vec3_acc_t * f_start = get_off_f();
    if (f_start != 0) {
      acc_t * ev_global = get_ev_global();
      if (ev_global != 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(x:alloc_if(0) free_if(1)) \
          nocopy(f_start:alloc_if(0) free_if(1)) \
          nocopy(ev_global:alloc_if(0) free_if(1))
      }

      if (q != 0) {
        #pragma offload_transfer target (mic:_cop) \
          nocopy(q:alloc_if(0) free_if(1))
      }
      if (quat != 0) {
        #pragma offload_transfer target (mic:_cop) \
          nocopy(quat:alloc_if(0) free_if(1))
      }
      lmp->memory->destroy(f_start);
    }

    if (_separate_buffers) {
      lmp->memory->destroy(_host_x);
      if (q != 0) lmp->memory->destroy(_host_q);
      if (quat != 0) lmp->memory->destroy(_host_quat);
    }
    #endif

    lmp->memory->destroy(x);
    if (q != nullptr) lmp->memory->destroy(q);
    if (quat != nullptr) lmp->memory->destroy(quat);
    lmp->memory->destroy(_f);
    _buf_size = _buf_local_size = 0;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow(const int nall, const int nlocal,
                                       const int nthreads,
                                       const int offload_end)
{
  free_buffers();
  _buf_size = static_cast<double>(nall) * 1.1 + 1;
  if (lmp->force->newton_pair)
    _buf_local_size = _buf_size;
  else
    _buf_local_size = static_cast<double>(nlocal) * 1.1 + 1;
  const int f_stride = get_stride(_buf_local_size);
  lmp->memory->create(_x, _buf_size,"intel_x");
  if (lmp->atom->q != nullptr)
    lmp->memory->create(_q, _buf_size, "intel_q");
  if (lmp->atom->ellipsoid != nullptr)
    lmp->memory->create(_quat, _buf_size, "intel_quat");
  #ifdef _LMP_INTEL_OFFLOAD
  if (lmp->force->newton_pair)
  #else
  if (lmp->force->newton_pair || lmp->atom->molecular)
  #endif
    lmp->memory->create(_f, f_stride * nthreads, "intel_f");
  else
    lmp->memory->create(_f, f_stride, "intel_f");

  #ifdef _LMP_INTEL_OFFLOAD
  if (_separate_buffers) {
    lmp->memory->create(_host_x, _buf_size,"intel_host_x");
    if (lmp->atom->q != nullptr)
      lmp->memory->create(_host_q, _buf_size, "intel_host_q");
    if (lmp->atom->ellipsoid != nullptr)
      lmp->memory->create(_host_quat, _buf_size, "intel_host_quat");
  }

  if (offload_end > 0) {
    int fm;
    if (lmp->force->newton_pair) fm = _off_threads;
    else fm = 1;
    lmp->memory->create(_off_f, f_stride * fm, "intel_off_f");
    const atom_t * const x = get_x();
    const flt_t * const q = get_q();
    const vec3_acc_t * f_start = get_off_f();
    acc_t * ev_global = get_ev_global();
    if (lmp->atom->q != nullptr) {
      if (x != nullptr && q != nullptr && f_start != nullptr && ev_global != nullptr) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(x,q:length(_buf_size) alloc_if(1) free_if(0)) \
          nocopy(f_start:length(f_stride*fm) alloc_if(1) free_if(0))\
          nocopy(ev_global:length(8) alloc_if(1) free_if(0))
      }
    } else {
      if (x != nullptr && f_start != nullptr && ev_global != nullptr) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(x:length(_buf_size) alloc_if(1) free_if(0)) \
          nocopy(f_start:length(f_stride*fm) alloc_if(1) free_if(0))\
          nocopy(ev_global:length(8) alloc_if(1) free_if(0))
      }
    }
    if (lmp->atom->ellipsoid != nullptr) {
      const quat_t * const quat = get_quat();
      if (quat != nullptr) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(quat:length(_buf_size) alloc_if(1) free_if(0))
      }
    }
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_nmax()
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_off_map_nmax > 0) {
    const tagint * tag = _off_map_tag;
    const tagint * special = _off_map_special;
    const int * nspecial = _off_map_nspecial;
    #pragma offload_transfer target(mic:_cop) \
      nocopy(tag:alloc_if(0) free_if(1)) \
      nocopy(special,nspecial:alloc_if(0) free_if(1))
    _off_map_nmax = 0;
    _host_nmax = 0;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_nmax(const int offload_end)
{
  #ifdef _LMP_INTEL_OFFLOAD
  free_nmax();
  int size = lmp->atom->nmax;
  _host_nmax = size;

  if (!offload_end) return;
  tagint *special;
  int *nspecial;
  int tag_length, special_length, nspecial_length;
  if (lmp->atom->molecular) {
    special = lmp->atom->special[0];
    nspecial = lmp->atom->nspecial[0];
    special_length = size * lmp->atom->maxspecial;
    nspecial_length = size * 3;
  } else {
    special = &_special_holder;
    nspecial = &_nspecial_holder;
    special_length = 1;
    nspecial_length = 1;
  }
  if (_need_tag)
    tag_length = size;
  else
    tag_length = 1;
  tagint *tag = lmp->atom->tag;
  #pragma offload_transfer target(mic:_cop) \
    nocopy(tag:length(tag_length) alloc_if(1) free_if(0)) \
    nocopy(special:length(special_length) alloc_if(1) free_if(0)) \
    nocopy(nspecial:length(nspecial_length) alloc_if(1) free_if(0))
  _off_map_tag = tag;
  _off_map_special = special;
  _off_map_nspecial = nspecial;
  _off_map_nmax = size;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_list_local()
{
  if (_off_map_listlocal > 0) {
    if (_neigh_list_ptrs[0].cnumneigh) {
      int * cnumneigh = _neigh_list_ptrs[0].cnumneigh;
      _neigh_list_ptrs[0].cnumneigh = nullptr;
      #ifdef _LMP_INTEL_OFFLOAD
      if (_off_map_ilist != nullptr) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cnumneigh:alloc_if(0) free_if(1))
      }
      #endif
      lmp->memory->destroy(cnumneigh);
    }

    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_map_ilist != nullptr) {
      const int * ilist = _off_map_ilist;
      const int * numneigh = _off_map_numneigh;
      const int ** firstneigh = (const int **)_off_map_firstneigh;
      _off_map_ilist = nullptr;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ilist,firstneigh,numneigh:alloc_if(0) free_if(1))
    }
    #endif
    _off_map_listlocal = 0;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_list_ptrs()
{
  for (int list_num = 1; list_num < _n_list_ptrs; list_num++) {
    if (_neigh_list_ptrs[list_num].size) {
      lmp->memory->destroy(_neigh_list_ptrs[list_num].cnumneigh);
      lmp->memory->destroy(_neigh_list_ptrs[list_num].numneighhalf);
    }
    _neigh_list_ptrs[list_num].size = 0;
    _neigh_list_ptrs[list_num].list_ptr = nullptr;
  }
  _n_list_ptrs = 1;
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::grow_data3(NeighList *list, int *&numneighhalf, int *&cnumneigh)
{
  const int size = list->get_maxlocal();
  int list_num;
  for (list_num = 0; list_num < _n_list_ptrs; list_num++)
    if (_neigh_list_ptrs[list_num].list_ptr == (void*)list) break;
  if (list_num == _n_list_ptrs) {
    if (_n_list_ptrs == _max_list_ptrs) {
      _max_list_ptrs *= 2;
      auto new_list = new IntelNeighListPtrs[_max_list_ptrs];
      for (int i = 0; i < _n_list_ptrs; i++) new_list[i] = _neigh_list_ptrs[i];
      delete []_neigh_list_ptrs;
      _neigh_list_ptrs = new_list;
    }
    _neigh_list_ptrs[list_num].list_ptr = (void *)list;
    _neigh_list_ptrs[list_num].size = 0;
    _n_list_ptrs++;
  }
  if (size > _neigh_list_ptrs[list_num].size) {
    if (_neigh_list_ptrs[list_num].size) {
      lmp->memory->destroy(_neigh_list_ptrs[list_num].cnumneigh);
      lmp->memory->destroy(_neigh_list_ptrs[list_num].numneighhalf);
    }
    lmp->memory->create(_neigh_list_ptrs[list_num].cnumneigh, size, "_cnumneigh");
    lmp->memory->create(_neigh_list_ptrs[list_num].numneighhalf, size, "_cnumneigh");
    _neigh_list_ptrs[list_num].size = size;
  }
  numneighhalf = _neigh_list_ptrs[list_num].numneighhalf;
  cnumneigh = _neigh_list_ptrs[list_num].cnumneigh;
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_list_local(NeighList *list, const int three_body,
                                                  const int offload_end)
{
  free_list_local();
  int size = list->get_maxlocal();
  _off_map_listlocal = size;
  if (three_body)
    lmp->memory->create(_neigh_list_ptrs[0].cnumneigh, size, "_cnumneigh");

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_end > 0) {
    int tb_size = size;
    if (three_body == 0) {
      lmp->memory->create(_neigh_list_ptrs[0].cnumneigh, 16, "_cnumneigh");
      tb_size = 16;
    }
    int ** firstneigh = list->firstneigh;
    int * numneigh = list->numneigh;
    int * ilist = list->ilist;
    int * cnumneigh = _neigh_list_ptrs[0].cnumneigh;
    #pragma offload_transfer target(mic:_cop) \
      nocopy(ilist:length(size) alloc_if(1) free_if(0)) \
      nocopy(firstneigh:length(size) alloc_if(1) free_if(0)) \
      nocopy(numneigh:length(size) alloc_if(1) free_if(0)) \
      nocopy(cnumneigh:length(tb_size) alloc_if(1) free_if(0))
    _off_map_ilist = ilist;
    _off_map_firstneigh = firstneigh;
    _off_map_numneigh = numneigh;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_nbor_list()
{
  if (_list_alloc_atoms > 0) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_list_alloc) {
      int * list_alloc = _list_alloc;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(list_alloc:alloc_if(0) free_if(1))
      _off_list_alloc = false;
    }
    #endif
    lmp->memory->destroy(_list_alloc);
    _list_alloc_atoms = 0;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_nbor_list(NeighList * /*list*/,
                                                 const int nlocal,
                                                 const int nthreads,
                                                 const int offload_end,
                                                 const int pack_width)
{
  free_nbor_list();
  _list_alloc_atoms = 1.10 * nlocal;
  int nt = MAX(nthreads, _off_threads);
  int list_alloc_size = (_list_alloc_atoms + nt * 2 + pack_width - 1) *
    get_max_nbors();
  lmp->memory->create(_list_alloc, list_alloc_size, "_list_alloc");
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_end > 0) {
    int * list_alloc =_list_alloc;

    if (list_alloc != nullptr) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(list_alloc:length(list_alloc_size) alloc_if(1) free_if(0))
      _off_list_alloc = true;
    }
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_ccache()
{
  if (_ccachex) {
    flt_t *ccachex = _ccachex;
    flt_t *ccachey = _ccachey;
    flt_t *ccachez = _ccachez;
    flt_t *ccachew = _ccachew;
    int *ccachei = _ccachei;
    int *ccachej = _ccachej;
    #ifdef LMP_USE_AVXCD
    acc_t *ccachef = _ccachef;
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_ccache) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ccachex,ccachey,ccachez,ccachew:alloc_if(0) free_if(1)) \
        nocopy(ccachei,ccachej:alloc_if(0) free_if(1))

      #ifdef LMP_USE_AVXCD
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ccachef:alloc_if(0) free_if(1))
      #endif
    }
    _off_ccache = 0;
    #endif

    lmp->memory->destroy(ccachex);
    lmp->memory->destroy(ccachey);
    lmp->memory->destroy(ccachez);
    lmp->memory->destroy(ccachew);
    lmp->memory->destroy(ccachei);
    lmp->memory->destroy(ccachej);
    #ifdef LMP_USE_AVXCD
    lmp->memory->destroy(ccachef);
    #endif

    _ccachex = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::grow_ccache(const int off_flag,
        const int nthreads,
        const int width)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_ccachex && off_flag && _off_ccache == 0)
    free_ccache();
  #endif
  if (_ccachex)
    return;

  const int nsize = get_max_nbors() * width;
  int esize = MIN(sizeof(int), sizeof(flt_t));
  IP_PRE_get_stride(_ccache_stride, nsize, esize, 0);
  int nt = MAX(nthreads, _off_threads);
  const int vsize = _ccache_stride * nt;

  lmp->memory->create(_ccachex, vsize , "_ccachex");
  lmp->memory->create(_ccachey, vsize, "_ccachey");
  lmp->memory->create(_ccachez, vsize, "_ccachez");
  lmp->memory->create(_ccachew, vsize, "_ccachew");
  lmp->memory->create(_ccachei, vsize, "_ccachei");
  lmp->memory->create(_ccachej, vsize, "_ccachej");
  #ifdef LMP_USE_AVXCD
  IP_PRE_get_stride(_ccache_stride3, nsize * 3, sizeof(acc_t), 0);
  lmp->memory->create(_ccachef, _ccache_stride3 * nt, "_ccachef");
  #endif
  memset(_ccachei, 0, vsize * sizeof(int));
  memset(_ccachej, 0, vsize * sizeof(int));

  #ifdef _LMP_INTEL_OFFLOAD
  if (off_flag) {
    flt_t *ccachex = _ccachex;
    flt_t *ccachey = _ccachey;
    flt_t *ccachez = _ccachez;
    flt_t *ccachew = _ccachew;
    int *ccachei = _ccachei;
    int *ccachej = _ccachej;

    if (ccachex != nullptr && ccachey !=nullptr && ccachez != nullptr &&
        ccachew != nullptr && ccachei != nullptr && ccachej !=nullptr) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ccachex,ccachey:length(vsize) alloc_if(1) free_if(0)) \
        nocopy(ccachez,ccachew:length(vsize) alloc_if(1) free_if(0)) \
        in(ccachei:length(vsize) alloc_if(1) free_if(0)) \
        in(ccachej:length(vsize) alloc_if(1) free_if(0))
    }
    #ifdef LMP_USE_AVXCD
    if (ccachef != nullptr) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ccachef:length(_ccache_stride3 * nt) alloc_if(1) free_if(0))
    }
    #endif
    _off_ccache = 1;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_ncache()
{
  if (_ncache_alloc) {
    flt_t *ncachex = _ncachex;
    flt_t *ncachey = _ncachey;
    flt_t *ncachez = _ncachez;
    int *ncachej = _ncachej;
    int *ncachejtype = _ncachejtype;
    tagint *ncachetag = _ncachetag;

    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_ncache) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ncachex,ncachey,ncachez,ncachej:alloc_if(0) free_if(1)) \
        nocopy(ncachejtype:alloc_if(0) free_if(1))
      if (ncachetag) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(ncachetag:alloc_if(0) free_if(1))
      }
    }
    _off_ncache = 0;
    #endif

    lmp->memory->destroy(ncachex);
    lmp->memory->destroy(ncachey);
    lmp->memory->destroy(ncachez);
    lmp->memory->destroy(ncachej);
    lmp->memory->destroy(ncachejtype);
    if (ncachetag)
      lmp->memory->destroy(ncachetag);
    _ncache_alloc = 0;
    _ncachetag = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::grow_ncache(const int off_flag,
                                             const int nthreads)
{
  const int nsize = get_max_nbors() * 3;
  int esize = MIN(sizeof(int), sizeof(flt_t));
  IP_PRE_get_stride(_ncache_stride, nsize, esize, 0);
  int nt = MAX(nthreads, _off_threads);
  const int vsize = _ncache_stride * nt;

  if (_ncache_alloc) {
    if (vsize > _ncache_alloc || (need_tag() && _ncachetag == nullptr))
      free_ncache();
    #ifdef _LMP_INTEL_OFFLOAD
    else if (off_flag && _off_ncache == 0)
      free_ncache();
    #endif
    else
      return;
  }

  lmp->memory->create(_ncachex, vsize, "_ncachex");
  lmp->memory->create(_ncachey, vsize, "_ncachey");
  lmp->memory->create(_ncachez, vsize, "_ncachez");
  lmp->memory->create(_ncachej, vsize, "_ncachej");
  lmp->memory->create(_ncachejtype, vsize, "_ncachejtype");
  if (need_tag())
    lmp->memory->create(_ncachetag, vsize, "_ncachetag");

  _ncache_alloc = vsize;

  #ifdef _LMP_INTEL_OFFLOAD
  if (off_flag) {
    flt_t *ncachex = _ncachex;
    flt_t *ncachey = _ncachey;
    flt_t *ncachez = _ncachez;
    int *ncachej = _ncachej;
    int *ncachejtype = _ncachejtype;

    if (ncachex != nullptr && ncachey !=nullptr && ncachez != nullptr &&
        ncachej != nullptr && ncachejtype != nullptr) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ncachex,ncachey:length(vsize) alloc_if(1) free_if(0)) \
        nocopy(ncachez,ncachej:length(vsize) alloc_if(1) free_if(0)) \
        nocopy(ncachejtype:length(vsize) alloc_if(1) free_if(0))
    }
    int tsize = vsize;
    if (!need_tag()) {
      tsize = 16;
      lmp->memory->create(_ncachetag, tsize, "_ncachetag");
    }
    tagint *ncachetag = _ncachetag;
    #pragma offload_transfer target(mic:_cop)                   \
      nocopy(ncachetag:length(tsize) alloc_if(1) free_if(0))
    _off_ncache = 1;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

#ifndef _LMP_INTEL_OFFLOAD
template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::fdotr_reduce_l5(const int lf, const int lt,
    const int nthreads, const int f_stride, acc_t &ov0, acc_t &ov1,
    acc_t &ov2, acc_t &ov3, acc_t &ov4, acc_t &ov5)
{
  IP_PRE_fdotr_acc_force_l5(lf, lt, 0, nthreads, _f, f_stride, _x, ov0,
                            ov1, ov2, ov3, ov4, ov5);
}
#endif

/* ---------------------------------------------------------------------- */

#ifndef _LMP_INTEL_OFFLOAD
template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::fdotr_reduce(const int nall,
    const int nthreads, const int f_stride, acc_t &ov0, acc_t &ov1,
    acc_t &ov2, acc_t &ov3, acc_t &ov4, acc_t &ov5)
{
  int iifrom, iito, tid;
  IP_PRE_fdotr_acc_force(nall, 0, nthreads, _f, f_stride, _x, 0, 2,
                         ov0, ov1, ov2, ov3, ov4, ov5);
}
#endif

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::set_ntypes(const int ntypes,
                                            const int use_ghost_cut)
{
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * cutneighsqo = _cutneighsq[0];
      if (_off_threads > 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cutneighsqo:alloc_if(0) free_if(1))
      }
      flt_t * cutneighghostsqo;
      if (_cutneighghostsq && _off_threads > 0) {
        cutneighghostsqo = _cutneighghostsq[0];
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cutneighghostsqo:alloc_if(0) free_if(1))
      }
      #endif
      lmp->memory->destroy(_cutneighsq);
      if (_cutneighghostsq != nullptr) lmp->memory->destroy(_cutneighghostsq);
    }
    if (ntypes > 0) {
      lmp->memory->create(_cutneighsq, ntypes, ntypes, "_cutneighsq");
      if (use_ghost_cut)
        lmp->memory->create(_cutneighghostsq, ntypes, ntypes,
                            "_cutneighghostsq");
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * cutneighsqo = _cutneighsq[0];
      const int ntypes2 = ntypes * ntypes;
      if (_off_threads > 0 && cutneighsqo != nullptr) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cutneighsqo:length(ntypes2) alloc_if(1) free_if(0))
      }
      if (use_ghost_cut) {
        flt_t * cutneighghostsqo = _cutneighghostsq[0];
        if (_off_threads > 0 && cutneighghostsqo != nullptr) {
          #pragma offload_transfer target(mic:_cop) \
            nocopy(cutneighghostsqo:length(ntypes2) alloc_if(1) free_if(0))
        }
      }
      #endif
    }
    _ntypes = ntypes;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
double IntelBuffers<flt_t, acc_t>::memory_usage(const int nthreads)
{
  double tmem = sizeof(atom_t);
  if (lmp->atom->q) tmem += sizeof(flt_t);
  if (lmp->atom->torque) tmem += sizeof(quat_t);
  #ifdef _LMP_INTEL_OFFLOAD
  if (_separate_buffers) tmem *= 2;
  #endif
  tmem *= _buf_size;

  const int fstride = get_stride(_buf_local_size);
  tmem += fstride * nthreads * sizeof(vec3_acc_t);
  #ifdef _LMP_INTEL_OFFLOAD
  if (_off_f) tmem += fstride*_off_threads * sizeof(vec3_acc_t);
  #endif

  tmem += (_list_alloc_atoms + _off_threads) * get_max_nbors() * sizeof(int);
  tmem += _ntypes * _ntypes * sizeof(int);

  tmem += _buf_local_size + (_n_list_ptrs - 1) * _buf_local_size * 2;

  return tmem;
}

/* ---------------------------------------------------------------------- */

template class LAMMPS_NS::IntelBuffers<float,float>;
template class LAMMPS_NS::IntelBuffers<float,double>;
template class LAMMPS_NS::IntelBuffers<double,double>;
