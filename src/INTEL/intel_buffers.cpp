// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
    lmp(lmp_in), _x(nullptr), _q(nullptr), _quat(nullptr), _f(nullptr),
    _n_list_ptrs(1), _max_list_ptrs(4), _buf_size(0), _buf_local_size(0) {
  _torque_flag = 0;
  _neigh_list_ptrs = new IntelNeighListPtrs[_max_list_ptrs];
  _neigh_list_ptrs[0].cnumneigh = nullptr;
  _list_alloc_atoms = 0;
  _ntypes = 0;
  _list_local_size = 0;
  _ccachex = nullptr;
  _ncache_alloc = 0;
  _ncachetag = nullptr;
  _cutneighsq = nullptr;
  _cutneighghostsq = nullptr;
  _need_tag = 0;
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
                                       const int nthreads)
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
  if (lmp->force->newton_pair || lmp->atom->molecular)
    lmp->memory->create(_f, f_stride * nthreads, "intel_f");
  else
    lmp->memory->create(_f, f_stride, "intel_f");
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_list_local()
{
  if (_list_local_size > 0) {
    if (_neigh_list_ptrs[0].cnumneigh) {
      int * cnumneigh = _neigh_list_ptrs[0].cnumneigh;
      _neigh_list_ptrs[0].cnumneigh = nullptr;
      lmp->memory->destroy(cnumneigh);
    }
    _list_local_size = 0;
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
void IntelBuffers<flt_t, acc_t>::_grow_list_local(NeighList *list, const int three_body)
{
  free_list_local();
  int size = list->get_maxlocal();
  _list_local_size = size;
  if (three_body)
    lmp->memory->create(_neigh_list_ptrs[0].cnumneigh, size, "_cnumneigh");
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_nbor_list()
{
  if (_list_alloc_atoms > 0) {
    lmp->memory->destroy(_list_alloc);
    _list_alloc_atoms = 0;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_nbor_list(NeighList * /*list*/,
                                                 const int nlocal,
                                                 const int nthreads,
                                                 const int pack_width)
{
  free_nbor_list();
  _list_alloc_atoms = 1.10 * nlocal;
  int nt = nthreads;

  bigint list_alloc_size =
    (bigint)(_list_alloc_atoms + nt * 2 + pack_width - 1) * (bigint)get_max_nbors();
  _list_alloc = (int *) lmp->memory->smalloc(list_alloc_size * sizeof(int), "_list_alloc");
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
void IntelBuffers<flt_t, acc_t>::grow_ccache(const int /*off_flag*/,
        const int nthreads,
        const int width)
{
  if (_ccachex)
    return;

  const int nsize = get_max_nbors() * width;
  int esize = MIN(sizeof(int), sizeof(flt_t));
  IP_PRE_get_stride(_ccache_stride, nsize, esize, 0);
  int nt = nthreads;
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
void IntelBuffers<flt_t, acc_t>::grow_ncache(const int /*off_flag*/,
                                             const int nthreads)
{
  const int nsize = get_max_nbors() * 3;
  int esize = MIN(sizeof(int), sizeof(flt_t));
  IP_PRE_get_stride(_ncache_stride, nsize, esize, 0);
  int nt = nthreads;
  const int vsize = _ncache_stride * nt;

  if (_ncache_alloc) {
    if (vsize > _ncache_alloc || (need_tag() && _ncachetag == nullptr))
      free_ncache();
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
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::fdotr_reduce_l5(const int lf, const int lt,
    const int nthreads, const int f_stride, acc_t &ov0, acc_t &ov1,
    acc_t &ov2, acc_t &ov3, acc_t &ov4, acc_t &ov5)
{
  IP_PRE_fdotr_acc_force_l5(lf, lt, 0, nthreads, _f, f_stride, _x, ov0,
                            ov1, ov2, ov3, ov4, ov5);
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::fdotr_reduce(const int nall,
    const int nthreads, const int f_stride, acc_t &ov0, acc_t &ov1,
    acc_t &ov2, acc_t &ov3, acc_t &ov4, acc_t &ov5)
{
  int iifrom, iito, tid;
  IP_PRE_fdotr_acc_force(nall, 0, nthreads, _f, f_stride, _x, 2,
                         ov0, ov1, ov2, ov3, ov4, ov5);
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::set_ntypes(const int ntypes,
                                            const int use_ghost_cut)
{
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      lmp->memory->destroy(_cutneighsq);
      if (_cutneighghostsq != nullptr) lmp->memory->destroy(_cutneighghostsq);
    }
    if (ntypes > 0) {
      lmp->memory->create(_cutneighsq, ntypes, ntypes, "_cutneighsq");
      if (use_ghost_cut)
        lmp->memory->create(_cutneighghostsq, ntypes, ntypes,
                            "_cutneighghostsq");
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
  if (_torque_flag) tmem += sizeof(quat_t);
  tmem *= _buf_size;

  const int fstride = get_stride(_buf_local_size);
  tmem += fstride * nthreads * sizeof(vec3_acc_t);

  tmem += (bigint)(_list_alloc_atoms) * (bigint)get_max_nbors() * sizeof(int);
  tmem += _ntypes * _ntypes * sizeof(int);

  tmem += _buf_local_size + (_n_list_ptrs - 1) * _buf_local_size * 2;

  return tmem;
}

/* ---------------------------------------------------------------------- */

template class LAMMPS_NS::IntelBuffers<float,float>;
template class LAMMPS_NS::IntelBuffers<float,double>;
template class LAMMPS_NS::IntelBuffers<double,double>;
