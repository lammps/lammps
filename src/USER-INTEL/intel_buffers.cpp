/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
IntelBuffers<flt_t, acc_t>::IntelBuffers(class LAMMPS *lmp_in) :
    lmp(lmp_in), _x(0), _q(0), _quat(0), _f(0), _off_threads(0),
    _buf_size(0), _buf_local_size(0) {
  _list_alloc_atoms = 0;
  _ntypes = 0;
  _off_map_maxlocal = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  _separate_buffers = 0;
  _off_f = 0;
  _off_map_ilist = 0;
  _off_map_nmax = 0;
  _off_map_maxhead = 0;
  _off_list_alloc = false;
  _off_threads = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
IntelBuffers<flt_t, acc_t>::~IntelBuffers()
{
  free_buffers();
  free_all_nbor_buffers();
  set_ntypes(0);
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
    if (q != 0) lmp->memory->destroy(q);
    if (quat != 0) lmp->memory->destroy(quat);
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
  if (lmp->atom->torque)
    _buf_local_size *= 2;
  const int f_stride = get_stride(_buf_local_size);
  lmp->memory->create(_x, _buf_size,"intel_x");
  if (lmp->atom->q != NULL)
    lmp->memory->create(_q, _buf_size, "intel_q");
  if (lmp->atom->ellipsoid != NULL)
    lmp->memory->create(_quat, _buf_size, "intel_quat");
  lmp->memory->create(_f, f_stride * nthreads, "intel_f");

  #ifdef _LMP_INTEL_OFFLOAD
  if (_separate_buffers) {
    lmp->memory->create(_host_x, _buf_size,"intel_host_x");
    if (lmp->atom->q != NULL)
      lmp->memory->create(_host_q, _buf_size, "intel_host_q");
    if (lmp->atom->ellipsoid != NULL)
      lmp->memory->create(_host_quat, _buf_size, "intel_host_quat");
  }
    
  if (offload_end > 0) {
    lmp->memory->create(_off_f, f_stride * _off_threads, "intel_off_f");
    const atom_t * const x = get_x();
    const flt_t * const q = get_q();
    const vec3_acc_t * f_start = get_off_f();
    acc_t * ev_global = get_ev_global();
    if (lmp->atom->q != NULL) {
      if (x != NULL && q != NULL && f_start != NULL && ev_global != NULL) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(x,q:length(_buf_size) alloc_if(1) free_if(0)) \
	  nocopy(f_start:length(f_stride*_off_threads) alloc_if(1) free_if(0))\
	  nocopy(ev_global:length(8) alloc_if(1) free_if(0))
      }
    } else {
      if (x != NULL && f_start != NULL && ev_global != NULL) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(x:length(_buf_size) alloc_if(1) free_if(0)) \
          nocopy(f_start:length(f_stride*_off_threads) alloc_if(1) free_if(0))\
	  nocopy(ev_global:length(8) alloc_if(1) free_if(0))
      }
    }
    if (lmp->atom->ellipsoid != NULL) {
      const quat_t * const quat = get_quat();
      if (quat != NULL) {
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
    const int * tag = _off_map_tag;
    const int * special = _off_map_special;
    const int * nspecial = _off_map_nspecial;
    const int * bins = _off_map_bins;
    if (tag != 0 && special != 0 && nspecial !=0 && bins != 0) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(tag:alloc_if(0) free_if(1)) \
	nocopy(special,nspecial:alloc_if(0) free_if(1)) \
	nocopy(bins:alloc_if(0) free_if(1))
    }
    _off_map_nmax = 0;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_nmax()
{
  #ifdef _LMP_INTEL_OFFLOAD
  free_nmax();
  int *special, *nspecial;
  int tag_length, special_length, nspecial_length;
  int size = lmp->atom->nmax;
  if (lmp->atom->molecular) {
    special = lmp->atom->special[0];
    nspecial = lmp->atom->nspecial[0];
    special_length = size * lmp->atom->maxspecial;
    nspecial_length = size * 3;
    tag_length = size;
  } else {
    special = &_special_holder;
    nspecial = &_nspecial_holder;
    special_length = 1;
    nspecial_length = 1;
    tag_length = 1;
  }
  int *tag = lmp->atom->tag;
  int *bins = lmp->neighbor->bins;
  #pragma offload_transfer target(mic:_cop) \
    nocopy(bins:length(size) alloc_if(1) free_if(0)) \
    nocopy(tag:length(tag_length) alloc_if(1) free_if(0)) \
    nocopy(special:length(special_length) alloc_if(1) free_if(0)) \
    nocopy(nspecial:length(nspecial_length) alloc_if(1) free_if(0))
  _off_map_tag = tag;
  _off_map_special = special;
  _off_map_nspecial = nspecial;
  _off_map_nmax = size;
  _off_map_bins = bins;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_local()
{
  if (_off_map_maxlocal > 0) {
    int * cnumneigh = _cnumneigh;
    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_map_ilist != NULL) {
      const int * ilist = _off_map_ilist;
      const int * numneigh = _off_map_numneigh;
      _off_map_ilist = NULL;
      if (numneigh != 0 && ilist != 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(ilist,numneigh,cnumneigh:alloc_if(0) free_if(1))
      }
    }
    #endif
    lmp->memory->destroy(cnumneigh);
    _off_map_maxlocal = 0;
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_local(NeighList *list, 
					     const int offload_end)
{
  free_local();
  int size = list->get_maxlocal();
  lmp->memory->create(_cnumneigh, size, "_cnumneigh");
  _off_map_maxlocal = size;

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_end > 0) {
    int * numneigh = list->numneigh;
    int * ilist = list->ilist;
    int * cnumneigh = _cnumneigh;
    if (cnumneigh != 0) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(ilist:length(size) alloc_if(1) free_if(0)) \
	nocopy(numneigh:length(size) alloc_if(1) free_if(0)) \
	nocopy(cnumneigh:length(size) alloc_if(1) free_if(0))
    }
    _off_map_ilist = ilist;
    _off_map_numneigh = numneigh;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_binhead()
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_off_map_maxhead > 0) {
    const int * binhead = _off_map_binhead;
    if (binhead !=0) {
      #pragma offload_transfer target(mic:_cop) \
        nocopy(binhead:alloc_if(0) free_if(1))
    }
    _off_map_maxhead = 0;
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_binhead()
{
  #ifdef _LMP_INTEL_OFFLOAD
  free_binhead();
  int * binhead = lmp->neighbor->binhead;
  const int maxhead = lmp->neighbor->maxhead;
  #pragma offload_transfer target(mic:_cop) \
    nocopy(binhead:length(maxhead) alloc_if(1) free_if(0))
  _off_map_binhead = binhead;
  _off_map_maxhead = maxhead;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::free_nbor_list()
{
  if (_list_alloc_atoms > 0) {
    lmp->memory->destroy(_list_alloc);
    _list_alloc_atoms = 0;

    #ifdef _LMP_INTEL_OFFLOAD
    if (_off_list_alloc) {
      int * list_alloc = _list_alloc;
      int * special_flag = lmp->neighbor->special_flag_alloc();
      int * stencil = _off_map_stencil;
      if (list_alloc != 0 && special_flag != 0 && stencil != 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(special_flag,stencil:alloc_if(0) free_if(1)) \
          nocopy(list_alloc:alloc_if(0) free_if(1))
      }
      _off_list_alloc = false;
    }
    #endif
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_nbor_list(NeighList *list, 
						 const int nlocal,
						 const int offload_end)
{
  free_nbor_list();
  _list_alloc_atoms = 1.10 * nlocal;
  int list_alloc_size = (_list_alloc_atoms + _off_threads) * get_max_nbors();
  lmp->memory->create(_list_alloc, list_alloc_size, "_list_alloc");
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_end > 0) {
    int * list_alloc =_list_alloc;
    int * special_flag = lmp->neighbor->special_flag;
    int * stencil = list->stencil;

    if (special_flag != NULL && list_alloc != NULL) {
      #pragma offload_transfer target(mic:_cop) \
        in(special_flag:length(4) alloc_if(1) free_if(0)) \
	in(stencil:length(list->maxstencil) alloc_if(1) free_if(0)) \
	nocopy(list_alloc:length(list_alloc_size) alloc_if(1) free_if(0))
      _off_map_stencil = stencil;
      _off_list_alloc = true;
    }
  }
  #endif
}

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::_grow_stencil(NeighList *list)
{
  #ifdef _LMP_INTEL_OFFLOAD
  int * stencil = _off_map_stencil;
  #pragma offload_transfer target(mic:_cop) \
    nocopy(stencil:alloc_if(0) free_if(1))
  stencil = list->stencil;
  #pragma offload_transfer target(mic:_cop) \
    in(stencil:length(list->maxstencil) alloc_if(1) free_if(0))
  _off_map_stencil = stencil;
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void IntelBuffers<flt_t, acc_t>::set_ntypes(const int ntypes)
{
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * cutneighsqo = _cutneighsq[0];
      if (cutneighsqo != 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cutneighsqo:alloc_if(0) free_if(1))
      }
      #endif
      lmp->memory->destroy(_cutneighsq);
    }
    if (ntypes > 0) {
      lmp->memory->create(_cutneighsq, ntypes, ntypes, "_cutneighsq");
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * cutneighsqo = _cutneighsq[0];
      if (cutneighsqo != NULL) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(cutneighsqo:length(ntypes * ntypes) alloc_if(1) free_if(0))
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

  tmem += _off_map_maxlocal * sizeof(int);
  tmem += (_list_alloc_atoms + _off_threads) * get_max_nbors() * sizeof(int);
  tmem += _ntypes * _ntypes * sizeof(int);

  return tmem;
}

/* ---------------------------------------------------------------------- */

template class IntelBuffers<float,float>;
template class IntelBuffers<float,double>;
template class IntelBuffers<double,double>;
