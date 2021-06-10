// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

// loop counters for doing a pack/unpack

#include "pack.h"

/* ----------------------------------------------------------------------
   Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation
     and # of values/element
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   pack from data -> buf
------------------------------------------------------------------------- */

#include "fftdata_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PackKokkos {
 public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

struct pack_3d_functor {
public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  pack_3d_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int in_orig = slow*nmid*nfast;
    const int plane = slow*nstride_plane;
    const int out_orig = plane + mid*nstride_line;
    const int in = in_orig + mid*nfast + fast;
    const int out = out_orig + fast;
    d_buf[buf_offset + in] = d_data[data_offset + out];
  }
};

static void pack_3d(typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  pack_3d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

struct unpack_3d_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast;
    const int plane = slow*nstride_plane;
    const int in_orig = plane + mid*nstride_line;
    const int in = in_orig + fast;
    const int out = out_orig + mid*nfast + fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
  }
};

static void unpack_3d(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
------------------------------------------------------------------------- */


struct unpack_3d_permute1_1_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute1_1_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast;
    const int plane = slow*nstride_line;
    const int in_orig = plane + mid;
    const int in = in_orig + fast*nstride_plane;
    const int out = out_orig + mid*nfast + fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
  }
};

static void unpack_3d_permute1_1(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute1_1_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}
/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 2 values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute1_2_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute1_2_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast*2;
    const int plane = slow*nstride_line;
    const int in_orig = plane + 2*mid;
    const int in = in_orig + fast*nstride_plane;
    const int out = out_orig + 2*mid*nfast + 2*fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
    d_data[data_offset + in+1] = d_buf[buf_offset + out+1];
  }
};

static void unpack_3d_permute1_2(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute1_2_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, nqty values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute1_n_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element

  unpack_3d_permute1_n_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
      nqty          = plan->nqty         ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast*nqty;
    const int plane = slow*nstride_line;
    const int instart = plane + nqty*mid;
    int in = instart + nqty*fast*nstride_plane;
    int out = out_orig + nqty*mid*nfast + nqty*fast;
    for (int iqty = 0; iqty < nqty; iqty++) d_data[data_offset + in++] = d_buf[buf_offset + out++];
  }
};

static void unpack_3d_permute1_n(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute1_n_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 1 value/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_1_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute2_1_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast;
    const int in_orig = slow + mid*nstride_plane;
    const int in = in_orig + fast*nstride_line;
    const int out = out_orig + mid*nfast + fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
  }
};

static void unpack_3d_permute2_1(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute2_1_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 2 values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_2_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute2_2_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast*2;
    const int in_orig = 2*slow + mid*nstride_plane;
    const int in = in_orig + fast*nstride_line;
    const int out = out_orig + 2*mid*nfast + 2*fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
    d_data[data_offset + in+1] = d_buf[buf_offset + out+1];
  }
};

static void unpack_3d_permute2_2(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute2_2_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}
/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, nqty values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_n_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element

  unpack_3d_permute2_n_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nslow         = plan->nslow        ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
      nqty          = plan->nqty         ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int fast = index / (nslow * nmid);
    const int index_new = index - (fast * nslow * nmid);
    const int mid = index_new / nslow;
    const int slow = index_new % nslow;
    const int out_orig = slow*nmid*nfast*nqty;
    const int instart = nqty*slow + mid*nstride_plane;
    int in = instart + nqty*fast*nstride_line;
    int out = out_orig + nqty*mid*nfast + nqty*fast;
    for (int iqty = 0; iqty < nqty; iqty++) d_data[data_offset + in++] = d_buf[buf_offset + out++];
  }
};

static void unpack_3d_permute2_n(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  const int nmid = plan->nmid;
  const int nfast = plan->nfast;
  unpack_3d_permute2_n_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nmid*nfast,f);
  DeviceType().fence();
}

};

}
