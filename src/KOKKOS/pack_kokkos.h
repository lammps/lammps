/* ----------------------------------------------------------------------
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

#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PackKokkos {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

struct pack_3d_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  pack_3d_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int in = slow*nmid*nfast;
    const int plane = slow*nstride_plane;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int out = plane + mid*nstride_line;
      for (int fast = 0; fast < nfast; fast++)
        d_buf[buf_offset + in++] = d_data[data_offset + out++];
    }
  }
};

static void pack_3d(typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  pack_3d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

struct unpack_3d_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast;
    const int plane = slow*nstride_plane;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int in = plane + mid*nstride_line;
      for (int fast = 0; fast < nfast; fast++)
        d_data[data_offset + in++] = d_buf[buf_offset + out++];
    }
  }
};

static void unpack_3d(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
------------------------------------------------------------------------- */


struct unpack_3d_permute1_1_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute1_1_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast;
    const int plane = slow*nstride_line;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int in = plane + mid;
      for (int fast = 0; fast < nfast; fast++, in += nstride_plane)
        d_data[data_offset + in] = d_buf[buf_offset + out++];
    }
  }
};

static void unpack_3d_permute1_1(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute1_1_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}
/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 2 values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute1_2_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute1_2_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast*2;
    const int plane = slow*nstride_line;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int in = plane + 2*mid;
      for (int fast = 0; fast < nfast; fast++, in += nstride_plane) {
        d_data[data_offset + in] = d_buf[buf_offset + out++];
        d_data[data_offset + in+1] = d_buf[buf_offset + out++];
      }
    }
  }
};

static void unpack_3d_permute1_2(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute1_2_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, nqty values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute1_n_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element

  unpack_3d_permute1_n_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
      nqty          = plan->nqty         ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast*nqty;
    const int plane = slow*nstride_line;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int instart = plane + nqty*mid;
      for (int fast = 0; fast < nfast; fast++, instart += nstride_plane) {
        int in = instart;
        for (int iqty = 0; iqty < nqty; iqty++) d_data[data_offset + in++] = d_buf[buf_offset + out++];
      }
    }
  }
};

static void unpack_3d_permute1_n(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute1_n_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 1 value/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_1_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute2_1_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int in = slow + mid*nstride_plane;
      for (int fast = 0; fast < nfast; fast++, in += nstride_line)
        d_data[data_offset + in] = d_buf[buf_offset + out++];
    }
  }
};

static void unpack_3d_permute2_1(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute2_1_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 2 values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_2_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices

  unpack_3d_permute2_2_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast*2;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int in = 2*slow + mid*nstride_plane;
      for (int fast = 0; fast < nfast; fast++, in += nstride_line) {
        d_data[data_offset + in] = d_buf[buf_offset + out++];
        d_data[data_offset + in+1] = d_buf[buf_offset + out++];
      }
    }
  }
};

static void unpack_3d_permute2_2(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute2_2_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}
/* ----------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, nqty values/element
------------------------------------------------------------------------- */

struct unpack_3d_permute2_n_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element

  unpack_3d_permute2_n_functor(typename AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_3d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nmid          = plan->nmid         ;
      nstride_line  = plan->nstride_line ;
      nstride_plane = plan->nstride_plane;
      nqty          = plan->nqty         ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &slow) const {
    int out = slow*nmid*nfast*nqty;
    for (int mid = 0; mid < nmid; mid++) { // TODO: use thread team or flatten loop
      int instart = nqty*slow + mid*nstride_plane;
      for (int fast = 0; fast < nfast; fast++, instart += nstride_line) {
        int in = instart;
        for (int iqty = 0; iqty < nqty; iqty++) d_data[data_offset + in++] = d_buf[buf_offset + out++];
      }
    }
  }
};

static void unpack_3d_permute2_n(typename AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_3d *plan)
{
  const int nslow = plan->nslow;
  unpack_3d_permute2_n_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow,f);
  DeviceType::fence();
}

};

}
