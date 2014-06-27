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

#ifndef LMP_NEIGHBOR_KOKKOS_H
#define LMP_NEIGHBOR_KOKKOS_H

#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class Device>
class NeighborKokkosExecute
{
  typedef ArrayTypes<Device> AT;

 public:
  NeighListKokkos<Device> neigh_list;
  const typename AT::t_xfloat_2d_randomread cutneighsq;
  const typename AT::t_int_1d bincount;
  const typename AT::t_int_1d_const c_bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;
  const typename AT::t_x_array_randomread x;
  const typename AT::t_int_1d_const type,mask;
  const typename AT::t_tagint_1d_const molecule;

  const int nbinx,nbiny,nbinz;
  const int mbinx,mbiny,mbinz;
  const int mbinxlo,mbinylo,mbinzlo;
  const X_FLOAT bininvx,bininvy,bininvz;
  X_FLOAT bboxhi[3],bboxlo[3];

  const int nlocal;

  typename AT::t_int_scalar resize;
  typename AT::t_int_scalar new_maxneighs;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_new_maxneighs;

  NeighborKokkosExecute(
    const NeighListKokkos<Device> &_neigh_list,
    const typename AT::t_xfloat_2d_randomread &_cutneighsq,
    const typename AT::t_int_1d &_bincount,
    const typename AT::t_int_2d &_bins,
    const int _nlocal,
        const typename AT::t_x_array_randomread &_x,
    const typename AT::t_int_1d_const &_type,
    const typename AT::t_int_1d_const &_mask,
    const typename AT::t_tagint_1d_const &_molecule,
    const int & _nbinx,const int & _nbiny,const int & _nbinz,
    const int & _mbinx,const int & _mbiny,const int & _mbinz,
    const int & _mbinxlo,const int & _mbinylo,const int & _mbinzlo,
    const X_FLOAT &_bininvx,const X_FLOAT &_bininvy,const X_FLOAT &_bininvz,
    const X_FLOAT *_bboxhi, const X_FLOAT* _bboxlo):
    neigh_list(_neigh_list), cutneighsq(_cutneighsq),
    bincount(_bincount),c_bincount(_bincount),bins(_bins),c_bins(_bins),
    nlocal(_nlocal),
    x(_x),type(_type),mask(_mask),molecule(_molecule),
    nbinx(_nbinx),nbiny(_nbiny),nbinz(_nbinz),
    mbinx(_mbinx),mbiny(_mbiny),mbinz(_mbinz),
    mbinxlo(_mbinxlo),mbinylo(_mbinylo),mbinzlo(_mbinzlo),
    bininvx(_bininvx),bininvy(_bininvy),bininvz(_bininvz) {

    bboxlo[0] = _bboxlo[0]; bboxlo[1] = _bboxlo[1]; bboxlo[2] = _bboxlo[2];
    bboxhi[0] = _bboxhi[0]; bboxhi[1] = _bboxhi[1]; bboxhi[2] = _bboxhi[2];
    
    resize = typename AT::t_int_scalar("NeighborKokkosFunctor::resize");
#ifndef KOKKOS_USE_UVM
    h_resize = Kokkos::create_mirror_view(resize);
#else
    h_resize = resize;
#endif
    h_resize() = 1;
    new_maxneighs = typename AT::
      t_int_scalar("NeighborKokkosFunctor::new_maxneighs");
#ifndef KOKKOS_USE_UVM
    h_new_maxneighs = Kokkos::create_mirror_view(new_maxneighs);
#else
    h_new_maxneighs = new_maxneighs;
#endif
    h_new_maxneighs() = neigh_list.maxneighs;
  };

  ~NeighborKokkosExecute() {neigh_list.clean_copy();};

  template<int HalfNeigh, int GhostNewton>
  KOKKOS_FUNCTION
  void build_Item(const int &i) const;

  template<int ClusterSize>
  KOKKOS_FUNCTION
  void build_cluster_Item(const int &i) const;

#if DEVICE==2
  template<int HalfNeigh>
  __device__ inline
  void build_ItemCuda(Device dev) const;
#endif

  KOKKOS_INLINE_FUNCTION
  void binatomsItem(const int &i) const;

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z) const
  {
    int ix,iy,iz;

    if (x >= bboxhi[0])
      ix = static_cast<int> ((x-bboxhi[0])*bininvx) + nbinx;
    else if (x >= bboxlo[0]) {
      ix = static_cast<int> ((x-bboxlo[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo[0])*bininvx) - 1;

    if (y >= bboxhi[1])
      iy = static_cast<int> ((y-bboxhi[1])*bininvy) + nbiny;
    else if (y >= bboxlo[1]) {
      iy = static_cast<int> ((y-bboxlo[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo[1])*bininvy) - 1;

    if (z >= bboxhi[2])
      iz = static_cast<int> ((z-bboxhi[2])*bininvz) + nbinz;
    else if (z >= bboxlo[2]) {
      iz = static_cast<int> ((z-bboxlo[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo[2])*bininvz) - 1;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }
};

template<class Device>
struct NeighborKokkosBinAtomsFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;

  NeighborKokkosBinAtomsFunctor(const NeighborKokkosExecute<Device> &_c):
    c(_c) {};
  ~NeighborKokkosBinAtomsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.binatomsItem(i);
  }
};

template<class Device,int HALF_NEIGH,int GHOST_NEWTON>
struct NeighborKokkosBuildFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;
  const size_t sharedsize;

  NeighborKokkosBuildFunctor(const NeighborKokkosExecute<Device> &_c, 
                             const size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {};
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_Item<HALF_NEIGH,GHOST_NEWTON>(i);
  }
#if DEVICE==2
  KOKKOS_INLINE_FUNCTION
  void operator() (Device dev) const {
    c.template build_ItemCuda<HALF_NEIGH>(dev);
  }
  size_t shmem_size() const { return sharedsize; }
#endif
};

template<class Device,int ClusterSize>
struct NeighborClusterKokkosBuildFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;
  const size_t sharedsize;

  NeighborClusterKokkosBuildFunctor(const NeighborKokkosExecute<Device> &_c,
                             const size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_cluster_Item<ClusterSize>(i);
  }
};

class NeighborKokkos : public Neighbor {
 public:
  class AtomKokkos *atomKK;

  int nlist_host;                       // pairwise neighbor lists on Host
  NeighListKokkos<LMPHostType> **lists_host;
  int nlist_device;                     // pairwise neighbor lists on Device
  NeighListKokkos<LMPDeviceType> **lists_device;

  NeighborKokkos(class LAMMPS *);
  ~NeighborKokkos();
  void init();

 private:
  int atoms_per_bin;
  DAT::tdual_xfloat_2d k_cutneighsq;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;

  void init_cutneighsq_kokkos(int);
  int init_lists_kokkos();
  void init_list_flags1_kokkos(int);
  void init_list_flags2_kokkos(int);
  void init_list_grow_kokkos(int);
  void choose_build(int, NeighRequest *);
  void build_kokkos(int);
  void setup_bins_kokkos(int);
  
  typedef void (NeighborKokkos::*PairPtrHost)
    (class NeighListKokkos<LMPHostType> *);
  PairPtrHost *pair_build_host;
  typedef void (NeighborKokkos::*PairPtrDevice)
    (class NeighListKokkos<LMPDeviceType> *);
  PairPtrDevice *pair_build_device;

  template<class DeviceType,int HALF_NEIGH>
  void full_bin_kokkos(NeighListKokkos<DeviceType> *list);
  template<class DeviceType>
  void full_bin_cluster_kokkos(NeighListKokkos<DeviceType> *list);

  typedef void (NeighborKokkos::*StencilPtrHost)
    (class NeighListKokkos<LMPHostType> *, int, int, int);
  StencilPtrHost *stencil_create_host;
  typedef void (NeighborKokkos::*StencilPtrDevice)
    (class NeighListKokkos<LMPDeviceType> *, int, int, int);
  StencilPtrDevice *stencil_create_device;
};

}

#endif

/* ERROR/WARNING messages:

*/
