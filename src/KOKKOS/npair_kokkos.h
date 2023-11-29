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

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairKokkos<LMPHostType,0,0,0,0,0> NPairKokkosFullBinHost;
NPairStyle(full/bin/kk/host,
           NPairKokkosFullBinHost,
           NP_BIN | NP_KOKKOS_HOST | NP_FULL | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPDeviceType,0,0,0,0,0> NPairKokkosFullBinDevice;
NPairStyle(full/bin/kk/device,
           NPairKokkosFullBinDevice,
           NP_BIN | NP_KOKKOS_DEVICE | NP_FULL | NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPHostType,0,0,1,0,0> NPairKokkosFullBinGhostHost;
NPairStyle(full/bin/ghost/kk/host,
           NPairKokkosFullBinGhostHost,
           NP_BIN | NP_KOKKOS_HOST | NP_FULL | NP_NEWTON | NP_NEWTOFF | NP_GHOST | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPDeviceType,0,0,1,0,0> NPairKokkosFullBinGhostDevice;
NPairStyle(full/bin/ghost/kk/device,
           NPairKokkosFullBinGhostDevice,
           NP_BIN | NP_KOKKOS_DEVICE | NP_FULL | NP_NEWTON | NP_NEWTOFF | NP_GHOST | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPHostType,1,1,0,0,0> NPairKokkosHalfBinNewtonHost;
NPairStyle(half/bin/newton/kk/host,
           NPairKokkosHalfBinNewtonHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTON | NP_ORTHO);

typedef NPairKokkos<LMPHostType,1,0,0,0,0> NPairKokkosHalfBinNewtoffHost;
NPairStyle(half/bin/newtoff/kk/host,
           NPairKokkosHalfBinNewtoffHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTOFF | NP_ORTHO);

typedef NPairKokkos<LMPDeviceType,1,1,0,0,0> NPairKokkosHalfBinNewtonDevice;
NPairStyle(half/bin/newton/kk/device,
           NPairKokkosHalfBinNewtonDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTON | NP_ORTHO);

typedef NPairKokkos<LMPDeviceType,1,0,0,0,0> NPairKokkosHalfBinNewtoffDevice;
NPairStyle(half/bin/newtoff/kk/device,
           NPairKokkosHalfBinNewtoffDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTOFF | NP_ORTHO);

typedef NPairKokkos<LMPHostType,1,1,0,1,0> NPairKokkosHalfBinNewtonTriHost;
NPairStyle(half/bin/newton/kk/host,
           NPairKokkosHalfBinNewtonTriHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTON | NP_TRI);

typedef NPairKokkos<LMPHostType,1,0,0,1,0> NPairKokkosHalfBinNewtoffTriHost;
NPairStyle(half/bin/newtoff/kk/host,
           NPairKokkosHalfBinNewtoffTriHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTOFF | NP_TRI);

typedef NPairKokkos<LMPDeviceType,1,1,0,1,0> NPairKokkosHalfBinNewtonTriDevice;
NPairStyle(half/bin/newton/kk/device,
           NPairKokkosHalfBinNewtonTriDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTON | NP_TRI);

typedef NPairKokkos<LMPDeviceType,1,0,0,1,0> NPairKokkosHalfBinNewtoffTriDevice;
NPairStyle(half/bin/newtoff/kk/device,
           NPairKokkosHalfBinNewtoffTriDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTOFF | NP_TRI);

typedef NPairKokkos<LMPHostType,1,0,1,0,0> NPairKokkosHalfBinNewtoffGhostHost;
NPairStyle(half/bin/newtoff/ghost/kk/host,
           NPairKokkosHalfBinNewtoffGhostHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTOFF | NP_GHOST | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPDeviceType,1,0,1,0,0> NPairKokkosHalfBinNewtoffGhostDevice;
NPairStyle(half/bin/newtoff/ghost/kk/device,
           NPairKokkosHalfBinNewtoffGhostDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTOFF | NP_GHOST | NP_ORTHO | NP_TRI);

typedef NPairKokkos<LMPHostType,1,1,0,0,1> NPairKokkosHalfBinNewtonSizeHost;
NPairStyle(half/bin/newton/size/kk/host,
           NPairKokkosHalfBinNewtonSizeHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTON | NP_SIZE | NP_ORTHO);

typedef NPairKokkos<LMPHostType,1,0,0,0,1> NPairKokkosHalfBinNewtoffSizeHost;
NPairStyle(half/bin/newtoff/size/kk/host,
           NPairKokkosHalfBinNewtoffSizeHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTOFF | NP_SIZE | NP_ORTHO);

typedef NPairKokkos<LMPDeviceType,1,1,0,0,1> NPairKokkosHalfBinNewtonSizeDevice;
NPairStyle(half/bin/newton/size/kk/device,
           NPairKokkosHalfBinNewtonSizeDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTON | NP_SIZE | NP_ORTHO);

typedef NPairKokkos<LMPDeviceType,1,0,0,0,1> NPairKokkosHalfBinNewtoffSizeDevice;
NPairStyle(half/bin/newtoff/size/kk/device,
           NPairKokkosHalfBinNewtoffSizeDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTOFF | NP_SIZE | NP_ORTHO);

typedef NPairKokkos<LMPHostType,1,1,0,1,1> NPairKokkosHalfBinNewtonSizeTriHost;
NPairStyle(half/bin/newton/size/kk/host,
           NPairKokkosHalfBinNewtonSizeTriHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTON | NP_SIZE | NP_TRI);

typedef NPairKokkos<LMPHostType,1,0,0,1,1> NPairKokkosHalfBinNewtoffSizeTriHost;
NPairStyle(half/bin/newtoff/size/kk/host,
           NPairKokkosHalfBinNewtoffSizeTriHost,
           NP_BIN | NP_KOKKOS_HOST | NP_HALF | NP_NEWTOFF | NP_SIZE | NP_TRI);

typedef NPairKokkos<LMPDeviceType,1,1,0,1,1> NPairKokkosHalfBinNewtonSizeTriDevice;
NPairStyle(half/bin/newton/size/kk/device,
           NPairKokkosHalfBinNewtonSizeTriDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTON | NP_SIZE | NP_TRI);

typedef NPairKokkos<LMPDeviceType,1,0,0,1,1> NPairKokkosHalfBinNewtoffSizeTriDevice;
NPairStyle(half/bin/newtoff/size/kk/device,
           NPairKokkosHalfBinNewtoffSizeTriDevice,
           NP_KOKKOS_DEVICE | NP_HALF | NP_BIN | NP_NEWTOFF | NP_SIZE | NP_TRI);
// clang-format on
#else

// clang-format off
#ifndef LMP_NPAIR_KOKKOS_H
#define LMP_NPAIR_KOKKOS_H

#include "npair.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType, int HALF_NEIGH, int NEWTON, int GHOST, int TRI, int SIZE>
class NPairKokkos : public NPair {
  typedef ArrayTypes<DeviceType> AT;

 public:
  NPairKokkos(class LAMMPS *);
  void copy_neighbor_info() override;
  void copy_bin_info() override;
  void copy_stencil_info() override;
  void build(class NeighList *) override;

 private:
  int newton_pair;
  typename AT::t_int_1d d_scalars;
  HAT::t_int_1d h_scalars;
  typename AT::t_int_scalar d_resize;
  typename AT::t_int_scalar d_new_maxneighs;
  HAT::t_int_scalar h_resize;
  HAT::t_int_scalar h_new_maxneighs;

  // data from Neighbor class

  DAT::tdual_xfloat_2d k_cutneighsq;

  // exclusion data from Neighbor class

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::tdual_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;
  DAT::tdual_int_1d k_ex_mol_intra;

  // data from NBin class

  int atoms_per_bin;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;
  DAT::tdual_int_1d k_atom2bin;

  // data from NStencil class

  int nstencil,last_stencil_old;
  DAT::tdual_int_1d k_stencil;  // # of J neighs for each I
  DAT::tdual_int_1d_3 k_stencilxyz;
};

template<class DeviceType>
class NeighborKokkosExecute
{
  typedef ArrayTypes<DeviceType> AT;

 public:
  NeighListKokkos<DeviceType> neigh_list;

  const double delta;

  // data from Neighbor class

  const typename AT::t_xfloat_2d_randomread cutneighsq;

  // exclusion data from Neighbor class

  const int exclude;

  const int nex_type;
  const typename AT::t_int_1d_const ex1_type,ex2_type;
  const typename AT::t_int_2d_const ex_type;

  const int nex_group;
  const typename AT::t_int_1d_const ex1_bit,ex2_bit;

  const int nex_mol;
  const typename AT::t_int_1d_const ex_mol_group;
  const typename AT::t_int_1d_const ex_mol_bit;
  const typename AT::t_int_1d_const ex_mol_intra;

  // data from NBin class

  const int mbins;
  const typename AT::t_int_1d bincount;
  const typename AT::t_int_1d_const c_bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;
  const typename AT::t_int_1d atom2bin;
  const typename AT::t_int_1d_const c_atom2bin;


  // data from NStencil class

  int nstencil;
  typename AT::t_int_1d d_stencil;  // # of J neighs for each I
  typename AT::t_int_1d_3 d_stencilxyz;

  // data from Atom class

  const typename AT::t_x_array_randomread x;
  const typename AT::t_float_1d radius;
  const typename AT::t_int_1d_const type,mask;
  const typename AT::t_tagint_1d_const molecule;
  const typename AT::t_tagint_1d_const tag;
  const typename AT::t_tagint_2d_const special;
  const typename AT::t_int_2d_const nspecial;
  const int molecular;
  int moltemplate;

  int special_flag[4];

  const int nbinx,nbiny,nbinz;
  const int mbinx,mbiny,mbinz;
  const int mbinxlo,mbinylo,mbinzlo;
  const X_FLOAT bininvx,bininvy,bininvz;
  X_FLOAT bboxhi[3],bboxlo[3];

  const int nlocal,nall,neigh_transpose;

  typename AT::t_int_scalar resize;
  typename AT::t_int_scalar new_maxneighs;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_new_maxneighs;

  const int xperiodic, yperiodic, zperiodic;
  const int xprd_half, yprd_half, zprd_half;

  // GRANULAR required member variables
  const X_FLOAT skin;

  NeighborKokkosExecute(
                        const NeighListKokkos<DeviceType> &_neigh_list,
                        const typename AT::t_xfloat_2d_randomread &_cutneighsq,
                        const typename AT::t_int_1d &_bincount,
                        const typename AT::t_int_2d &_bins,
                        const typename AT::t_int_1d &_atom2bin,
                        const int _mbins,const int _nstencil,
                        const typename AT::t_int_1d &_d_stencil,
                        const typename AT::t_int_1d_3 &_d_stencilxyz,
                        const int _nlocal,const int _nall,const int _neigh_transpose,
                        const typename AT::t_x_array_randomread &_x,
                        const typename AT::t_float_1d &_radius,
                        const typename AT::t_int_1d_const &_type,
                        const typename AT::t_int_1d_const &_mask,
                        const typename AT::t_tagint_1d_const &_molecule,
                        const typename AT::t_tagint_1d_const &_tag,
                        const typename AT::t_tagint_2d_const &_special,
                        const typename AT::t_int_2d_const &_nspecial,
                        const int &_molecular,
                        const int & _nbinx,const int & _nbiny,const int & _nbinz,
                        const int & _mbinx,const int & _mbiny,const int & _mbinz,
                        const int & _mbinxlo,const int & _mbinylo,const int & _mbinzlo,
                        const X_FLOAT &_bininvx,const X_FLOAT &_bininvy,const X_FLOAT &_bininvz,
                        const double &_delta,const int & _exclude,const int & _nex_type,
                        const typename AT::t_int_1d_const & _ex1_type,
                        const typename AT::t_int_1d_const & _ex2_type,
                        const typename AT::t_int_2d_const & _ex_type,
                        const int & _nex_group,
                        const typename AT::t_int_1d_const & _ex1_bit,
                        const typename AT::t_int_1d_const & _ex2_bit,
                        const int & _nex_mol,
                        const typename AT::t_int_1d_const & _ex_mol_group,
                        const typename AT::t_int_1d_const & _ex_mol_bit,
                        const typename AT::t_int_1d_const & _ex_mol_intra,
                        const X_FLOAT *_bboxhi, const X_FLOAT* _bboxlo,
                        const int & _xperiodic, const int & _yperiodic, const int & _zperiodic,
                        const int & _xprd_half, const int & _yprd_half, const int & _zprd_half,
                        const X_FLOAT _skin,
                        const typename AT::t_int_scalar _resize,
                        const typename ArrayTypes<LMPHostType>::t_int_scalar _h_resize,
                        const typename AT::t_int_scalar _new_maxneighs,
                        const typename ArrayTypes<LMPHostType>::t_int_scalar _h_new_maxneighs):
    neigh_list(_neigh_list), cutneighsq(_cutneighsq),delta(_delta),exclude(_exclude),
    nex_type(_nex_type),ex1_type(_ex1_type),ex2_type(_ex2_type),
    ex_type(_ex_type),nex_group(_nex_group),
    ex1_bit(_ex1_bit),ex2_bit(_ex2_bit),
    nex_mol(_nex_mol),ex_mol_group(_ex_mol_group),ex_mol_bit(_ex_mol_bit),
    ex_mol_intra(_ex_mol_intra),mbins(_mbins),
    bincount(_bincount),c_bincount(_bincount),bins(_bins),c_bins(_bins),
    atom2bin(_atom2bin),c_atom2bin(_atom2bin),
    nstencil(_nstencil),d_stencil(_d_stencil),d_stencilxyz(_d_stencilxyz),
    x(_x),radius(_radius),type(_type),mask(_mask),molecule(_molecule),
    tag(_tag),special(_special),nspecial(_nspecial),molecular(_molecular),
    nbinx(_nbinx),nbiny(_nbiny),nbinz(_nbinz),
    mbinx(_mbinx),mbiny(_mbiny),mbinz(_mbinz),
    mbinxlo(_mbinxlo),mbinylo(_mbinylo),mbinzlo(_mbinzlo),
    bininvx(_bininvx),bininvy(_bininvy),bininvz(_bininvz),
    nlocal(_nlocal),nall(_nall),neigh_transpose(_neigh_transpose),
    xperiodic(_xperiodic),yperiodic(_yperiodic),zperiodic(_zperiodic),
    xprd_half(_xprd_half),yprd_half(_yprd_half),zprd_half(_zprd_half),
    skin(_skin),resize(_resize),h_resize(_h_resize),
    new_maxneighs(_new_maxneighs),h_new_maxneighs(_h_new_maxneighs) {

    if (molecular == 2) moltemplate = 1;
    else moltemplate = 0;

    bboxlo[0] = _bboxlo[0]; bboxlo[1] = _bboxlo[1]; bboxlo[2] = _bboxlo[2];
    bboxhi[0] = _bboxhi[0]; bboxhi[1] = _bboxhi[1]; bboxhi[2] = _bboxhi[2];

    h_resize() = 1;
    h_new_maxneighs() = neigh_list.maxneighs;
  };

  ~NeighborKokkosExecute() {neigh_list.copymode = 1;};

  template<int HalfNeigh, int Newton, int Tri>
  KOKKOS_FUNCTION
  void build_Item(const int &i) const;

  template<int HalfNeigh>
  KOKKOS_FUNCTION
  void build_ItemGhost(const int &i) const;

  template<int HalfNeigh, int Newton, int Tri>
  KOKKOS_FUNCTION
  void build_ItemSize(const int &i) const;

#ifdef LMP_KOKKOS_GPU
  template<int HalfNeigh, int Newton, int Tri>
  LAMMPS_DEVICE_FUNCTION inline
  void build_ItemGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                     size_t sharedsize) const;

  template<int HalfNeigh>
  LAMMPS_DEVICE_FUNCTION inline
  void build_ItemGhostGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                     size_t sharedsize) const;

  template<int HalfNeigh, int Newton, int Tri>
  LAMMPS_DEVICE_FUNCTION inline
  void build_ItemSizeGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                         size_t sharedsize) const;
#endif

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z, int* i) const
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

    i[0] = ix - mbinxlo;
    i[1] = iy - mbinylo;
    i[2] = iz - mbinzlo;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

  KOKKOS_INLINE_FUNCTION
  int exclusion(const int &i,const int &j, const int &itype,const int &jtype) const;

  KOKKOS_INLINE_FUNCTION
  int find_special(const int &i, const int &j) const;

  KOKKOS_INLINE_FUNCTION
  int minimum_image_check(double dx, double dy, double dz) const {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

};

template<class DeviceType, int HALF_NEIGH, int NEWTON, int TRI>
struct NPairKokkosBuildFunctor {
  typedef DeviceType device_type;

  const NeighborKokkosExecute<DeviceType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctor(const NeighborKokkosExecute<DeviceType> &_c,
                             size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_Item<HALF_NEIGH,NEWTON,TRI>(i);
  }
#ifdef LMP_KOKKOS_GPU
  LAMMPS_DEVICE_FUNCTION inline
  void operator() (typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    c.template build_ItemGPU<HALF_NEIGH,NEWTON,TRI>(dev, sharedsize);
  }
  size_t team_shmem_size(const int team_size) const { (void) team_size; return sharedsize; }
#endif
};

template<int HALF_NEIGH, int NEWTON, int TRI>
struct NPairKokkosBuildFunctor<LMPHostType,HALF_NEIGH,NEWTON,TRI> {
  typedef LMPHostType device_type;

  const NeighborKokkosExecute<LMPHostType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctor(const NeighborKokkosExecute<LMPHostType> &_c,
                             size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_Item<HALF_NEIGH,NEWTON,TRI>(i);
  }

  void operator() (typename Kokkos::TeamPolicy<LMPHostType>::member_type /*dev*/) const {} // Should error out
};

template<class DeviceType,int HALF_NEIGH>
struct NPairKokkosBuildFunctorGhost {
  typedef DeviceType device_type;

  const NeighborKokkosExecute<DeviceType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctorGhost(const NeighborKokkosExecute<DeviceType> &_c,
                             size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_ItemGhost<HALF_NEIGH>(i);
  }

#ifdef LMP_KOKKOS_GPU
  LAMMPS_DEVICE_FUNCTION inline
  void operator() (typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    c.template build_ItemGhostGPU<HALF_NEIGH>(dev, sharedsize);
  }
  size_t team_shmem_size(const int team_size) const { (void) team_size; return sharedsize; }
#endif
};

template<int HALF_NEIGH>
struct NPairKokkosBuildFunctorGhost<LMPHostType,HALF_NEIGH> {
  typedef LMPHostType device_type;

  const NeighborKokkosExecute<LMPHostType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctorGhost(const NeighborKokkosExecute<LMPHostType> &_c,
                             size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_ItemGhost<HALF_NEIGH>(i);
  }

  void operator() (typename Kokkos::TeamPolicy<LMPHostType>::member_type /*dev*/) const {} // Should error out
};

template <class DeviceType, int HALF_NEIGH, int NEWTON, int TRI>
struct NPairKokkosBuildFunctorSize {
  typedef DeviceType device_type;

  const NeighborKokkosExecute<DeviceType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctorSize(const NeighborKokkosExecute<DeviceType> &_c,
                              size_t _sharedsize): c(_c), sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_ItemSize<HALF_NEIGH,NEWTON,TRI>(i);
  }

#ifdef LMP_KOKKOS_GPU
  LAMMPS_DEVICE_FUNCTION inline
  void operator() (typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    c.template build_ItemSizeGPU<HALF_NEIGH,NEWTON,TRI>(dev, sharedsize);
  }
  size_t team_shmem_size(const int team_size) const { (void) team_size; return sharedsize; }
#endif
};

template <int HALF_NEIGH, int NEWTON, int TRI>
struct NPairKokkosBuildFunctorSize<LMPHostType,HALF_NEIGH,NEWTON,TRI> {
  typedef LMPHostType device_type;

  const NeighborKokkosExecute<LMPHostType> c;
  size_t sharedsize;

  NPairKokkosBuildFunctorSize(const NeighborKokkosExecute<LMPHostType> &_c,
                              size_t _sharedsize): c(_c), sharedsize(_sharedsize) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_ItemSize<HALF_NEIGH,NEWTON,TRI>(i);
  }

  void operator() (typename Kokkos::TeamPolicy<LMPHostType>::member_type /*dev*/) const {} // Should error out
};

}

#endif
#endif

