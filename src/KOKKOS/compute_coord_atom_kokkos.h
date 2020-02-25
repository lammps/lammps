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

#ifdef COMPUTE_CLASS

ComputeStyle(coord/atom/kk,ComputeCoordAtomKokkos<LMPDeviceType>)
ComputeStyle(coord/atom/kk/device,ComputeCoordAtomKokkos<LMPDeviceType>)
ComputeStyle(coord/atom/kk/host,ComputeCoordAtomKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_COORD_ATOM_KOKKOS_H
#define LMP_COMPUTE_COORD_ATOM_KOKKOS_H

#include "compute_coord_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int CSTYLE, int NCOL>
struct TagComputeCoordAtom{};

template<class DeviceType>
class ComputeCoordAtomKokkos : public ComputeCoordAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeCoordAtomKokkos(class LAMMPS *, int, char **);
  virtual ~ComputeCoordAtomKokkos();
  void init();
  void compute_peratom();
  enum {NONE,CUTOFF,ORIENT};

  template<int CSTYLE, int NCOL>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeCoordAtom<CSTYLE,NCOL>, const int&) const;

 private:
  int inum;

  typename AT::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_1d d_typelo;
  typename AT::t_int_1d d_typehi;

  DAT::tdual_float_1d k_cvec;
  typename AT::t_float_1d d_cvec;
  DAT::tdual_float_2d k_carray;
  typename AT::t_float_2d d_carray;

  typename AT::t_float_2d d_normv;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
