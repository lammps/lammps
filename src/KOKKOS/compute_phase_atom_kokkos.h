/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(phase/atom/kk,ComputePhaseAtomKokkos<LMPDeviceType>)
ComputeStyle(phase/atom/kk/device,ComputePhaseAtomKokkos<LMPDeviceType>)
ComputeStyle(phase/atom/kk/host,ComputePhaseAtomKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_PHASE_KOKKOS_ATOM_H
#define LMP_COMPUTE_PHASE_KOKKOS_ATOM_H

#include "compute_phase_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputePhaseAtom{};

template<class DeviceType>
class ComputePhaseAtomKokkos : public ComputePhaseAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputePhaseAtomKokkos(class LAMMPS *, int, char **);
  virtual ~ComputePhaseAtomKokkos();
  void init();
  void compute_peratom();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePhaseAtom, const int&) const;

 private:
  typename AT::t_x_array_randomread x;
  typename AT::t_v_array_randomread v;
  typename ArrayTypes<DeviceType>::t_int_1d mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_float_2d k_phase;
  typename AT::t_float_2d d_phase;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
