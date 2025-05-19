/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(SPECIES/ATOM/kk,ComputeSpeciesAtomKokkos<LMPDeviceType>);
ComputeStyle(SPECIES/ATOM/kk/device,ComputeSpeciesAtomKokkos<LMPDeviceType>);
ComputeStyle(SPECIES/ATOM/kk/host,ComputeSpeciesAtomKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPECIES_ATOM_KOKKOS_H
#define LMP_COMPUTE_SPECIES_ATOM_KOKKOS_H

#include "compute_species_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeSpeciesAtomKokkos : public ComputeSpeciesAtom {
 public:
  ComputeSpeciesAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeSpeciesAtomKokkos() override;
  
  void compute_peratom() override;
  
 private:

  DAT::tdual_double_1d k_vector;   // DualView for vector data
  DAT::tdual_double_2d k_array;    // DualView for array data
  DAT::tdual_int_1d k_properties;  // DualView for property types mapping

  // ReaxFF data accessor
  typename ArrayTypes<DeviceType>::t_float_2d d_tmpbo;

};

}    // namespace LAMMPS_NS

#endif
#endif
