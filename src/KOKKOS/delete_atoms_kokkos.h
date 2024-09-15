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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(delete_atoms/kk,DeleteAtomsKokkos<LMPDeviceType>);
CommandStyle(delete_atoms/kk/device,DeleteAtomsKokkos<LMPDeviceType>);
CommandStyle(delete_atoms/kk/host,DeleteAtomsKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_DELETE_ATOMS_KOKKOS_H
#define LMP_DELETE_ATOMS_KOKKOS_H

#include "delete_atoms.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class DeleteAtomsKokkos : public DeleteAtoms {
 public:
  DeleteAtomsKokkos(class LAMMPS *);

  void command(int, char **) override;

  void delete_overlap(int, char **);

  //KOKKOS_INLINE_FUNCTION
  //void operator()(const int &i) const;

 protected:

  DAT::tdual_int_1d k_dlist;

};
}    // namespace LAMMPS_NS

#endif    //LMP_DELETE_ATOMS_KOKKOS_H
#endif
