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

#ifdef FIX_CLASS
// clang-format off
FixStyle(metatomic/kk, FixMetatomicKokkos<LMPDeviceType>);
// clang-format on
#else

#ifndef LMP_FIX_METATOMIC_KOKKOS_H
#define LMP_FIX_METATOMIC_KOKKOS_H

#include "fix_metatomic.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class MetatomicSystemAdaptorKokkos;

template<class DeviceType>
class FixMetatomicKokkos : public FixMetatomic {
public:
   typedef ArrayTypes<DeviceType> AT;

   FixMetatomicKokkos(class LAMMPS *, int, char **);
   ~FixMetatomicKokkos();

   void init() override;
   void initial_integrate(int) override;
   void post_force(int) override;
   void final_integrate() override;

private:
   void pick_device(c10::Device& device, const char* requested) override;

   // Kokkos view for type mapping
   Kokkos::View<int32_t*, Kokkos::LayoutRight, DeviceType> type_mapping_kk;
};

}    // namespace LAMMPS_NS

#endif
#endif
