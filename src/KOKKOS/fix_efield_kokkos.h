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
FixStyle(efield/kk,FixFixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/device,FixFixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/host,FixFixEfieldKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EFIELD_KOKKOS_H
#define LMP_FIX_EFIELD_KOKKOS_H

#include "fix_efield.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixEfieldKokkos : public FixEfield {
  public:
    FixEfieldKokkos(class LAMMPS *, int, char **);

    void post_force(int) override;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int, double &) const;

  private:
    typename ArrayTypes<DeviceType>::t_x_array x;
    typename ArrayTypes<DeviceType>::t_f_array f;
    typename ArrayTypes<DeviceType>::t_int_1d type;
    typename ArrayTypes<DeviceType>::t_int_1d mask;
    typename ArrayTypes<DeviceType>::t_float_1d_randomread q;
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_EFIELD_KOKKOS_H
#endif // FIX_CLASS
