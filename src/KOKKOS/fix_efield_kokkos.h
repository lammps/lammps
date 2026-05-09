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

FixStyle(efield/kk,        FixEfieldKokkos<LMPDeviceType,FixEfield>);
FixStyle(efield/kk/device, FixEfieldKokkos<LMPDeviceType,FixEfield>);
FixStyle(efield/kk/host,   FixEfieldKokkos<LMPHostType,FixEfield>);

  FixStyle(efield/tip4p/kk,        FixEfieldKokkos<LMPDeviceType,FixEfieldTIP4P>);
  FixStyle(efield/tip4p/kk/device, FixEfieldKokkos<LMPDeviceType,FixEfieldTIP4P>);
  FixStyle(efield/tip4p/kk/host,   FixEfieldKokkos<LMPHostType,FixEfieldTIP4P>);

// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EFIELD_KOKKOS_H
#define LMP_FIX_EFIELD_KOKKOS_H

#include "fix_efield.h"
#include "fix_efield_tip4p.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

// Curiously Recurring Template Pattern (CRTP):
// FixEfieldBase is either FixEfield (plain) or FixEfieldTIP4P.
// TIP4P = true when the base is FixEfieldTIP4P or derived from it.

template<class DeviceType, class FixEfieldBase>
class FixEfieldKokkos : public FixEfieldBase
{
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixEfieldKokkos(class LAMMPS *, int, char **);
  ~FixEfieldKokkos() override;
  void init() override;
  void post_force(int) override;

  static constexpr bool TIP4P = std::is_base_of_v<FixEfieldTIP4P, FixEfieldBase>;

 protected:

  using Pointers::atom;
  using Pointers::atomKK;
  using Pointers::memory;
  using Pointers::memoryKK;
  using Pointers::update;
  using Pointers::error;
  using Pointers::domain;

  using Fix::style;
  using Fix::kokkosable;
  using Fix::copymode;
  using Fix::execution_space;
  using Fix::datamask_read;
  using Fix::datamask_modify;
  using Fix::evflag;
  using Fix::vflag_global;
  using Fix::vflag_atom;
  using Fix::virial;
  using Fix::vatom;
  using Fix::maxvatom;
  using Fix::v_init;
  using Fix::groupbit;
  using Fix::igroup;

  using FixEfieldBase::force_flag;
  using FixEfieldBase::maxatom;
  using FixEfieldBase::region;
  using FixEfieldBase::varflag;
  using FixEfieldBase::qflag;
  using FixEfieldBase::muflag;
  using FixEfieldBase::ex;
  using FixEfieldBase::ey;
  using FixEfieldBase::ez;
  using FixEfieldBase::qe2f;
  using FixEfieldBase::xstyle;
  using FixEfieldBase::ystyle;
  using FixEfieldBase::zstyle;
  using FixEfieldBase::pstyle;
  using FixEfieldBase::estyle;
  using FixEfieldBase::xvar;
  using FixEfieldBase::yvar;
  using FixEfieldBase::zvar;
  using FixEfieldBase::evar;
  using FixEfieldBase::update_efield_variables;
  using FixEfieldBase::fsum;
  using FixEfieldBase::efield;

 private:

  DAT::ttransform_kkfloat_2d k_efield;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;
};

} // namespace LAMMPS_NS

#endif // !LMP_FIX_EFIELD_KOKKOS_H
#endif
