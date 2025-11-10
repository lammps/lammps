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
FixStyle(spring/self/kk,FixSpringSelfKokkos<LMPDeviceType>);
FixStyle(spring/self/kk/device,FixSpringSelfKokkos<LMPDeviceType>);
FixStyle(spring/self/kk/host,FixSpringSelfKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SPRING_SELF_KOKKOS_H
#define LMP_FIX_SPRING_SELF_KOKKOS_H

#include "fix_spring_self.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

struct TagFixSpringSelfUnpackExchange{};

template<class DeviceType>
class FixSpringSelfKokkos : public FixSpringSelf, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type;

  FixSpringSelfKokkos(class LAMMPS *, int, char **);
  ~FixSpringSelfKokkos() override;
  void init() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void pack_exchange_item(const int&, int &, const bool &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringSelfUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nrecv1,int nrecv1extra,
                              ExecutionSpace space) override;


  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  int nrecv1,nextrarecv1;

  DAT::ttransform_kkfloat_1d_3_lr k_xoriginal;
  typename AT::t_kkfloat_1d_3_lr d_xoriginal;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;

  int nsend;

  typename AT::t_int_2d d_sendlist;
  typename AT::t_double_1d_um d_buf;

  typename AT::t_int_1d d_exchange_sendlist;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  typename AT::t_int_scalar d_count;
  HAT::t_int_scalar h_count;

  double **xoriginal_tmp;    // original coords of atoms

};

template <class DeviceType>
struct FixSpringSelfKokkosPackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;
  FixSpringSelfKokkos<DeviceType> c;
  FixSpringSelfKokkosPackExchangeFunctor(FixSpringSelfKokkos<DeviceType>* c_ptr):c(*c_ptr) {};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, int &offset, const bool &final) const {
    c.pack_exchange_item(i, offset, final);
  }
};

}

#endif
#endif

