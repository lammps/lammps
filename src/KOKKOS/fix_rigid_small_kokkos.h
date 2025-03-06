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
FixStyle(rigid/small/kk,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/device,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/host,FixRigidSmallKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_SMALL_KOKKOS_H

#include "fix_rigid_small.h"
#include "kokkos_base.h"
#include "comm_kokkos.h"
#include <map>

struct TagInitialIntegrate{};
struct TagPackForwardInitial{};
struct TagUnpackForwardInitial{};
template<int SETXFLAG>
struct TagSetXV{};
struct TagUpdateXGC{};

namespace LAMMPS_NS {


template<class DeviceType>
class FixRigidSmallKokkos : public FixRigidSmall, public KokkosBase {

 public:
  typedef EV_FLOAT value_type;
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;


  FixRigidSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidSmallKokkos();
  //int setmask() override; use super
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void write_restart_file(const char *) override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void grow_body() override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int, int,
                              ExecutionSpace space) override;
  int pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                               DAT::tdual_xfloat_1d &k_buf,
                               int pbc_flag, int* pbc) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;

  int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d,
                                          DAT::tdual_xfloat_1d &) override;
  // reverse comm handled by host,
  // only happens when body and bodyown
  // are already on host

  void setup_pre_neighbor() override;
  bigint dof(int) override;
  void deform(int) override;
  void enforce2d() override;
  void zero_momentum() override;
  void zero_rotation() override;
  void *extract(const char *, int &) override;
  double extract_ke();
  double extract_erotational();
  double compute_scalar() override;
  void reset_atom2body() override;
  void image_shift() override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagInitialIntegrate, const int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPackForwardInitial, const int) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUnpackForwardInitial, const int) const;

  template<int SETXFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagSetXV<SETXFLAG>, const int, EV_FLOAT &ev) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateXGC, const int) const;

  void compute_forces_and_torques_kokkos();

 protected:

  void set_xv_kokkos(int);
  void apply_langevin_thermostat_kokkos();

  // TODO: Use AT stuff
  using ImageIntView1D = Kokkos::View<imageint*, Kokkos::LayoutRight, DeviceType>;
  using TagIntView1D = Kokkos::View<tagint*, Kokkos::LayoutRight, DeviceType>;
  using IntView1D = Kokkos::View<int*, Kokkos::LayoutRight, DeviceType>;
  using IntView2D = Kokkos::View<int**, Kokkos::LayoutRight, DeviceType>;
  using View1D = Kokkos::View<F_FLOAT*, Kokkos::LayoutRight, DeviceType>;
  using View2D = Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType>;

  using Range1D = Kokkos::RangePolicy<DeviceType>;

  void copy_body_host();
  void copy_body_device();
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT&, int, double[6], double[3], double[3], double[3]) const;
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT&, int, double[6]) const;

  IntView1D d_bodyown;
  TagIntView1D d_bodytag;
  IntView1D d_atom2body;
  ImageIntView1D d_xcmimage;
  View2D d_displace, d_vatom, d_langextra;

  int max_body_sent=0;
  std::map<int,int> n_body_recv, first_body;
  std::map<int*,int> n_body_sent;
  std::map<int*,IntView1D> d_body_sendlists;


  IntView1D d_sendlist;
  View1D d_buf;
  int first;

  CommKokkos *commKK;

  Kokkos::View<Body*, DeviceType> d_body;

  double xbox, ybox, zbox, xprd, yprd, zprd, xy, xz, yz;
  typename AT::t_x_array d_x;
  typename AT::t_v_array d_v;
  typename AT::t_f_array d_f;

  View1D d_rmass, d_mass;
  IntView1D d_type;
};

KOKKOS_INLINE_FUNCTION
void copy_body(Body *dest, Body *src){
  memcpy(dest, src, sizeof(Body));
}

}    // namespace LAMMPS_NS


#endif
#endif
