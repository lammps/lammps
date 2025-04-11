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
FixStyle(cmap/kk,FixCMAPKokkos<LMPDeviceType>);
FixStyle(cmap/kk/device,FixCMAPKokkos<LMPDeviceType>);
FixStyle(cmap/kk/host,FixCMAPKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_CMAP_KOKKOS_H
#define LMP_FIX_CMAP_KOKKOS_H

#include "fix_cmap.h"

#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixCmapPreNeighbor{};
struct TagFixCmapPostForce{};

template<class DeviceType>
class FixCMAPKokkos : public FixCMAP, public KokkosBase {
  typedef ArrayTypes<DeviceType> AT;

  public:
    FixCMAPKokkos(class LAMMPS *, int, char **);
    ~FixCMAPKokkos() override;

    void init() override;
    void pre_neighbor() override;
    void post_force(int) override;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixCmapPreNeighbor, const int, int&, const bool) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixCmapPostForce, const int, double&) const;

    void grow_arrays(int) override;
    void copy_arrays(int, int, int) override;
    void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
    void set_arrays(int) override;
    int pack_exchange(int, double *) override;
    int unpack_exchange(int, double *) override;

    int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

    void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nrecv1,int nrecv1extra,
                              ExecutionSpace space) override;

  protected:

    int nlocal;

    typename AT::t_x_array d_x;
    typename AT::t_f_array d_f;

    DAT::tdual_int_1d k_sametag;
    typename AT::t_int_1d d_sametag;
    int map_style;
    DAT::tdual_int_1d k_map_array;
    dual_hash_type k_map_hash;

    typename AT::t_int_scalar d_count;
    HAT::t_int_scalar h_count;

    DAT::tdual_int_1d k_num_crossterm;
    typename AT::t_int_1d d_num_crossterm;

    DAT::tdual_int_2d k_crossterm_type;
    typename AT::t_int_2d d_crosstermlist, d_crossterm_type;

    DAT::tdual_tagint_2d k_crossterm_atom1, k_crossterm_atom2, k_crossterm_atom3;
    DAT::tdual_tagint_2d k_crossterm_atom4, k_crossterm_atom5;
    typename AT::t_tagint_2d d_crossterm_atom1, d_crossterm_atom2, d_crossterm_atom3;
    typename AT::t_tagint_2d d_crossterm_atom4, d_crossterm_atom5;

    DAT::tdual_float_1d k_g_axis;
    typename AT::t_float_1d d_g_axis;

    DAT::tdual_float_3d k_cmapgrid, k_d1cmapgrid, k_d2cmapgrid, k_d12cmapgrid;
    typename AT::t_float_3d d_cmapgrid, d_d1cmapgrid, d_d2cmapgrid, d_d12cmapgrid;

    // calculate dihedral angles
    KOKKOS_INLINE_FUNCTION
    double dihedral_angle_atan2(double, double, double, double, double, double, double, double,
      double, double) const;

    // perform bicubic interpolation at point of interest
    KOKKOS_INLINE_FUNCTION
    void bc_interpol(double, double, int, int, double *, double *, double *, double *,
      double &, double &, double &) const;

    // copied from Domain
    KOKKOS_INLINE_FUNCTION
    int closest_image(const int, int) const;

};

} // namespace LAMMPS_NS

#endif // LMP_FIX_CMAP_KOKKOS_H
#endif // FIX_CLASS
