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

/* ----------------------------------------------------------------------
   Contributing authors: Vladislav Galigerov (HSE),
                         Daniil Pavlov (MIPT)
------------------------------------------------------------------------- */

#include "fix_wall_flow_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "math_const.h"
#include "memory_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

template <class DeviceType>
FixWallFlowKokkos<DeviceType>::FixWallFlowKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixWallFlow(lmp, narg, arg), rand_pool(rndseed + comm->me)
{
  kokkosable = 1;
  exchange_comm_device = sort_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | RMASS_MASK | TYPE_MASK | MASK_MASK;
  datamask_modify = V_MASK;

  memory->destroy(current_segment);
  current_segment = nullptr;
  grow_arrays(atomKK->nmax);

  d_walls = d_walls_t("FixWallFlowKokkos::walls", walls.size());
  auto h_walls = Kokkos::create_mirror_view(d_walls);
  for (int i = 0; i < (int) walls.size(); ++i) h_walls(i) = walls[i];
  Kokkos::deep_copy(d_walls, h_walls);
}

template <class DeviceType> FixWallFlowKokkos<DeviceType>::~FixWallFlowKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_current_segment, current_segment);
}

template <class DeviceType> void FixWallFlowKokkos<DeviceType>::init()
{
  atomKK->sync(execution_space, datamask_read);
  k_current_segment.template sync<DeviceType>();
  d_x = atomKK->k_x.template view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixWallFlowInit>(0, atom->nlocal), *this);
  copymode = 0;

  k_current_segment.template modify<DeviceType>();
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void FixWallFlowKokkos<DeviceType>::operator()(TagFixWallFlowInit,
                                                                      const int &i) const
{
  double pos = d_x(i, flowax);
  d_current_segment(i) = compute_current_segment_kk(pos);
}

template <class DeviceType> void FixWallFlowKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(execution_space, datamask_read);
  k_current_segment.template sync<DeviceType>();

  d_x = atomKK->k_x.template view<DeviceType>();
  d_v = atomKK->k_v.template view<DeviceType>();
  d_type = atomKK->k_type.template view<DeviceType>();
  d_mask = atomKK->k_mask.template view<DeviceType>();
  d_mass = atomKK->k_mass.template view<DeviceType>();
  d_rmass = atomKK->k_rmass.template view<DeviceType>();

  copymode = 1;
  if (d_rmass.data()) {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType, TagFixWallFlowEndOfStep<RMassTag>>(0, atom->nlocal), *this);
  } else {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType, TagFixWallFlowEndOfStep<MassTag>>(0, atom->nlocal), *this);
  }
  copymode = 0;
  atomKK->modified(execution_space, datamask_modify);
  k_current_segment.template modify<DeviceType>();
}

template <class DeviceType>
template <class MTag>
KOKKOS_INLINE_FUNCTION void FixWallFlowKokkos<DeviceType>::operator()(TagFixWallFlowEndOfStep<MTag>,
                                                                      const int &atom_i) const
{
  if (d_mask[atom_i] & groupbit) {
    double pos = d_x(atom_i, flowax);
    int prev_segment = d_current_segment(atom_i);
    d_current_segment(atom_i) = compute_current_segment_kk(pos);
    if (prev_segment != d_current_segment(atom_i)) { generate_velocity_kk<MTag>(atom_i); }
  }
}

template <class DeviceType>
template <class MTag>
KOKKOS_INLINE_FUNCTION void FixWallFlowKokkos<DeviceType>::generate_velocity_kk(int atom_i) const
{
  const int newton_iteration_count = 10;
  double mass = get_mass(MTag(), atom_i);
  const double gamma = 1.0 / std::sqrt(2.0 * kT / mass);
  double delta = gamma * flowvel;

  const double edd = std::exp(-delta * delta) / MathConst::MY_PIS + delta * std::erf(delta);
  const double probability_threshold = 0.5 * (1. + delta / edd);

  double direction = 1.0;

  rand_type_t rand_gen = rand_pool.get_state();

  if (/*random->uniform()*/ rand_gen.drand() > probability_threshold) {
    delta = -delta;
    direction = -direction;
  }

  const double xi_0 = rand_gen.drand();    //random->uniform();
  const double F_inf = edd + delta;
  const double xi = xi_0 * F_inf;
  const double x_0 = (std::sqrt(delta * delta + 2) - delta) * 0.5;
  double x = x_0;
  for (int i = 0; i < newton_iteration_count; ++i) {
    x -= (std::exp(x * x) * MathConst::MY_PIS * (xi - delta * std::erfc(x)) - 1.0) / (x + delta) *
        0.5;
  }

  const double nu = x + delta;
  const double v = nu / gamma;

  d_v(atom_i, flowax) = v * direction;
  d_v(atom_i, (flowax + 1) % 3) =
      /*random->gaussian()*/ rand_gen.normal() / (gamma * MathConst::MY_SQRT2);
  d_v(atom_i, (flowax + 2) % 3) =
      /*random->gaussian()*/ rand_gen.normal() / (gamma * MathConst::MY_SQRT2);

  rand_pool.free_state(rand_gen);
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION int
FixWallFlowKokkos<DeviceType>::compute_current_segment_kk(double pos) const
{
  int result = 0;
  for (; result < (int) d_walls.extent(0) - 1; ++result) {
    if (pos >= d_walls[result] && pos < d_walls[result + 1]) { return result; }
  }
  return -1;    // -1 is "out of box" region
}

template <class DeviceType> void FixWallFlowKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_current_segment.template sync<DeviceType>();
  memoryKK->grow_kokkos(k_current_segment, current_segment, nmax, "WallFlowKK::current_segment");
  k_current_segment.template modify<DeviceType>();

  d_current_segment = k_current_segment.template view<DeviceType>();
  h_current_segment = k_current_segment.template view<LMPHostType>();
}

template <class DeviceType> void FixWallFlowKokkos<DeviceType>::copy_arrays(int i, int j, int)
{
  k_current_segment.template sync<LMPHostType>();
  h_current_segment(j) = h_current_segment(i);
  k_current_segment.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template <class DeviceType>
void FixWallFlowKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_current_segment.sync_device();

  Sorter.sort(LMPDeviceType(), k_current_segment.d_view);

  k_current_segment.modify_device();
}

template <class DeviceType> int FixWallFlowKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_current_segment.sync_host();
  buf[0] = static_cast<double>(h_current_segment(i));
  return 1;
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void FixWallFlowKokkos<DeviceType>::operator()(TagFixWallFlowPackExchange,
                                                                      const int &mysend) const
{
  const int send_i = d_sendlist(mysend);
  const int segment = d_current_segment(send_i);
  d_buf(mysend) = static_cast<double>(segment);

  const int copy_i = d_copylist(mysend);
  if (copy_i > -1) { d_current_segment(send_i) = d_current_segment(copy_i); }
}

template <class DeviceType>
int FixWallFlowKokkos<DeviceType>::pack_exchange_kokkos(const int &nsend,
                                                        DAT::tdual_xfloat_2d &k_buf,
                                                        DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_int_1d k_copylist,
                                                        ExecutionSpace /*space*/)
{
  k_current_segment.template sync<DeviceType>();

  k_buf.template sync<DeviceType>();
  k_sendlist.template sync<DeviceType>();
  k_copylist.template sync<DeviceType>();

  d_sendlist = k_sendlist.view<DeviceType>();
  d_copylist = k_copylist.view<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(k_buf.template view<DeviceType>().data(),
                                                          k_buf.extent(0) * k_buf.extent(1));

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixWallFlowPackExchange>(0, nsend),
                       *this);

  copymode = 0;

  k_buf.template modify<DeviceType>();
  k_current_segment.template modify<DeviceType>();

  return nsend;
}

template <class DeviceType> int FixWallFlowKokkos<DeviceType>::unpack_exchange(int i, double *buf)
{
  k_current_segment.sync_host();
  h_current_segment(i) = static_cast<int>(buf[0]);
  k_current_segment.modify_host();
  return 1;
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void FixWallFlowKokkos<DeviceType>::operator()(TagFixWallFlowUnpackExchange,
                                                                      const int &i) const
{
  int index = d_indices(i);
  if (index > -1) { d_current_segment(index) = static_cast<int>(d_buf(i)); }
}

template <class DeviceType>
void FixWallFlowKokkos<DeviceType>::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                                                           DAT::tdual_int_1d &k_indices, int nrecv,
                                                           int /*nrecv1*/, int /*nextrarecv1*/,
                                                           ExecutionSpace /*space*/)
{
  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(k_buf.template view<DeviceType>().data(),
                                                          k_buf.extent(0) * k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixWallFlowUnpackExchange>(0, nrecv),
                       *this);
  copymode = 0;

  k_current_segment.template modify<DeviceType>();
}

namespace LAMMPS_NS {
template class FixWallFlowKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallFlowKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
