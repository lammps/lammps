// clang-format off
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
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "fix_spring_self_kokkos.h"

#include "atom_kokkos.h"
#include "update.h"
#include "domain_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringSelfKokkos<DeviceType>::FixSpringSelfKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSpringSelf(lmp, narg, arg)
{
  kokkosable = 1;
  exchange_comm_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  xoriginal_tmp = xoriginal;
  xoriginal = nullptr;

  int nmax = atom->nmax;
  grow_arrays(nmax);

  for (int i = 0; i < atom->nlocal; i++) {
    k_xoriginal.h_view(i,0) = xoriginal_tmp[i][0];
    k_xoriginal.h_view(i,1) = xoriginal_tmp[i][1];
    k_xoriginal.h_view(i,2) = xoriginal_tmp[i][2];
  }

  k_xoriginal.modify_host();

  d_count = typename AT::t_int_scalar("spring/self:count");
  h_count = Kokkos::create_mirror_view(d_count);

  memory->destroy(xoriginal_tmp);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringSelfKokkos<DeviceType>::~FixSpringSelfKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_xoriginal,xoriginal);
  xoriginal = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::init()
{
  FixSpringSelf::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double espring_kk;

  k_xoriginal.modify<LMPHostType>();
  k_xoriginal.sync<DeviceType>();

  copymode = 1;

  {
  // local variables for lambda capture
  auto prd = Few<double,3>(domain->prd);
  auto h = Few<double,6>(domain->h);
  auto triclinic = domain->triclinic;
  auto l_k = k;
  auto l_xoriginal = d_xoriginal;

  auto l_x = x;
  auto l_f = f;
  auto l_mask = mask;
  auto l_image = image;
  auto l_groupbit = groupbit;
  auto l_xflag = xflag;
  auto l_yflag = yflag;
  auto l_zflag = zflag;

  Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(const int& i, double& espring_kk) {
    if (l_mask[i] & l_groupbit) {
      Few<double,3> x_i;
      x_i[0] = l_x(i,0);
      x_i[1] = l_x(i,1);
      x_i[2] = l_x(i,2);
      auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,l_image(i));
      auto dx = unwrap[0] - l_xoriginal(i, 0);
      auto dy = unwrap[1] - l_xoriginal(i, 1);
      auto dz = unwrap[2] - l_xoriginal(i, 2);
      if (!l_xflag) dx = 0.0;
      if (!l_yflag) dy = 0.0;
      if (!l_zflag) dz = 0.0;
      l_f(i,0) -= l_k*dx;
      l_f(i,1) -= l_k*dy;
      l_f(i,2) -= l_k*dz;
      espring_kk += l_k * (dx*dx + dy*dy + dz*dz);
    }
  },espring_kk);
  }

  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  espring = 0.5*espring_kk;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::grow_arrays(int nmax)
{
  memoryKK->grow_kokkos(k_xoriginal,xoriginal,nmax,"spring/self:xoriginal");
  d_xoriginal = k_xoriginal.view<DeviceType>();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_xoriginal.sync_host();

  FixSpringSelf::copy_arrays(i,j,delflag);

  k_xoriginal.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixSpringSelfKokkos<DeviceType>::pack_exchange_item(const int &mysend, int &offset, const bool &/*final*/) const
{
  const int i = d_exchange_sendlist(mysend);

  int m = nsend + offset;
  d_buf[mysend] = m;
  d_buf[m++] = d_xoriginal(i,0);
  d_buf[m++] = d_xoriginal(i,1);
  d_buf[m++] = d_xoriginal(i,2);
  if (mysend == nsend-1) d_count() = m;
  offset = m - nsend;

  const int j = d_copylist(mysend);
  if (j > -1) {
    d_xoriginal(i,0) = d_xoriginal(j,0);
    d_xoriginal(i,1) = d_xoriginal(j,1);
    d_xoriginal(i,2) = d_xoriginal(j,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixSpringSelfKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{

  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_exchange_sendlist.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_copylist = k_copylist.view<DeviceType>();
  d_exchange_sendlist = k_exchange_sendlist.view<DeviceType>();
  this->nsend = nsend;


  k_xoriginal.template sync<DeviceType>();

  Kokkos::deep_copy(d_count,0);

  copymode = 1;

  FixSpringSelfKokkosPackExchangeFunctor<DeviceType> pack_exchange_functor(this);
  Kokkos::parallel_scan(nsend,pack_exchange_functor);

  copymode = 0;

  k_buf.modify<DeviceType>();

  if (space == Host) k_buf.sync<LMPHostType>();
  else k_buf.sync<LMPDeviceType>();

  k_xoriginal.template modify<DeviceType>();

  Kokkos::deep_copy(h_count,d_count);

  return h_count();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixSpringSelfKokkos<DeviceType>::operator()(TagFixSpringSelfUnpackExchange, const int &i) const
{
  int index = d_indices(i);

  if (index > -1) {
    int m = d_buf[i];
    if (i >= nrecv1)
      m = nextrarecv1 + d_buf[nextrarecv1 + i - nrecv1];

    d_xoriginal(index,0) = static_cast<tagint> (d_buf[m++]);
    d_xoriginal(index,1) = static_cast<tagint> (d_buf[m++]);
    d_xoriginal(index,2) = static_cast<tagint> (d_buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixSpringSelfKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int nrecv1, int nextrarecv1,
  ExecutionSpace /*space*/)
{
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  this->nrecv1 = nrecv1;
  this->nextrarecv1 = nextrarecv1;

  k_xoriginal.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringSelfUnpackExchange>(0,nrecv),*this);

  copymode = 0;

  k_xoriginal.template modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixSpringSelfKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_xoriginal.sync_host();

  int m = FixSpringSelf::pack_exchange(i,buf);

  k_xoriginal.modify_host();

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixSpringSelfKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_xoriginal.sync_host();

  int m = FixSpringSelf::unpack_exchange(nlocal,buf);

  k_xoriginal.modify_host();

  return m;
}

namespace LAMMPS_NS {
template class FixSpringSelfKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringSelfKokkos<LMPHostType>;
#endif
}

