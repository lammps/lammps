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

#include "compute_temp_profile_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempProfileKokkos<DeviceType>::ComputeTempProfileKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempProfile(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbin = 0;

  // per-bin accumulators (nbins, ncount are set by the base constructor).
  // the sums are accumulated in double regardless of the precision of the
  // per-atom data, and the host mirrors are persistent so no per-step
  // allocations happen on GPU builds.

  d_vbin = t_double_2d("temp/profile/kk:vbin", nbins, ncount);
  d_binave = t_double_2d("temp/profile/kk:binave", nbins, ncount);
  h_vbin = Kokkos::create_mirror_view(d_vbin);
  h_binave = Kokkos::create_mirror_view(d_binave);
}

/* ----------------------------------------------------------------------
   compute average COM velocity in each bin (device), mirrors
   ComputeTempProfile::bin_average but with the per-atom work on device
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::bin_average_kk()
{
  if (box_change) bin_setup();

  // copy the binning frame into device-friendly scalars.  for triclinic
  // boxes the base class init() pointed boxlo/boxhi/prd at the lamda-space
  // values, so the same bin formula applies to lamda coordinates.

  for (int d = 0; d < 3; d++) {
    m_boxlo[d] = static_cast<KK_FLOAT>(boxlo[d]);
    m_boxhi[d] = static_cast<KK_FLOAT>(boxhi[d]);
    m_prd[d] = static_cast<KK_FLOAT>(prd[d]);
    m_invdelta[d] = static_cast<KK_FLOAT>(invdelta[d]);
    m_periodicity[d] = periodicity[d];
  }

  int nlocal = atom->nlocal;

  if (atom->nmax > maxbin) {
    maxbin = atom->nmax;
    d_bin = typename AT::t_int_1d("temp/profile/kk:bin", maxbin);
  }

  // assign each atom to a bin on the device.  for triclinic boxes convert
  // the coordinates to lamda space around the kernel; DomainKokkos does the
  // conversions with device kernels that sync X themselves.

  if (triclinic) domain->x2lamda(nlocal);

  atomKK->sync(execution_space, X_MASK | MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileBin>(0,nlocal),*this);
  copymode = 0;

  if (triclinic) domain->lamda2x(nlocal);

  // keep the host bin[] array current: the inherited per-atom
  // remove_bias()/restore_bias() (used by non-KOKKOS thermostats pointed at
  // this compute) and compute_array() index into it.  this copies nlocal
  // ints to the host, which is cheap next to the per-atom velocity data.

  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(bin);
    memory->create(bin,maxatom,"temp/profile:bin");
  }
  if (nlocal > 0) {
    Kokkos::View<int*,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>> h_bin_um(bin,nlocal);
    Kokkos::deep_copy(h_bin_um, Kokkos::subview(d_bin, Kokkos::pair<int,int>(0,nlocal)));
  }

  // sum each atom's mass-weighted velocity, mass, and count into its bin

  atomKK->sync(execution_space, V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  atomKK->k_mass.sync<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  if (atomKK->rmass) rmass = atomKK->k_rmass.view<DeviceType>();
  else mass = atomKK->k_mass.view<DeviceType>();

  Kokkos::deep_copy(d_vbin, 0.0);

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileScatter<1> >(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileScatter<0> >(0,nlocal),*this);
  copymode = 0;

  // sum bins across procs on the host (bit-faithful to the CPU reduction)

  Kokkos::deep_copy(h_vbin, d_vbin);
  for (int i = 0; i < nbins; i++)
    for (int j = 0; j < ncount; j++) vbin[i][j] = h_vbin(i,j);

  MPI_Allreduce(vbin[0],binave[0],nbins*ncount,MPI_DOUBLE,MPI_SUM,world);

  // compute ave COM velocity in each bin, checking for no particles

  int nc2 = ncount-2;
  int nc1 = ncount-1;
  for (int i = 0; i < nbins; i++)
    if (binave[i][nc1] > 0.0)
      for (int j = 0; j < nc2; j++) binave[i][j] /= binave[i][nc2];

  for (int i = 0; i < nbins; i++)
    for (int j = 0; j < ncount; j++) h_binave(i,j) = binave[i][j];
  Kokkos::deep_copy(d_binave, h_binave);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileBin, const int &i) const
{
  if (!(mask[i] & groupbit)) { d_bin[i] = 0; return; }

  int ibinx = 0, ibiny = 0, ibinz = 0;
  KK_FLOAT coord;

  if (nbinx > 1) {
    coord = x(i,0);
    if (m_periodicity[0]) {
      if (coord < m_boxlo[0]) coord += m_prd[0];
      if (coord >= m_boxhi[0]) coord -= m_prd[0];
    }
    ibinx = static_cast<int>((coord - m_boxlo[0]) * m_invdelta[0]);
    if (ibinx < 0) ibinx = 0;
    if (ibinx > nbinx-1) ibinx = nbinx-1;
  }
  if (nbiny > 1) {
    coord = x(i,1);
    if (m_periodicity[1]) {
      if (coord < m_boxlo[1]) coord += m_prd[1];
      if (coord >= m_boxhi[1]) coord -= m_prd[1];
    }
    ibiny = static_cast<int>((coord - m_boxlo[1]) * m_invdelta[1]);
    if (ibiny < 0) ibiny = 0;
    if (ibiny > nbiny-1) ibiny = nbiny-1;
  }
  if (nbinz > 1) {
    coord = x(i,2);
    if (m_periodicity[2]) {
      if (coord < m_boxlo[2]) coord += m_prd[2];
      if (coord >= m_boxhi[2]) coord -= m_prd[2];
    }
    ibinz = static_cast<int>((coord - m_boxlo[2]) * m_invdelta[2]);
    if (ibinz < 0) ibinz = 0;
    if (ibinz > nbinz-1) ibinz = nbinz-1;
  }

  d_bin[i] = nbinx*nbiny*ibinz + nbinx*ibiny + ibinx;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileScatter<RMASS>, const int &i) const
{
  if (mask[i] & groupbit) {
    const int ibin = d_bin[i];
    double massone;
    if (RMASS) massone = static_cast<double>(rmass[i]);
    else massone = static_cast<double>(mass[type[i]]);
    if (xflag) Kokkos::atomic_add(&d_vbin(ibin,ivx), massone*(double)v(i,0));
    if (yflag) Kokkos::atomic_add(&d_vbin(ibin,ivy), massone*(double)v(i,1));
    if (zflag) Kokkos::atomic_add(&d_vbin(ibin,ivz), massone*(double)v(i,2));
    Kokkos::atomic_add(&d_vbin(ibin,ncount-2), massone);
    Kokkos::atomic_add(&d_vbin(ibin,ncount-1), 1.0);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempProfileKokkos<DeviceType>::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  bin_average_kk();

  int nlocal = atom->nlocal;

  CTEMP t_kk;
  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  double t = t_kk.t0;
  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileScalar<RMASS>, const int &i, CTEMP& t_kk) const
{
  if (mask[i] & groupbit) {
    const int ibin = d_bin[i];
    double vt0 = static_cast<double>(v(i,0)), vt1 = static_cast<double>(v(i,1)), vt2 = static_cast<double>(v(i,2));
    if (xflag) vt0 -= d_binave(ibin,ivx);
    if (yflag) vt1 -= d_binave(ibin,ivy);
    if (zflag) vt2 -= d_binave(ibin,ivz);
    double massone;
    if (RMASS) massone = static_cast<double>(rmass[i]);
    else massone = static_cast<double>(mass[type[i]]);
    t_kk.t0 += (vt0*vt0 + vt1*vt1 + vt2*vt2) * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  bin_average_kk();

  int nlocal = atom->nlocal;

  CTEMP t_kk;
  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileVector<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  double t[6];
  t[0] = t_kk.t0; t[1] = t_kk.t1; t[2] = t_kk.t2;
  t[3] = t_kk.t3; t[4] = t_kk.t4; t[5] = t_kk.t5;

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileVector<RMASS>, const int &i, CTEMP& t_kk) const
{
  if (mask[i] & groupbit) {
    const int ibin = d_bin[i];
    double vt0 = static_cast<double>(v(i,0)), vt1 = static_cast<double>(v(i,1)), vt2 = static_cast<double>(v(i,2));
    if (xflag) vt0 -= d_binave(ibin,ivx);
    if (yflag) vt1 -= d_binave(ibin,ivy);
    if (zflag) vt2 -= d_binave(ibin,ivz);
    double massone;
    if (RMASS) massone = static_cast<double>(rmass[i]);
    else massone = static_cast<double>(mass[type[i]]);
    t_kk.t0 += massone * vt0*vt0;
    t_kk.t1 += massone * vt1*vt1;
    t_kk.t2 += massone * vt2*vt2;
    t_kk.t3 += massone * vt0*vt1;
    t_kk.t4 += massone * vt0*vt2;
    t_kk.t5 += massone * vt1*vt2;
  }
}

/* ----------------------------------------------------------------------
   per-bin temperature output (rare diagnostic path): mirrors the base
   compute_array(), but bins via bin_average_kk() -- the base bin_assign()
   cannot be used here because DomainKokkos::x2lamda(int) converts the
   coordinates on the device, which the host loop would not see
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::compute_array()
{
  int i,ibin;
  double vthermal[3];

  invoked_array = update->ntimestep;

  bin_average_kk();

  atomKK->sync(Host, V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  atomKK->k_mass.sync<LMPHostType>();

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nbins; i++) tbin[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) vthermal[0] = v[i][0] - binave[ibin][ivx];
      else vthermal[0] = v[i][0];
      if (yflag) vthermal[1] = v[i][1] - binave[ibin][ivy];
      else vthermal[1] = v[i][1];
      if (zflag) vthermal[2] = v[i][2] - binave[ibin][ivz];
      else vthermal[2] = v[i][2];

      if (rmass)
        tbin[ibin] += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
                       vthermal[2]*vthermal[2]) * rmass[i];
      else
        tbin[ibin] += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
                       vthermal[2]*vthermal[2]) * mass[type[i]];
    }

  MPI_Allreduce(tbin,tbinall,nbins,MPI_DOUBLE,MPI_SUM,world);

  double totcount = 0.0;
  for (i = 0; i < nbins; i++) {
    array[i][0] = binave[i][ncount-1];
    totcount += array[i][0];
  }
  double nper = (totcount > 0.0) ? domain->dimension - (extra_dof + fix_dof)/totcount : 0.0;
  double dofbin, tfactorbin;
  for (i = 0; i < nbins; i++) {
    if (array[i][0] > 0.0) {
      dofbin = nper*array[i][0] - nstreaming;
      tfactorbin = (dofbin > 0.0) ? force->mvv2e / (dofbin * force->boltz) : 0.0;
      array[i][1] = tfactorbin*tbinall[i];
    } else array[i][1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   bias is the binned streaming (COM) velocity computed by the last
   compute_scalar()/compute_vector() (mirrors the CPU assumption)
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::remove_bias_all()
{
  remove_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::remove_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileRemoveBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileRemoveBias, const int &i) const
{
  if (mask[i] & groupbit) {
    const int ibin = d_bin[i];
    if (xflag) v(i,0) -= static_cast<KK_FLOAT>(d_binave(ibin,ivx));
    if (yflag) v(i,1) -= static_cast<KK_FLOAT>(d_binave(ibin,ivy));
    if (zflag) v(i,2) -= static_cast<KK_FLOAT>(d_binave(ibin,ivz));
  }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::restore_bias_all()
{
  restore_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

template<class DeviceType>
void ComputeTempProfileKokkos<DeviceType>::restore_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempProfileRestoreBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempProfileKokkos<DeviceType>::operator()(TagComputeTempProfileRestoreBias, const int &i) const
{
  if (mask[i] & groupbit) {
    const int ibin = d_bin[i];
    if (xflag) v(i,0) += static_cast<KK_FLOAT>(d_binave(ibin,ivx));
    if (yflag) v(i,1) += static_cast<KK_FLOAT>(d_binave(ibin,ivy));
    if (zflag) v(i,2) += static_cast<KK_FLOAT>(d_binave(ibin,ivz));
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempProfileKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempProfileKokkos<LMPHostType>;
#endif
}
