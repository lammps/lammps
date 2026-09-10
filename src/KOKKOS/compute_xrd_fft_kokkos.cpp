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
   KOKKOS version of compute xrd/fft.

   The whole non-uniform transform runs on the device: the extent of a rank's
   atoms in mesh coordinates, the grouping of those atoms by element, the
   spreading, the packing and unpacking of the exchange, the transform itself,
   and the accumulation of the structure factor over the modes a rank owns.
   Only the decisions about which pieces of the mesh have to move, which are
   integer arithmetic on a handful of numbers, and the final assembly of the
   output array stay on the host.

   Contributing author: derived from compute xrd/fft
------------------------------------------------------------------------- */

#include "compute_xrd_fft_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "fft3d_kokkos.h"
#include "group.h"
#include "kokkos.h"
#include "memory.h"
#include "platform.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeXRDFFTKokkos<DeviceType>::ComputeXRDFFTKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeXRDFFT(lmp, narg, arg), fftkk(nullptr)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | TYPE_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;

  device_ready = 0;
  gpu_aware = lmp->kokkos->gpu_aware_flag;

  // ExecutionSpaceFromDevice yields HostKK rather than the legacy Host, so the
  // test is against Device.  without it the kk/host variant of a GPU build
  // would stage its messages through a host mirror of the buffers it is
  // already computing in, and the copy back would overwrite them

  mpi_direct = gpu_aware || (execution_space != Device);
  bucket_maxatoms = 0;
  spread_lo = spread_hi = 0;
  nfoot_kk = 0;
  order_kk = order;
  nlower_kk = nlower;
  bufoff = 0;
  slot_off = 0;

  if (order > MAXORDER)
    error->all(FLERR,"Compute xrd/fft/kk: order must be at most {}; compute xrd/fft without "
               "the kk suffix has no such limit",MAXORDER);

  k_overrun = DAT::tdual_int_scalar("xrd/fft/kk:overrun");
  d_overrun = k_overrun.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeXRDFFTKokkos<DeviceType>::~ComputeXRDFFTKokkos()
{
  if (copymode) return;

  deallocate();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::deallocate()
{
  delete fftkk;
  fftkk = nullptr;

  k_density_own = FFT_DAT::tdual_FFT_SCALAR_1d();
  k_density_slab = FFT_DAT::tdual_FFT_SCALAR_1d();
  k_work1 = FFT_DAT::tdual_FFT_SCALAR_1d();
  k_sendbuf = FFT_DAT::tdual_FFT_SCALAR_1d();
  k_recvbuf = FFT_DAT::tdual_FFT_SCALAR_1d();
  d_density_own = typename FFT_AT::t_FFT_SCALAR_1d();
  d_density_slab = typename FFT_AT::t_FFT_SCALAR_1d();
  d_work1 = typename FFT_AT::t_FFT_SCALAR_1d();
  d_sendbuf = typename FFT_AT::t_FFT_SCALAR_1d();
  d_recvbuf = typename FFT_AT::t_FFT_SCALAR_1d();
  d_unpacksrc = typename FFT_AT::t_FFT_SCALAR_1d();

  k_slot_atoms = DAT::tdual_int_1d();
  d_slot_atoms = typename AT::t_int_1d();

  bucket_maxatoms = 0;
  device_ready = 0;

  ComputeXRDFFT::deallocate();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::init()
{
  const int was_setup = setup_done;

  ComputeXRDFFT::init();

  if (!was_setup) {
    setup_device();
    copy_scaling();
  }
}

/* ----------------------------------------------------------------------
   the mesh and the transform, in place of the host buffers and FFT3d
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::allocate_mesh()
{
  k_density_slab = FFT_DAT::tdual_FFT_SCALAR_1d("xrd/fft/kk:density_slab",MAX(nfft,1));
  d_density_slab = k_density_slab.template view<DeviceType>();

  k_work1 = FFT_DAT::tdual_FFT_SCALAR_1d("xrd/fft/kk:work1",MAX(2*nfft,1));
  d_work1 = k_work1.template view<DeviceType>();

  if (nfft > 0) {
    int tmp;
    fftkk = new FFT3dKokkos<DeviceType>(lmp,fft_comm,nmesh[0],nmesh[1],nmesh[2],
                                        fftlo[0],ffthi[0],fftlo[1],ffthi[1],fftlo[2],ffthi[2],
                                        fftlo[0],ffthi[0],fftlo[1],ffthi[1],fftlo[2],ffthi[2],
                                        0,0,&tmp,0,0,gpu_aware);
  }
}

/* ----------------------------------------------------------------------
   device copies of the tables the mesh geometry fixes
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::setup_device()
{
  const int ncheb = 3*order*KB_NCHEB;
  k_kb_cheb = DAT::tdual_double_1d("xrd/fft/kk:kb_cheb",ncheb);
  auto h_kb_cheb = k_kb_cheb.view_host();
  for (int i = 0; i < ncheb; i++) h_kb_cheb(i) = kb_cheb[i];
  k_kb_cheb.modify_host();
  k_kb_cheb.template sync<DeviceType>();
  d_kb_cheb = k_kb_cheb.template view<DeviceType>();

  k_slot_of_type = DAT::tdual_int_1d("xrd/fft/kk:slot_of_type",MAX(ntypes,1));
  auto h_slot_of_type = k_slot_of_type.view_host();
  for (int i = 0; i < ntypes; i++) h_slot_of_type(i) = slot_of_type[i];
  k_slot_of_type.modify_host();
  k_slot_of_type.template sync<DeviceType>();
  d_slot_of_type = k_slot_of_type.template view<DeviceType>();

  const int n = MAX(nown,1);

  k_own_idx = DAT::tdual_int_1d("xrd/fft/kk:own_idx",n);
  k_own_deconv = DAT::tdual_double_1d("xrd/fft/kk:own_deconv",n);
  auto h_own_idx = k_own_idx.view_host();
  auto h_own_deconv = k_own_deconv.view_host();
  for (int a = 0; a < nown; a++) {
    h_own_idx(a) = own_idx[a];
    h_own_deconv(a) = own_deconv[a];
  }
  k_own_idx.modify_host();
  k_own_idx.template sync<DeviceType>();
  k_own_deconv.modify_host();
  k_own_deconv.template sync<DeviceType>();
  d_own_idx = k_own_idx.template view<DeviceType>();
  d_own_deconv = k_own_deconv.template view<DeviceType>();

  k_own_asf = DAT::tdual_double_1d("xrd/fft/kk:own_asf",MAX(nslot*nown,1));
  d_own_asf = k_own_asf.template view<DeviceType>();

  k_Fre = DAT::tdual_double_1d("xrd/fft/kk:Fre",n);
  k_Fim = DAT::tdual_double_1d("xrd/fft/kk:Fim",n);
  d_Fre = k_Fre.template view<DeviceType>();
  d_Fim = k_Fim.template view<DeviceType>();

  device_ready = 1;
}

/* ----------------------------------------------------------------------
   the scattering factors change with the box, the tables above do not
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::copy_scaling()
{
  if (!device_ready) return;

  const int n = nslot*nown;
  auto h_own_asf = k_own_asf.view_host();
  for (int i = 0; i < n; i++) h_own_asf(i) = own_asf[i];
  k_own_asf.modify_host();
  k_own_asf.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::refresh_scaling()
{
  ComputeXRDFFT::refresh_scaling();
  copy_scaling();
}

/* ----------------------------------------------------------------------
   everything a kernel reads that is not a view
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::set_kernel_state()
{
  for (int d = 0; d < 3; d++) {
    nmesh_kk[d] = nmesh[d];
    foot_lo_kk[d] = foot_lo[d];
    foot_n_kk[d] = foot_n[d];
    fftn_kk[d] = fftn[d];
    for (int e = 0; e < 3; e++) mesh_vec_kk[d][e] = mesh_vec[d][e];
  }
  nfoot_kk = (int) nfoot;
  order_kk = order;
  nlower_kk = nlower;
}

/* ----------------------------------------------------------------------
   the extent of this rank's atoms in mesh coordinates
------------------------------------------------------------------------- */

template<class DeviceType>
int ComputeXRDFFTKokkos<DeviceType>::minmax_u(double *umin, double *umax)
{
  for (int d = 0; d < 3; d++) {
    umin[d] = 0.0;
    umax[d] = 0.0;
  }

  const int nlocal = atom->nlocal;
  if (nlocal == 0) return 0;

  XRDFFTBoundsFunctor<DeviceType> f;
  f.x = d_x;
  f.mask = d_mask;
  f.groupbit = groupbit;
  for (int d = 0; d < 3; d++)
    for (int e = 0; e < 3; e++) f.mv[d][e] = mesh_vec[d][e];

  XRDFFTBounds b;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal),f,b);

  if (!b.nany) return 0;

  for (int d = 0; d < 3; d++) {
    umin[d] = b.umin[d];
    umax[d] = b.umax[d];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::grow_density_own()
{
  const bigint n = MAX(nfoot,(bigint)1);
  if ((bigint) d_density_own.extent(0) >= n) return;

  k_density_own = FFT_DAT::tdual_FFT_SCALAR_1d("xrd/fft/kk:density_own",n);
  d_density_own = k_density_own.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   list the group atoms of each element together, as the host style does,
   with a scan per element so that the order within one is the atom order
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::bucket_atoms()
{
  const int nlocal = atom->nlocal;

  memory->grow(slot_start,nslot+1,"xrd/fft:slot_start");

  if (nlocal > bucket_maxatoms) {
    bucket_maxatoms = nlocal;
    k_slot_atoms = DAT::tdual_int_1d("xrd/fft/kk:slot_atoms",bucket_maxatoms);
    d_slot_atoms = k_slot_atoms.template view<DeviceType>();
  }

  slot_start[0] = 0;

  XRDFFTBucketFunctor<DeviceType> f;
  f.type = d_type;
  f.mask = d_mask;
  f.slot_of_type = d_slot_of_type;
  f.slot_atoms = d_slot_atoms;
  f.groupbit = groupbit;

  for (int s = 0; s < nslot; s++) {
    int total = 0;
    if (nlocal > 0) {
      f.slot = s;
      f.offset = slot_start[s];
      Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType>(0,nlocal),f,total);
    }
    slot_start[s+1] = slot_start[s] + total;
  }
}

/* ----------------------------------------------------------------------
   spread all group atoms of element slot s onto this rank's footprint
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::spread(int s)
{
  if (!nfoot) return;

  spread_lo = slot_start[s];
  spread_hi = slot_start[s+1];
  if (spread_hi <= spread_lo) return;

#ifdef LMP_KOKKOS_GPU
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTSpreadAtomic>(spread_lo,spread_hi),
                       *this);
  copymode = 0;
#else

  // atomics are slow on a CPU, so each thread takes a slice of the footprint
  // and only writes into that, which needs none

  copymode = 1;
  Kokkos::TeamPolicy<DeviceType, TagXRDFFTSpreadTiled> config(lmp->kokkos->nthreads,1);
  Kokkos::parallel_for(config,*this);
  copymode = 0;
#endif
}

/* ----------------------------------------------------------------------
   the mesh index of one atom along each dimension, as an index into the
   footprint of this rank.  returns nonzero when the stencil leaves it, which
   set_footprint() is supposed to make impossible
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int ComputeXRDFFTKokkos<DeviceType>::stencil_base(const int ii, double *delta, int *g) const
{
  const double xi = d_x(ii,0), yi = d_x(ii,1), zi = d_x(ii,2);
  int bad = 0;

  for (int d = 0; d < 3; d++) {
    const int n = nmesh_kk[d];
    const int flen = foot_n_kk[d];

    // fold into the diffraction cell.  atoms may sit far outside the box in a
    // non-periodic direction, and the cell may be much smaller than the box

    double u = xi*mesh_vec_kk[d][0] + yi*mesh_vec_kk[d][1] + zi*mesh_vec_kk[d][2];
    u -= n*Kokkos::floor(u/n);
    if ((u >= n) || (u < 0.0)) u = 0.0;

    const int ngrid = (int) Kokkos::floor(u + 0.5);
    delta[d] = ngrid - u;

    int gg = (ngrid + nlower_kk - foot_lo_kk[d]) % n;
    if (gg < 0) gg += n;

    // the stencil runs over consecutive footprint indices, wrapping only where
    // the footprint is the whole mesh.  anywhere else it has to fit

    if ((flen < n) && (gg + order_kk > flen)) bad = 1;
    g[d] = gg;
  }

  if (bad) {
    Kokkos::atomic_fetch_max(&d_overrun(),1);
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   the spreading window at the stencil points of one dimension, from its
   Chebyshev expansion in the offset of the atom from its nearest mesh point
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::stencil_rho(const int d, const double delta,
                                                  double *rho) const
{
  double t[KB_NCHEB];

  const double uu = 2.0*delta;
  const double uu2 = 2.0*uu;
  t[0] = 1.0;
  t[1] = uu;
  for (int i = 2; i < KB_NCHEB; i++) t[i] = uu2*t[i-1] - t[i-2];

  for (int j = 0; j < order_kk; j++) {
    const int off = (d*order_kk + j)*KB_NCHEB;
    double w = 0.0;
    for (int i = 0; i < KB_NCHEB; i++) w += d_kb_cheb(off+i)*t[i];
    rho[j] = w;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTSpreadAtomic, const int &b) const
{
  Kokkos::View<FFT_SCALAR*, Kokkos::LayoutRight, typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> >
    a_density_own = d_density_own;

  const int ii = d_slot_atoms(b);

  double delta[3];
  int g[3];
  if (stencil_base(ii,delta,g)) return;

  double rho0[MAXORDER], rho1[MAXORDER], rho2[MAXORDER];
  stencil_rho(0,delta[0],rho0);
  stencil_rho(1,delta[1],rho1);
  stencil_rho(2,delta[2],rho2);

  const int nx = nmesh_kk[0], ny = nmesh_kk[1], nz = nmesh_kk[2];
  const int fx = foot_n_kk[0], fy = foot_n_kk[1];

  for (int kk = 0; kk < order_kk; kk++) {
    int mz = g[2] + kk;
    if (mz >= nz) mz -= nz;
    const double z0 = rho2[kk];
    for (int jj = 0; jj < order_kk; jj++) {
      int my = g[1] + jj;
      if (my >= ny) my -= ny;
      const int row = (mz*fy + my)*fx;
      const double y0 = z0*rho1[jj];
      for (int i = 0; i < order_kk; i++) {
        int mx = g[0] + i;
        if (mx >= nx) mx -= nx;
        a_density_own(row + mx) += (FFT_SCALAR) (y0*rho0[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTSpreadTiled,
  typename Kokkos::TeamPolicy<DeviceType, TagXRDFFTSpreadTiled>::member_type dev) const
{
  const int tid = dev.league_rank();
  const int nthreads = dev.league_size();
  const int idelta = 1 + nfoot_kk/nthreads;
  const int ifrom = tid*idelta;
  const int ito = ((ifrom + idelta) > nfoot_kk) ? nfoot_kk : ifrom + idelta;
  if (ifrom >= ito) return;

  const int nx = nmesh_kk[0], ny = nmesh_kk[1], nz = nmesh_kk[2];
  const int fx = foot_n_kk[0], fy = foot_n_kk[1];

  for (int b = spread_lo; b < spread_hi; b++) {
    const int ii = d_slot_atoms(b);

    double delta[3];
    int g[3];
    if (stencil_base(ii,delta,g)) continue;

    // the footprint points this atom can reach, widened to the whole extent of
    // a dimension its stencil wraps in.  most atoms miss this thread's slice
    // entirely and are dropped before their window is ever evaluated

    int lo[3], hi[3];
    for (int d = 0; d < 3; d++) {
      if (g[d] + order_kk > foot_n_kk[d]) {
        lo[d] = 0;
        hi[d] = foot_n_kk[d] - 1;
      } else {
        lo[d] = g[d];
        hi[d] = g[d] + order_kk - 1;
      }
    }
    if (((hi[2]*fy + hi[1])*fx + hi[0] < ifrom) ||
        ((lo[2]*fy + lo[1])*fx + lo[0] >= ito)) continue;

    double rho0[MAXORDER], rho1[MAXORDER], rho2[MAXORDER];
    stencil_rho(0,delta[0],rho0);
    stencil_rho(1,delta[1],rho1);
    stencil_rho(2,delta[2],rho2);

    for (int kk = 0; kk < order_kk; kk++) {
      int mz = g[2] + kk;
      if (mz >= nz) mz -= nz;
      const double z0 = rho2[kk];
      for (int jj = 0; jj < order_kk; jj++) {
        int my = g[1] + jj;
        if (my >= ny) my -= ny;
        const int row = (mz*fy + my)*fx;
        const double y0 = z0*rho1[jj];
        for (int i = 0; i < order_kk; i++) {
          int mx = g[0] + i;
          if (mx >= nx) mx -= nx;
          const int il = row + mx;
          if ((il < ifrom) || (il >= ito)) continue;
          d_density_own(il) += (FFT_SCALAR) (y0*rho0[i]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTZeroFoot, const int &i) const
{
  d_density_own(i) = (FFT_SCALAR) 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTZeroSlab, const int &i) const
{
  d_density_slab(i) = (FFT_SCALAR) 0.0;
}

/* ----------------------------------------------------------------------
   the overlap of a footprint and a brick, walked in the order the host style
   walks it: the segments of a dimension end to end, slowest dimension first
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTPack, const int &m) const
{
  const int lenx = seglen_kk[0], leny = seglen_kk[1];
  const int nxy = lenx*leny;

  const int k = m/nxy;
  const int r = m - k*nxy;
  const int j = r/lenx;
  const int i = r - j*lenx;

  const int sz = (k < seg_kk[2][0][2]) ? seg_kk[2][0][0] + k
                                       : seg_kk[2][1][0] + (k - seg_kk[2][0][2]);
  const int sy = (j < seg_kk[1][0][2]) ? seg_kk[1][0][0] + j
                                       : seg_kk[1][1][0] + (j - seg_kk[1][0][2]);
  const int sx = (i < seg_kk[0][0][2]) ? seg_kk[0][0][0] + i
                                       : seg_kk[0][1][0] + (i - seg_kk[0][0][2]);

  d_sendbuf(bufoff + m) = d_density_own((sz*segsrc_n[1] + sy)*segsrc_n[0] + sx);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTUnpack, const int &m) const
{
  const int lenx = seglen_kk[0], leny = seglen_kk[1];
  const int nxy = lenx*leny;

  const int k = m/nxy;
  const int r = m - k*nxy;
  const int j = r/lenx;
  const int i = r - j*lenx;

  const int mz = (k < seg_kk[2][0][2]) ? seg_kk[2][0][1] + k
                                       : seg_kk[2][1][1] + (k - seg_kk[2][0][2]);
  const int my = (j < seg_kk[1][0][2]) ? seg_kk[1][0][1] + j
                                       : seg_kk[1][1][1] + (j - seg_kk[1][0][2]);
  const int mx = (i < seg_kk[0][0][2]) ? seg_kk[0][0][1] + i
                                       : seg_kk[0][1][1] + (i - seg_kk[0][0][2]);

  const int dst = ((mz - segdst_lo[2])*fftn_kk[1] + (my - segdst_lo[1]))*fftn_kk[0] +
                  (mx - segdst_lo[0]);

  d_density_slab(dst) += d_unpacksrc(bufoff + m);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTWork, const int &i) const
{
  d_work1(2*i) = d_density_slab(i);
  d_work1(2*i+1) = (FFT_SCALAR) 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeXRDFFTKokkos<DeviceType>::operator()(TagXRDFFTModes, const int &a) const
{
  const int idx = d_own_idx(a);
  const double w = d_own_deconv(a)*d_own_asf(slot_off + a);
  d_Fre(a) += w*(double) d_work1(2*idx);
  d_Fim(a) += w*(double) d_work1(2*idx+1);
}

/* ----------------------------------------------------------------------
   sum every rank's footprint into the bricks of the mesh, as the host style
   does, with the pieces packed and unpacked by kernels and the messages sent
   straight from device memory where MPI can read it
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::fold_reduce(int tag)
{
  if (all_full && is_slab) {
    FFT_SCALAR *sbuf, *rbuf;
    if (mpi_direct) {
      Kokkos::fence();
      sbuf = d_density_own.data();
      rbuf = d_density_slab.data();
    } else {
      k_density_own.modify_device();
      k_density_own.sync_host();
      sbuf = k_density_own.view_host().data();
      rbuf = k_density_slab.view_host().data();
    }

    MPI_Reduce_scatter(sbuf,rbuf,recvcounts,MPI_FFT_SCALAR,MPI_SUM,world);

    if (!mpi_direct) {
      k_density_slab.modify_host();
      k_density_slab.sync_device();
    }
    return;
  }

  if (nfft > 0) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTZeroSlab>(0,nfft),*this);
    copymode = 0;
  }

  if (nfoot > maxsend) {
    maxsend = nfoot;
    k_sendbuf = FFT_DAT::tdual_FFT_SCALAR_1d("xrd/fft/kk:sendbuf",maxsend);
    d_sendbuf = k_sendbuf.template view<DeviceType>();
  }

  // the ranks whose bricks this footprint reaches

  int ndest = 0;

  for (int qz = 0; qz < pgrid[2]; qz++)
    for (int qy = 0; qy < pgrid[1]; qy++)
      for (int qx = 0; qx < pgrid[0]; qx++) {
        const int q = (qz*pgrid[1] + qy)*pgrid[0] + qx;
        const bigint n = brick_count(q);
        if (n == 0) continue;
        destlist[ndest++] = q;
        scount[q] = (int) n;
      }

  int nout = 0;
  for (int i = 0; i < ndest; i++) {
    const int q = destlist[i];
    sstart[q] = nout;
    nout += scount[q];
  }

  // pack every piece before any of them is sent, so that a build without
  // GPU-aware MPI copies the whole buffer to the host once

  for (int d = 0; d < 3; d++) segsrc_n[d] = foot_n[d];

  for (int i = 0; i < ndest; i++) {
    const int q = destlist[i];
    for (int d = 0; d < 3; d++) {
      seg_kk[d][0][0] = seg_kk[d][0][1] = seg_kk[d][0][2] = 0;
      seg_kk[d][1][0] = seg_kk[d][1][1] = seg_kk[d][1][2] = 0;
      const int ns = segments(foot_lo[d],foot_n[d],nmesh[d],qgrid_lo(q,d),qgrid_hi(q,d),
                              seg_kk[d]);
      int len = 0;
      for (int t = 0; t < ns; t++) len += seg_kk[d][t][2];
      seglen_kk[d] = len;
    }
    bufoff = sstart[q];
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTPack>(0,scount[q]),*this);
    copymode = 0;
  }

  FFT_SCALAR *hsend = nullptr;
  if (ndest && mpi_direct) Kokkos::fence();
  if (ndest && !mpi_direct) {
    k_sendbuf.modify_device();
    k_sendbuf.sync_host();
    hsend = k_sendbuf.view_host().data();
  }

  int nsend = 0;
  int selfoff = -1;

  for (int i = 0; i < ndest; i++) {
    const int q = destlist[i];
    if (q == me) {
      selfoff = sstart[q];
      continue;
    }
    FFT_SCALAR *buf = mpi_direct ? &d_sendbuf.data()[sstart[q]] : &hsend[sstart[q]];
    MPI_Isend(buf,scount[q],MPI_FFT_SCALAR,q,tag,world,&requests[nsend++]);
  }

  for (int i = 0; i < ndest; i++) scount[destlist[i]] = 0;

  for (int d = 0; d < 3; d++) segdst_lo[d] = fftlo[d];

  // this rank's own contribution needs no message

  if (selfoff >= 0) {
    d_unpacksrc = d_sendbuf;
    bufoff = selfoff;
    unpack_kernel(me);
  }

  // count the messages to expect, then take them as they arrive

  int nrecv = 0;
  bigint biggest = 0;

  if (nfft > 0) {
    for (int q = 0; q < nprocs; q++) {
      if (q == me) continue;
      const int *lo = &foot_all[6*q];
      const int *ln = &foot_all[6*q+3];
      if ((ln[0] == 0) || (ln[1] == 0) || (ln[2] == 0)) continue;
      bigint n = 1;
      for (int d = 0; d < 3; d++) {
        int seg[2][3];
        const int ns = segments(lo[d],ln[d],nmesh[d],fftlo[d],ffthi[d],seg);
        int len = 0;
        for (int t = 0; t < ns; t++) len += seg[t][2];
        n *= len;
      }
      if (n == 0) continue;
      nrecv++;
      if (n > biggest) biggest = n;
    }
  }

  if (biggest > maxrecv) {
    maxrecv = biggest;
    k_recvbuf = FFT_DAT::tdual_FFT_SCALAR_1d("xrd/fft/kk:recvbuf",maxrecv);
    d_recvbuf = k_recvbuf.template view<DeviceType>();
  }

  for (int i = 0; i < nrecv; i++) {
    MPI_Status status;

    // every message is taken into the same buffer, so the unpack of the
    // previous one has to have read it before this one lands

    Kokkos::fence();

    FFT_SCALAR *buf = mpi_direct ? d_recvbuf.data() : k_recvbuf.view_host().data();
    MPI_Recv(buf,(int)maxrecv,MPI_FFT_SCALAR,MPI_ANY_SOURCE,tag,world,&status);

    if (!mpi_direct) {
      k_recvbuf.modify_host();
      k_recvbuf.sync_device();
    }
    d_unpacksrc = d_recvbuf;
    bufoff = 0;
    unpack_kernel(status.MPI_SOURCE);
  }

  if (nsend) MPI_Waitall(nsend,requests,MPI_STATUSES_IGNORE);
}

/* ----------------------------------------------------------------------
   add the piece of the mesh that came from rank q into this rank's brick
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::unpack_kernel(int q)
{
  const int *lo = &foot_all[6*q];
  const int *ln = &foot_all[6*q+3];
  if ((ln[0] == 0) || (ln[1] == 0) || (ln[2] == 0)) return;

  bigint total = 1;
  for (int d = 0; d < 3; d++) {
    seg_kk[d][0][0] = seg_kk[d][0][1] = seg_kk[d][0][2] = 0;
    seg_kk[d][1][0] = seg_kk[d][1][1] = seg_kk[d][1][2] = 0;
    const int ns = segments(lo[d],ln[d],nmesh[d],fftlo[d],ffthi[d],seg_kk[d]);
    int len = 0;
    for (int t = 0; t < ns; t++) len += seg_kk[d][t][2];
    seglen_kk[d] = len;
    total *= len;
  }

  if (total == 0) return;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTUnpack>(0,(int)total),*this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeXRDFFTKokkos<DeviceType>::compute_array()
{
  invoked_array = update->ntimestep;

  if (update_reciprocal()) refresh_scaling();

  double t0 = platform::walltime();

  bigint natoms = group->count(igroup);
  if (natoms == 0) natoms = 1;

  atomKK->sync(execution_space,datamask_read);
  d_x = atomKK->k_x.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  Kokkos::deep_copy(d_Fre,0.0);
  Kokkos::deep_copy(d_Fim,0.0);

  k_overrun.view_host()() = 0;
  k_overrun.modify_host();
  k_overrun.template sync<DeviceType>();

  // which part of the mesh this rank spreads into depends on where its atoms
  // are, so it is found once per invocation and shared by all the elements

  set_footprint();
  set_kernel_state();

  for (int s = 0; s < nslot; s++) {

    if (nfoot) {
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTZeroFoot>(0,nfoot_kk),*this);
      copymode = 0;
    }

    spread(s);

    fold_reduce(s);

    if (nfft > 0) {
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTWork>(0,nfft),*this);
      copymode = 0;

      fftkk->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::FORWARD);
    }

    // the FFT uses the exp(-i...) convention while compute xrd is defined with
    // exp(+i...), so this yields the complex conjugate of F.  only |F|^2 is
    // reported, so the conjugation is not observable.

    if (nown > 0) {
      slot_off = s*nown;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagXRDFFTModes>(0,nown),*this);
      copymode = 0;
    }
  }

  k_overrun.template modify<DeviceType>();
  k_overrun.sync_host();
  if (k_overrun.view_host()())
    error->one(FLERR,"Compute XRD/FFT: atom outside the mesh footprint of its MPI rank; "
               "this is an internal error, please report it");

  k_Fre.template modify<DeviceType>();
  k_Fre.sync_host();
  k_Fim.template modify<DeviceType>();
  k_Fim.sync_host();

  // every mode is owned by exactly one rank, so each row is contributed once

  auto h_Fre = k_Fre.view_host();
  auto h_Fim = k_Fim.view_host();

  for (int n = 0; n < size_array_rows; n++) Iloc[n] = 0.0;
  for (int a = 0; a < nown; a++) {
    const double re = h_Fre(a), im = h_Fim(a);
    Iloc[own_row[a]] = own_lp[a]*(re*re + im*im)/natoms;
  }

  MPI_Allreduce(Iloc,Iall,size_array_rows,MPI_DOUBLE,MPI_SUM,world);

  for (int n = 0; n < size_array_rows; n++) array[n][1] = Iall[n];

  if ((me == 0) && echo)
    utils::logmesg(lmp,"-----\nCompute XRD/FFT id:{} Elapsed time: {:0.2f} s\n-----\n",
                   id,platform::walltime()-t0);
}

/* ----------------------------------------------------------------------
   the mesh, the exchange buffers and the tables the kernels read are dual
   views, so on a device they are held on both sides of it
------------------------------------------------------------------------- */

template<class DeviceType>
double ComputeXRDFFTKokkos<DeviceType>::memory_usage()
{
  double bytes = ComputeXRDFFT::memory_usage();

  if (!setup_done || (execution_space != Device)) return bytes;

  bytes += (double)(nfoot + maxsend + maxrecv + 3*nfft) * sizeof(FFT_SCALAR);
  bytes += (double)(nown*(3+nslot) + 3*order*KB_NCHEB) * sizeof(double);
  bytes += (double)(nown + ntypes + bucket_maxatoms) * sizeof(int);
  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeXRDFFTKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeXRDFFTKokkos<LMPHostType>;
#endif
}
