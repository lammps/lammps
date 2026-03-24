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
   TIP4P PPPM Kokkos: subclass of PPPMKokkos that overrides the five
   virtual hooks that differ from standard PPPM:

     particle_map     - use M-site grid positions (from d_xM) instead of x
     make_rho         - same grid positions as particle_map
     fieldforce_ik    - split long-range force on O to O+(2×H) via alpha
     fieldforce_peratom - split per-atom energy/virial on O to O+(2×H)
     slabcorr         - TIP4P dipole (uses M for z-moment, real z for r2)

   Two-phase kernel strategy (performance and correctness):
     Phase 1 - every atom writes only to its own index (no races).
     Phase 2 - only O atoms write to their H atoms' indices.
               Since each H belongs to exactly one O, there is no
               write conflict between O iterations, so no atomics are
               needed anywhere in the hot path.

   The host-side find_M_host / tip4p_preprocess_host runs once per
   compute() call (after x2lamda if triclinic) to fill d_xM, d_ih1,
   d_ih2, which are then read-only on the device.
------------------------------------------------------------------------- */

#include "pppm_tip4p_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int OFFSET = 16384;
static constexpr FFT_SCALAR ZEROF = 0.0;
static constexpr double SMALL = 0.00001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PPPMTIP4PKokkos<DeviceType>::PPPMTIP4PKokkos(LAMMPS *lmp) : PPPMKokkos<DeviceType>(lmp)
{
  this->tip4pflag = 1;
  this->triclinic_support = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PPPMTIP4PKokkos<DeviceType>::~PPPMTIP4PKokkos() = default;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::init()
{
  if (this->force->newton == 0)
    this->error->all(FLERR, "Kspace style pppm/tip4p/kk requires newton on");
  PPPMKokkos<DeviceType>::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PPPMTIP4PKokkos<DeviceType>::memory_usage()
{
  double bytes = PPPMKokkos<DeviceType>::memory_usage();
  const int n = (int) k_xM.view_host().extent(0);
  if (n > 0) {
    bytes += (double) n * 3 * sizeof(KK_FLOAT);  // k_xM
    bytes += (double) 2 * n * sizeof(int);        // k_ih1 + k_ih2
  }
  return bytes;
}

/* ======================================================================
   Host-side M-site geometry: ported from PPPMTIP4P::find_M
   (private in the host class, so must be duplicated here)
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::find_M_host(int i, int &iH1, int &iH2, double *xM)
{
  double **x = this->atom->x;

  iH1 = this->atom->map(this->atom->tag[i] + 1);
  iH2 = this->atom->map(this->atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1)
    this->error->one(FLERR, "TIP4P hydrogen is missing");
  if (this->atom->type[iH1] != this->typeH || this->atom->type[iH2] != this->typeH)
    this->error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

  if (this->triclinic) {

    int *sametag = this->atom->sametag;
    double xo[3], xh1[3], xh2[3], xm[3];
    const int nlocal = this->atom->nlocal;

    for (int d = 0; d < 3; ++d) {
      xo[d]  = x[i][d];
      xh1[d] = x[iH1][d];
      xh2[d] = x[iH2][d];
    }

    // local atoms are in lamda coords, ghosts are not
    if (i   < nlocal) this->domain->lamda2x(x[i],   xo);
    if (iH1 < nlocal) this->domain->lamda2x(x[iH1], xh1);
    if (iH2 < nlocal) this->domain->lamda2x(x[iH2], xh2);

    // find closest periodic image of H1 to O
    {
      double delx = xo[0] - xh1[0];
      double dely = xo[1] - xh1[1];
      double delz = xo[2] - xh1[2];
      double rsqmin = delx*delx + dely*dely + delz*delz;
      int closest = iH1;

      while (sametag[iH1] >= 0) {
        iH1 = sametag[iH1];
        delx = xo[0] - x[iH1][0];
        dely = xo[1] - x[iH1][1];
        delz = xo[2] - x[iH1][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < rsqmin) {
          rsqmin = rsq;
          closest = iH1;
          xh1[0] = x[iH1][0];
          xh1[1] = x[iH1][1];
          xh1[2] = x[iH1][2];
        }
      }
      iH1 = closest;
    }

    // find closest periodic image of H2 to O
    {
      double delx = xo[0] - xh2[0];
      double dely = xo[1] - xh2[1];
      double delz = xo[2] - xh2[2];
      double rsqmin = delx*delx + dely*dely + delz*delz;
      int closest = iH2;

      while (sametag[iH2] >= 0) {
        iH2 = sametag[iH2];
        delx = xo[0] - x[iH2][0];
        dely = xo[1] - x[iH2][1];
        delz = xo[2] - x[iH2][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < rsqmin) {
          rsqmin = rsq;
          closest = iH2;
          xh2[0] = x[iH2][0];
          xh2[1] = x[iH2][1];
          xh2[2] = x[iH2][2];
        }
      }
      iH2 = closest;
    }

    // compute M in real coords, then convert to lamda for PPPM
    xm[0] = xo[0] + this->alpha * 0.5 * ((xh1[0] - xo[0]) + (xh2[0] - xo[0]));
    xm[1] = xo[1] + this->alpha * 0.5 * ((xh1[1] - xo[1]) + (xh2[1] - xo[1]));
    xm[2] = xo[2] + this->alpha * 0.5 * ((xh1[2] - xo[2]) + (xh2[2] - xo[2]));
    this->domain->x2lamda(xm, xM);

  } else {

    iH1 = this->domain->closest_image(i, iH1);
    iH2 = this->domain->closest_image(i, iH2);

    xM[0] = x[i][0] + this->alpha * 0.5 * ((x[iH1][0] - x[i][0]) + (x[iH2][0] - x[i][0]));
    xM[1] = x[i][1] + this->alpha * 0.5 * ((x[iH1][1] - x[i][1]) + (x[iH2][1] - x[i][1]));
    xM[2] = x[i][2] + this->alpha * 0.5 * ((x[iH1][2] - x[i][2]) + (x[iH2][2] - x[i][2]));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::tip4p_preprocess_host()
{
  alpha_kk = static_cast<KK_FLOAT>(this->alpha);

  const int nlocal     = this->atom->nlocal;
  const int nmax_atom  = this->atom->nmax;

  // grow views if needed (nmax grows monotonically, so this is rare)
  if (k_xM.view_host().extent(0) < (size_t) nmax_atom) {
    k_xM = DAT::tdual_kkfloat_1d_3("pppm/tip4p/kk:xM", nmax_atom);
    k_ih1 = DAT::tdual_int_1d("pppm/tip4p/kk:ih1", nmax_atom);
    k_ih2 = DAT::tdual_int_1d("pppm/tip4p/kk:ih2", nmax_atom);
  }

  // Ensure host copies of x and type are current.
  // PPPMKokkos::compute() calls atomKK->sync(execution_space, datamask_read)
  // which syncs host→device when execution_space is GPU.  For GPU builds the
  // host C-arrays (atom->x, atom->type) remain the source and are current.
  // For safety against fixes that modify positions on device, sync to host.
  this->atomKK->k_x.sync_host();
  this->atomKK->k_type.sync_host();

  double **x   = this->atom->x;
  int    *type = this->atom->type;

  auto h_xM  = k_xM.view_host();
  auto h_ih1 = k_ih1.view_host();
  auto h_ih2 = k_ih2.view_host();

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == this->typeO) {
      int iH1, iH2;
      double xM[3];
      find_M_host(i, iH1, iH2, xM);
      h_xM(i, 0) = static_cast<KK_FLOAT>(xM[0]);
      h_xM(i, 1) = static_cast<KK_FLOAT>(xM[1]);
      h_xM(i, 2) = static_cast<KK_FLOAT>(xM[2]);
      h_ih1(i) = iH1;
      h_ih2(i) = iH2;
    } else {
      h_xM(i, 0) = static_cast<KK_FLOAT>(x[i][0]);
      h_xM(i, 1) = static_cast<KK_FLOAT>(x[i][1]);
      h_xM(i, 2) = static_cast<KK_FLOAT>(x[i][2]);
      h_ih1(i) = -1;
      h_ih2(i) = -1;
    }
  }

  k_xM.modify_host();
  k_xM.template sync<DeviceType>();
  k_ih1.modify_host();
  k_ih1.template sync<DeviceType>();
  k_ih2.modify_host();
  k_ih2.template sync<DeviceType>();

  d_xM   = k_xM.template view<DeviceType>();
  d_ih1  = k_ih1.template view<DeviceType>();
  d_ih2  = k_ih2.template view<DeviceType>();
  d_type = this->atomKK->k_type.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::pp_pre_particle_map()
{
  tip4p_preprocess_host();
}

/* ======================================================================
   particle_map: map each atom's "effective" position (M-site for O,
   actual position for H/others) onto the grid.
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::particle_map()
{
  const int nlocal = this->atomKK->nlocal;

  this->shift_kk    = static_cast<KK_FLOAT>(this->shift);
  this->delxinv_kk  = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk  = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk  = static_cast<KK_FLOAT>(this->delzinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->k_flag.view_host()() = 0;
  this->k_flag.modify_host();
  this->k_flag.template sync<DeviceType>();

  if (!std::isfinite(this->boxlo[0]) || !std::isfinite(this->boxlo[1]) ||
      !std::isfinite(this->boxlo[2]))
    this->error->one(FLERR,
        "Non-numeric box dimensions - simulation unstable" + utils::errorurl(6));

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_particle_map>(0, nlocal), *this);
  this->copymode = 0;

  this->k_flag.template modify<DeviceType>();
  this->k_flag.sync_host();
  if (this->k_flag.view_host()())
    this->error->one(FLERR, Error::NOLASTLINE,
        "Out of range atoms - cannot compute PPPM" + utils::errorurl(4));
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPM_particle_map, const int &i) const
{
  const int nx = static_cast<int>(
      (d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk + this->shift_kk) - OFFSET;
  const int ny = static_cast<int>(
      (d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk + this->shift_kk) - OFFSET;
  const int nz = static_cast<int>(
      (d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk + this->shift_kk) - OFFSET;

  this->d_part2grid(i, 0) = nx;
  this->d_part2grid(i, 1) = ny;
  this->d_part2grid(i, 2) = nz;

  if (nx + this->nlower < this->nxlo_out || nx + this->nupper > this->nxhi_out ||
      ny + this->nlower < this->nylo_out || ny + this->nupper > this->nyhi_out ||
      nz + this->nlower < this->nzlo_out || nz + this->nupper > this->nzhi_out)
    this->k_flag.template view<DeviceType>()() = 1;
}

/* ======================================================================
   make_rho: assign charge to grid using M-site positions.
   The zero-pass kernel is dispatched to the base class (no TIP4P logic).
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::make_rho()
{
  this->numz_out = this->nzhi_out - this->nzlo_out + 1;
  this->numy_out = this->nyhi_out - this->nylo_out + 1;
  this->numx_out = this->nxhi_out - this->nxlo_out + 1;
  const int inum_out = this->numz_out * this->numy_out * this->numx_out;

  this->shiftone_kk  = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk   = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk   = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk   = static_cast<KK_FLOAT>(this->delzinv);
  this->delvolinv_kk = static_cast<KK_FLOAT>(this->delvolinv);
  this->boxlo_kk[0]  = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1]  = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2]  = static_cast<KK_FLOAT>(this->boxlo[2]);

  // zero the density brick using the base class functor (no TIP4P specifics)
  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_make_rho_zero>(0, inum_out),
                       static_cast<PPPMKokkos<DeviceType> &>(*this));
  this->copymode = 0;

  this->nlocal = this->atomKK->nlocal;

#ifdef LMP_KOKKOS_GPU
  this->copymode = 1;
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_make_rho_atomic>(0, this->nlocal), *this);
  this->copymode = 0;
#else
  this->ix = this->nxhi_out - this->nxlo_out + 1;
  this->iy = this->nyhi_out - this->nylo_out + 1;
  this->copymode = 1;
  Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho> config(
      this->lmp->kokkos->nthreads, 1);
  Kokkos::parallel_for(config, *this);
  this->copymode = 0;
#endif
}

/* GPU path: atomic scatter into density_brick via M-site positions */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_make_rho_atomic,
                                             const int &i) const
{
  Kokkos::View<FFT_SCALAR ***, Kokkos::LayoutRight,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      a_density_brick = this->d_density_brick;

  int nx = this->d_part2grid(i, 0);
  int ny = this->d_part2grid(i, 1);
  int nz = this->d_part2grid(i, 2);

  const FFT_SCALAR dx = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(nx) + this->shiftone_kk -
      (d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk);
  const FFT_SCALAR dy = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(ny) + this->shiftone_kk -
      (d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk);
  const FFT_SCALAR dz = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(nz) + this->shiftone_kk -
      (d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk);

  nz -= this->nzlo_out;
  ny -= this->nylo_out;
  nx -= this->nxlo_out;

  this->compute_rho1d(i, dx, dy, dz);

  const FFT_SCALAR z0 = static_cast<FFT_SCALAR>(this->delvolinv_kk * this->q[i]);
  for (int n = this->nlower; n <= this->nupper; n++) {
    const int mz = n + nz;
    const FFT_SCALAR y0 = z0 * this->d_rho1d(i, n + this->order / 2, 2);
    for (int m = this->nlower; m <= this->nupper; m++) {
      const int my = m + ny;
      const FFT_SCALAR x0 = y0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (int l = this->nlower; l <= this->nupper; l++) {
        a_density_brick(mz, my, l + nx) += x0 * this->d_rho1d(i, l + this->order / 2, 0);
      }
    }
  }
}

/* CPU path: team-based (thread-per-grid-slice) scatter using M-site positions */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(
    TagPPPMTIP4P_make_rho,
    typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type dev) const
{
  const int tid      = dev.league_rank();
  const int nthreads = dev.league_size();
  const int idelta   = 1 + this->ngrid / nthreads;
  const int ifrom    = tid * idelta;
  const int ito      = ((ifrom + idelta) > this->ngrid) ? this->ngrid : ifrom + idelta;

  for (int i = 0; i < this->nlocal; i++) {

    int nx = this->d_part2grid(i, 0);
    int ny = this->d_part2grid(i, 1);
    int nz = this->d_part2grid(i, 2);

    // early exit if atom's stencil is entirely outside this thread's grid slice
    if (((nz + this->nlower - this->nzlo_out) * this->ix * this->iy >= ito) ||
        ((nz + this->nupper - this->nzlo_out + 1) * this->ix * this->iy < ifrom))
      continue;

    const FFT_SCALAR dx = static_cast<FFT_SCALAR>(
        static_cast<KK_FLOAT>(nx) + this->shiftone_kk -
        (d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk);
    const FFT_SCALAR dy = static_cast<FFT_SCALAR>(
        static_cast<KK_FLOAT>(ny) + this->shiftone_kk -
        (d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk);
    const FFT_SCALAR dz = static_cast<FFT_SCALAR>(
        static_cast<KK_FLOAT>(nz) + this->shiftone_kk -
        (d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk);

    nz -= this->nzlo_out;
    ny -= this->nylo_out;
    nx -= this->nxlo_out;

    this->compute_rho1d(i, dx, dy, dz);

    const FFT_SCALAR z0 = static_cast<FFT_SCALAR>(this->delvolinv_kk * this->q[i]);
    for (int n = this->nlower; n <= this->nupper; n++) {
      const int mz = n + nz;
      const int in = mz * this->ix * this->iy;
      const FFT_SCALAR y0 = z0 * this->d_rho1d(i, n + this->order / 2, 2);
      for (int m = this->nlower; m <= this->nupper; m++) {
        const int my = m + ny;
        const int im = in + my * this->ix;
        const FFT_SCALAR x0 = y0 * this->d_rho1d(i, m + this->order / 2, 1);
        for (int l = this->nlower; l <= this->nupper; l++) {
          const int mx = l + nx;
          const int il = im + mx;
          if (il >= ito) break;
          if (il < ifrom) continue;
          this->d_density_brick(mz, my, mx) += x0 * this->d_rho1d(i, l + this->order / 2, 0);
        }
      }
    }
  }
}

/* ======================================================================
   fieldforce_ik: two-phase kernel.

   Phase 1 (all atoms): each atom reads its stencil from d_rho1d (filled
     during make_rho using M-site positions) and writes force to its OWN
     index only.  No writes to other atoms' indices → no race conditions.

   Phase 2 (O atoms only): O atoms re-read their stencil and scatter
     0.5*alpha fraction to each H atom.  Since each H belongs to exactly
     one O molecule, no two O iterations write to the same H → no race.
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_ik()
{
  const int nlocal = this->atomKK->nlocal;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  this->copymode = 1;
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_fieldforce_ik_phase1>(0, nlocal), *this);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_fieldforce_ik_phase2>(0, nlocal), *this);
  this->copymode = 0;
}

/* phase1: each atom writes force to its own index */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_fieldforce_ik_phase1,
                                             const int &i) const
{
  int nx = this->d_part2grid(i, 0) - this->nxlo_out;
  int ny = this->d_part2grid(i, 1) - this->nylo_out;
  int nz = this->d_part2grid(i, 2) - this->nzlo_out;

  FFT_SCALAR ekx = ZEROF, eky = ZEROF, ekz = ZEROF;
  for (int n = this->nlower; n <= this->nupper; n++) {
    const int mz = n + nz;
    const FFT_SCALAR z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (int m = this->nlower; m <= this->nupper; m++) {
      const int my = m + ny;
      const FFT_SCALAR y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (int l = this->nlower; l <= this->nupper; l++) {
        const FFT_SCALAR x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
        ekx -= x0 * this->d_vdx_brick(mz, my, l + nx);
        eky -= x0 * this->d_vdy_brick(mz, my, l + nx);
        ekz -= x0 * this->d_vdz_brick(mz, my, l + nx);
      }
    }
  }

  const KK_FLOAT qfactor = this->qscale_kk * static_cast<KK_FLOAT>(this->q[i]);

  if (d_type(i) != this->typeO) {
    // H and other atoms: full force on self
    this->f(i, 0) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekx));
    this->f(i, 1) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(eky));
    if (this->slabflag != 2)
      this->f(i, 2) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekz));
  } else {
    // O atom: apply (1 - alpha) fraction to O itself
    const KK_FLOAT one_minus_alpha = KK_FLOAT(1) - alpha_kk;
    this->f(i, 0) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekx) * one_minus_alpha);
    this->f(i, 1) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(eky) * one_minus_alpha);
    if (this->slabflag != 2)
      this->f(i, 2) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekz) * one_minus_alpha);
  }
}

/* phase2: O atoms scatter 0.5*alpha fraction to each H (no races) */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_fieldforce_ik_phase2,
                                             const int &i) const
{
  if (d_type(i) != this->typeO) return;

  int nx = this->d_part2grid(i, 0) - this->nxlo_out;
  int ny = this->d_part2grid(i, 1) - this->nylo_out;
  int nz = this->d_part2grid(i, 2) - this->nzlo_out;

  FFT_SCALAR ekx = ZEROF, eky = ZEROF, ekz = ZEROF;
  for (int n = this->nlower; n <= this->nupper; n++) {
    const int mz = n + nz;
    const FFT_SCALAR z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (int m = this->nlower; m <= this->nupper; m++) {
      const int my = m + ny;
      const FFT_SCALAR y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (int l = this->nlower; l <= this->nupper; l++) {
        const FFT_SCALAR x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
        ekx -= x0 * this->d_vdx_brick(mz, my, l + nx);
        eky -= x0 * this->d_vdy_brick(mz, my, l + nx);
        ekz -= x0 * this->d_vdz_brick(mz, my, l + nx);
      }
    }
  }

  const KK_FLOAT half_alpha = KK_FLOAT(0.5) * alpha_kk;
  const KK_FLOAT qfactor    = this->qscale_kk * static_cast<KK_FLOAT>(this->q[i]);
  const KK_FLOAT hfx = qfactor * static_cast<KK_FLOAT>(ekx) * half_alpha;
  const KK_FLOAT hfy = qfactor * static_cast<KK_FLOAT>(eky) * half_alpha;
  const KK_FLOAT hfz = qfactor * static_cast<KK_FLOAT>(ekz) * half_alpha;

  const int iH1 = d_ih1(i);
  const int iH2 = d_ih2(i);

  // no race: each H is the child of exactly one O
  this->f(iH1, 0) += static_cast<KK_ACC_FLOAT>(hfx);
  this->f(iH1, 1) += static_cast<KK_ACC_FLOAT>(hfy);
  this->f(iH2, 0) += static_cast<KK_ACC_FLOAT>(hfx);
  this->f(iH2, 1) += static_cast<KK_ACC_FLOAT>(hfy);
  if (this->slabflag != 2) {
    this->f(iH1, 2) += static_cast<KK_ACC_FLOAT>(hfz);
    this->f(iH2, 2) += static_cast<KK_ACC_FLOAT>(hfz);
  }
}

/* ======================================================================
   fieldforce_peratom: two-phase kernel.

   Phase 1 (all atoms, dispatched as TagPPPM_fieldforce_peratom):
     Each atom interpolates u_brick / v*_brick at its M-site position and
     writes per-atom energy / virial to its OWN index:
       non-O: full q[i]*u contribution
       O:     (1-alpha)*q[i]*u contribution

   Phase 2 (O atoms only, TagPPPMTIP4P_fieldforce_peratom_phase2):
     O atoms re-read their already-set d_rho1d and scatter 0.5*alpha
     contribution to each H.  No races — same argument as fieldforce_ik.
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_peratom()
{
  const int nlocal = this->atomKK->nlocal;

  this->shiftone_kk  = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk   = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk   = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk   = static_cast<KK_FLOAT>(this->delzinv);
  this->boxlo_kk[0]  = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1]  = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2]  = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->copymode = 1;
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPM_fieldforce_peratom>(0, nlocal), *this);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_fieldforce_peratom_phase2>(0, nlocal), *this);
  this->copymode = 0;
}

/* phase1: all atoms, own eatom/vatom only */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPM_fieldforce_peratom,
                                             const int &i) const
{
  int nx = this->d_part2grid(i, 0);
  int ny = this->d_part2grid(i, 1);
  int nz = this->d_part2grid(i, 2);

  // fractional offset within the stencil (uses M-site for O, actual pos for H)
  const FFT_SCALAR dx = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(nx) + this->shiftone_kk -
      (d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk);
  const FFT_SCALAR dy = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(ny) + this->shiftone_kk -
      (d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk);
  const FFT_SCALAR dz = static_cast<FFT_SCALAR>(
      static_cast<KK_FLOAT>(nz) + this->shiftone_kk -
      (d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk);

  nz -= this->nzlo_out;
  ny -= this->nylo_out;
  nx -= this->nxlo_out;

  this->compute_rho1d(i, dx, dy, dz);

  FFT_SCALAR u_pa = ZEROF, v0 = ZEROF, v1 = ZEROF, v2 = ZEROF;
  FFT_SCALAR v3 = ZEROF, v4 = ZEROF, v5 = ZEROF;

  for (int n = this->nlower; n <= this->nupper; n++) {
    const int mz = n + nz;
    const FFT_SCALAR z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (int m = this->nlower; m <= this->nupper; m++) {
      const int my = m + ny;
      const FFT_SCALAR y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (int l = this->nlower; l <= this->nupper; l++) {
        const int mx = l + nx;
        const FFT_SCALAR x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
        if (this->eflag_atom) u_pa += x0 * this->d_u_brick(mz, my, mx);
        if (this->vflag_atom) {
          v0 += x0 * this->d_v0_brick(mz, my, mx);
          v1 += x0 * this->d_v1_brick(mz, my, mx);
          v2 += x0 * this->d_v2_brick(mz, my, mx);
          v3 += x0 * this->d_v3_brick(mz, my, mx);
          v4 += x0 * this->d_v4_brick(mz, my, mx);
          v5 += x0 * this->d_v5_brick(mz, my, mx);
        }
      }
    }
  }

  // per-atom weight: (1-alpha) for O, full for others
  const bool is_O = (d_type(i) == this->typeO);
  const KK_FLOAT wt = is_O ? (KK_FLOAT(1) - alpha_kk) : KK_FLOAT(1);
  const KK_ACC_FLOAT q_wt = static_cast<KK_ACC_FLOAT>(this->q[i] * wt);

  if (this->eflag_atom)
    this->d_eatom[i] += q_wt * static_cast<KK_ACC_FLOAT>(u_pa);
  if (this->vflag_atom) {
    this->d_vatom(i, 0) += q_wt * static_cast<KK_ACC_FLOAT>(v0);
    this->d_vatom(i, 1) += q_wt * static_cast<KK_ACC_FLOAT>(v1);
    this->d_vatom(i, 2) += q_wt * static_cast<KK_ACC_FLOAT>(v2);
    this->d_vatom(i, 3) += q_wt * static_cast<KK_ACC_FLOAT>(v3);
    this->d_vatom(i, 4) += q_wt * static_cast<KK_ACC_FLOAT>(v4);
    this->d_vatom(i, 5) += q_wt * static_cast<KK_ACC_FLOAT>(v5);
  }
}

/* phase2: O atoms scatter 0.5*alpha contribution to each H (no races) */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_fieldforce_peratom_phase2,
                                             const int &i) const
{
  if (d_type(i) != this->typeO) return;

  // d_rho1d[i] was filled in phase1 and is still valid (same M-site stencil)
  int nx = this->d_part2grid(i, 0) - this->nxlo_out;
  int ny = this->d_part2grid(i, 1) - this->nylo_out;
  int nz = this->d_part2grid(i, 2) - this->nzlo_out;

  FFT_SCALAR u_pa = ZEROF, v0 = ZEROF, v1 = ZEROF, v2 = ZEROF;
  FFT_SCALAR v3 = ZEROF, v4 = ZEROF, v5 = ZEROF;

  for (int n = this->nlower; n <= this->nupper; n++) {
    const int mz = n + nz;
    const FFT_SCALAR z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (int m = this->nlower; m <= this->nupper; m++) {
      const int my = m + ny;
      const FFT_SCALAR y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (int l = this->nlower; l <= this->nupper; l++) {
        const int mx = l + nx;
        const FFT_SCALAR x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
        if (this->eflag_atom) u_pa += x0 * this->d_u_brick(mz, my, mx);
        if (this->vflag_atom) {
          v0 += x0 * this->d_v0_brick(mz, my, mx);
          v1 += x0 * this->d_v1_brick(mz, my, mx);
          v2 += x0 * this->d_v2_brick(mz, my, mx);
          v3 += x0 * this->d_v3_brick(mz, my, mx);
          v4 += x0 * this->d_v4_brick(mz, my, mx);
          v5 += x0 * this->d_v5_brick(mz, my, mx);
        }
      }
    }
  }

  const KK_ACC_FLOAT q_halpha = static_cast<KK_ACC_FLOAT>(
      this->q[i] * KK_FLOAT(0.5) * alpha_kk);

  const int iH1 = d_ih1(i);
  const int iH2 = d_ih2(i);

  if (this->eflag_atom) {
    const KK_ACC_FLOAT contrib = q_halpha * static_cast<KK_ACC_FLOAT>(u_pa);
    this->d_eatom[iH1] += contrib;
    this->d_eatom[iH2] += contrib;
  }
  if (this->vflag_atom) {
    this->d_vatom(iH1, 0) += q_halpha * static_cast<KK_ACC_FLOAT>(v0);
    this->d_vatom(iH1, 1) += q_halpha * static_cast<KK_ACC_FLOAT>(v1);
    this->d_vatom(iH1, 2) += q_halpha * static_cast<KK_ACC_FLOAT>(v2);
    this->d_vatom(iH1, 3) += q_halpha * static_cast<KK_ACC_FLOAT>(v3);
    this->d_vatom(iH1, 4) += q_halpha * static_cast<KK_ACC_FLOAT>(v4);
    this->d_vatom(iH1, 5) += q_halpha * static_cast<KK_ACC_FLOAT>(v5);
    this->d_vatom(iH2, 0) += q_halpha * static_cast<KK_ACC_FLOAT>(v0);
    this->d_vatom(iH2, 1) += q_halpha * static_cast<KK_ACC_FLOAT>(v1);
    this->d_vatom(iH2, 2) += q_halpha * static_cast<KK_ACC_FLOAT>(v2);
    this->d_vatom(iH2, 3) += q_halpha * static_cast<KK_ACC_FLOAT>(v3);
    this->d_vatom(iH2, 4) += q_halpha * static_cast<KK_ACC_FLOAT>(v4);
    this->d_vatom(iH2, 5) += q_halpha * static_cast<KK_ACC_FLOAT>(v5);
  }
}

/* ======================================================================
   slabcorr: TIP4P dipole correction.

   Key differences from base PPPM slabcorr:
     - dipole uses M-site z for O atoms  (dipole_all)
     - dipole_r2 uses actual atom z      (matches host PPPMTIP4P::slabcorr)
     - per-atom slab energy uses actual atom z (host code)
     - force on O split to O and 2×H via alpha (two-phase for safety)
     - energy added to energy_1 (matching host PPPMTIP4P::slabcorr line 526)
   ====================================================================== */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::slabcorr()
{
  this->zprd_slab = this->domain->zprd * this->slab_volfactor;
  const int nlocal = this->atomKK->nlocal;

  // dipole sum using M-site z for O atoms
  double dipole = 0.0;
  this->copymode = 1;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr1>(0, nlocal), *this, dipole);
  this->copymode = 0;
  MPI_Allreduce(&dipole, &this->dipole_all, 1, MPI_DOUBLE, MPI_SUM, this->world);

  // dipole_r2 uses actual atom z (same as host)
  this->dipole_r2 = 0.0;
  if (this->eflag_atom || fabs(this->qsum) > SMALL) {
    this->copymode = 1;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr2>(0, nlocal), *this,
        this->dipole_r2);
    this->copymode = 0;
    double tmp;
    MPI_Allreduce(&this->dipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, this->world);
    this->dipole_r2 = tmp;
  }

  const double e_slabcorr =
      MY_2PI * (this->dipole_all * this->dipole_all
               - this->qsum * this->dipole_r2
               - this->qsum * this->qsum * this->zprd_slab * this->zprd_slab / 12.0)
      / this->volume;
  this->qscale    = this->force->qqrd2e * this->scale;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  // match host: PPPMTIP4P::slabcorr() uses energy_1 (not energy)
  if (this->eflag_global) this->energy_1 += this->qscale * e_slabcorr;

  if (this->eflag_atom) {
    this->efact = this->qscale * MY_2PI / this->volume;
    this->copymode = 1;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr3>(0, nlocal), *this);
    this->copymode = 0;
  }

  this->ffact = this->qscale * (-4.0 * MY_PI / this->volume);

  // two-phase z-force: phase1 = own force, phase2 = H-atom scatter (no races)
  this->copymode = 1;
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr4_phase1>(0, nlocal), *this);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr4_phase2>(0, nlocal), *this);
  this->copymode = 0;
}

/* dipole: M-site z for O, actual z for others */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr1,
                                             const int &i, double &dipole) const
{
  const double z = (d_type(i) == this->typeO)
      ? static_cast<double>(d_xM(i, 2))
      : static_cast<double>(this->x(i, 2));
  dipole += static_cast<double>(this->q[i]) * z;
}

/* dipole_r2: actual atom z (host code uses x[i][2] not xM[2]) */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr2,
                                             const int &i, double &dipole_r2) const
{
  const double z = static_cast<double>(this->x(i, 2));
  dipole_r2 += static_cast<double>(this->q[i]) * z * z;
}

/* per-atom slab energy: actual atom z (host uses x[i][2]) */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr3, const int &i) const
{
  const double z_i = static_cast<double>(this->x(i, 2));
  const double q_i = static_cast<double>(this->q[i]);
  this->d_eatom[i] += static_cast<KK_ACC_FLOAT>(
      this->efact * q_i *
      (z_i * this->dipole_all
       - 0.5 * (this->dipole_r2 + this->qsum * z_i * z_i)
       - this->qsum * this->zprd_slab * this->zprd_slab / 12.0));
}

/* z-force phase1: each atom writes only its own force */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr4_phase1,
                                             const int &i) const
{
  const double z_i     = static_cast<double>(this->x(i, 2));
  const double fzi_all = this->ffact * static_cast<double>(this->q[i])
                         * (this->dipole_all - this->qsum * z_i);

  const KK_FLOAT wt = (d_type(i) == this->typeO)
      ? (KK_FLOAT(1) - alpha_kk)
      : KK_FLOAT(1);
  this->f(i, 2) += static_cast<KK_ACC_FLOAT>(fzi_all * static_cast<double>(wt));
}

/* z-force phase2: O atoms scatter 0.5*alpha to each H (no races) */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr4_phase2,
                                             const int &i) const
{
  if (d_type(i) != this->typeO) return;

  const double z_i      = static_cast<double>(this->x(i, 2));
  const double fzi_corr = this->ffact * static_cast<double>(this->q[i])
                          * (this->dipole_all - this->qsum * z_i);
  const KK_ACC_FLOAT hf = static_cast<KK_ACC_FLOAT>(
      fzi_corr * static_cast<double>(KK_FLOAT(0.5) * alpha_kk));

  const int iH1 = d_ih1(i);
  const int iH2 = d_ih2(i);
  this->f(iH1, 2) += hf;
  this->f(iH2, 2) += hf;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PPPMTIP4PKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PPPMTIP4PKokkos<LMPHostType>;
#endif
}
