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
    this->error->all(FLERR,"Kspace style pppm/tip4p/kk requires newton on");
  PPPMKokkos<DeviceType>::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PPPMTIP4PKokkos<DeviceType>::memory_usage()
{
  double bytes = PPPMKokkos<DeviceType>::memory_usage();
  const int n = (int) k_xM.extent(0);
  if (n > 0) {
    bytes += (double) n * 3 * sizeof(KK_FLOAT);
    bytes += (double) 2 * n * sizeof(int);
  }
  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::find_M_host(int i, int &iH1, int &iH2, double *xM)
{
  double **x = this->atom->x;

  iH1 = this->atom->map(this->atom->tag[i] + 1);
  iH2 = this->atom->map(this->atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) this->error->one(FLERR,"TIP4P hydrogen is missing");
  if (this->atom->type[iH1] != this->typeH || this->atom->type[iH2] != this->typeH)
    this->error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  if (this->triclinic) {

    int *sametag = this->atom->sametag;
    double xo[3], xh1[3], xh2[3], xm[3];
    const int nlocal = this->atom->nlocal;

    for (int ii = 0; ii < 3; ++ii) {
      xo[ii] = x[i][ii];
      xh1[ii] = x[iH1][ii];
      xh2[ii] = x[iH2][ii];
    }

    if (i < nlocal) this->domain->lamda2x(x[i], xo);
    if (iH1 < nlocal) this->domain->lamda2x(x[iH1], xh1);
    if (iH2 < nlocal) this->domain->lamda2x(x[iH2], xh2);

    double delx = xo[0] - xh1[0];
    double dely = xo[1] - xh1[1];
    double delz = xo[2] - xh1[2];
    double rsqmin = delx * delx + dely * dely + delz * delz;
    double rsq;
    int closest = iH1;

    while (sametag[iH1] >= 0) {
      iH1 = sametag[iH1];
      delx = xo[0] - x[iH1][0];
      dely = xo[1] - x[iH1][1];
      delz = xo[2] - x[iH1][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH1;
        xh1[0] = x[iH1][0];
        xh1[1] = x[iH1][1];
        xh1[2] = x[iH1][2];
      }
    }
    iH1 = closest;

    closest = iH2;
    delx = xo[0] - xh2[0];
    dely = xo[1] - xh2[1];
    delz = xo[2] - xh2[2];
    rsqmin = delx * delx + dely * dely + delz * delz;

    while (sametag[iH2] >= 0) {
      iH2 = sametag[iH2];
      delx = xo[0] - x[iH2][0];
      dely = xo[1] - x[iH2][1];
      delz = xo[2] - x[iH2][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH2;
        xh2[0] = x[iH2][0];
        xh2[1] = x[iH2][1];
        xh2[2] = x[iH2][2];
      }
    }
    iH2 = closest;

    double delx1 = xh1[0] - xo[0];
    double dely1 = xh1[1] - xo[1];
    double delz1 = xh1[2] - xo[2];

    double delx2 = xh2[0] - xo[0];
    double dely2 = xh2[1] - xo[1];
    double delz2 = xh2[2] - xo[2];

    xm[0] = xo[0] + this->alpha * 0.5 * (delx1 + delx2);
    xm[1] = xo[1] + this->alpha * 0.5 * (dely1 + dely2);
    xm[2] = xo[2] + this->alpha * 0.5 * (delz1 + delz2);

    this->domain->x2lamda(xm, xM);

  } else {

    iH1 = this->domain->closest_image(i, iH1);
    iH2 = this->domain->closest_image(i, iH2);

    double delx1 = x[iH1][0] - x[i][0];
    double dely1 = x[iH1][1] - x[i][1];
    double delz1 = x[iH1][2] - x[i][2];

    double delx2 = x[iH2][0] - x[i][0];
    double dely2 = x[iH2][1] - x[i][1];
    double delz2 = x[iH2][2] - x[i][2];

    xM[0] = x[i][0] + this->alpha * 0.5 * (delx1 + delx2);
    xM[1] = x[i][1] + this->alpha * 0.5 * (dely1 + dely2);
    xM[2] = x[i][2] + this->alpha * 0.5 * (delz1 + delz2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::tip4p_preprocess_host()
{
  alpha_kk = static_cast<KK_FLOAT>(this->alpha);

  const int nlocal = this->atom->nlocal;
  const int nmax_atom = this->atom->nmax;

  if (k_xM.extent(0) < (size_t) nmax_atom) {
    k_xM = DAT::tdual_kkfloat_1d_3("pppm/tip4p/kk:xM", nmax_atom);
    k_ih1 = DAT::tdual_int_1d("pppm/tip4p/kk:ih1", nmax_atom);
    k_ih2 = DAT::tdual_int_1d("pppm/tip4p/kk:ih2", nmax_atom);
  }

  double **x = this->atom->x;
  int *type = this->atom->type;
  auto h_xM = k_xM.view_host();
  auto h_ih1 = k_ih1.view_host();
  auto h_ih2 = k_ih2.view_host();

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == this->typeO) {
      int iH1, iH2;
      double xM[3];
      find_M_host(i, iH1, iH2, xM);
      h_xM(i, 0) = xM[0];
      h_xM(i, 1) = xM[1];
      h_xM(i, 2) = xM[2];
      h_ih1(i) = iH1;
      h_ih2(i) = iH2;
    } else {
      h_xM(i, 0) = x[i][0];
      h_xM(i, 1) = x[i][1];
      h_xM(i, 2) = x[i][2];
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

  d_xM = k_xM.template view<DeviceType>();
  d_ih1 = k_ih1.template view<DeviceType>();
  d_ih2 = k_ih2.template view<DeviceType>();
  d_type = this->atomKK->k_type.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::pp_pre_particle_map()
{
  tip4p_preprocess_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::particle_map()
{
  int nlocal = this->atomKK->nlocal;

  this->shift_kk = static_cast<KK_FLOAT>(this->shift);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPM_particle_map,
                                                                    const int &i) const
{
  const int nx =
      static_cast<int>((d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk + this->shift_kk) -
      OFFSET;
  const int ny =
      static_cast<int>((d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk + this->shift_kk) -
      OFFSET;
  const int nz =
      static_cast<int>((d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk + this->shift_kk) -
      OFFSET;

  this->d_part2grid(i, 0) = nx;
  this->d_part2grid(i, 1) = ny;
  this->d_part2grid(i, 2) = nz;

  if (nx + this->nlower < this->nxlo_out || nx + this->nupper > this->nxhi_out ||
      ny + this->nlower < this->nylo_out || ny + this->nupper > this->nyhi_out ||
      nz + this->nlower < this->nzlo_out || nz + this->nupper > this->nzhi_out)
    this->k_flag.template view<DeviceType>()() = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::make_rho()
{
  this->numz_out = this->nzhi_out - this->nzlo_out + 1;
  this->numy_out = this->nyhi_out - this->nylo_out + 1;
  this->numx_out = this->nxhi_out - this->nxlo_out + 1;
  const int inum_out = this->numz_out * this->numy_out * this->numx_out;

  this->shiftone_kk = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
  this->delvolinv_kk = static_cast<KK_FLOAT>(this->delvolinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_make_rho_zero>(0, inum_out),
                        static_cast<PPPMKokkos<DeviceType> &>(*this));
  this->copymode = 0;

  this->nlocal = this->atomKK->nlocal;

#ifdef LMP_KOKKOS_GPU
  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_make_rho_atomic>(0, this->nlocal),
                       *this);
  this->copymode = 0;
#else
  this->ix = this->nxhi_out - this->nxlo_out + 1;
  this->iy = this->nyhi_out - this->nylo_out + 1;

  this->copymode = 1;
  Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho> config(this->lmp->kokkos->nthreads, 1);
  Kokkos::parallel_for(config, *this);
  this->copymode = 0;
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(
    TagPPPMTIP4P_make_rho_atomic, const int &i) const
{
  Kokkos::View<FFT_SCALAR ***, Kokkos::LayoutRight, typename KKDevice<DeviceType>::value,
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
        const int mx = l + nx;
        a_density_brick(mz, my, mx) += x0 * this->d_rho1d(i, l + this->order / 2, 0);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(
    TagPPPMTIP4P_make_rho,
    typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type dev) const
{
  int tid = dev.league_rank();
  const int nthreads = dev.league_size();
  const int idelta = 1 + this->ngrid / nthreads;
  int ifrom = tid * idelta;
  int ito = ((ifrom + idelta) > this->ngrid) ? this->ngrid : ifrom + idelta;

  for (int i = 0; i < this->nlocal; i++) {

    int nx = this->d_part2grid(i, 0);
    int ny = this->d_part2grid(i, 1);
    int nz = this->d_part2grid(i, 2);

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_ik()
{
  int nlocal = this->atomKK->nlocal;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_fieldforce_ik>(0, nlocal), *this);
  this->copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPM_fieldforce_ik,
                                                                    const int &i) const
{
  int l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR x0, y0, z0;
  FFT_SCALAR ekx, eky, ekz;

  nx = this->d_part2grid(i, 0);
  ny = this->d_part2grid(i, 1);
  nz = this->d_part2grid(i, 2);

  nz -= this->nzlo_out;
  ny -= this->nylo_out;
  nx -= this->nxlo_out;

  ekx = eky = ekz = ZEROF;
  for (n = this->nlower; n <= this->nupper; n++) {
    mz = n + nz;
    z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (m = this->nlower; m <= this->nupper; m++) {
      my = m + ny;
      y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (l = this->nlower; l <= this->nupper; l++) {
        mx = l + nx;
        x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
        ekx -= x0 * this->d_vdx_brick(mz, my, mx);
        eky -= x0 * this->d_vdy_brick(mz, my, mx);
        ekz -= x0 * this->d_vdz_brick(mz, my, mx);
      }
    }
  }

  const KK_FLOAT qfactor = this->qscale_kk * this->q[i];
  const int itype = d_type(i);

  if (itype != this->typeO) {
    this->f(i, 0) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekx));
    this->f(i, 1) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(eky));
    if (this->slabflag != 2)
      this->f(i, 2) += static_cast<KK_ACC_FLOAT>(qfactor * static_cast<KK_FLOAT>(ekz));
  } else {
    const int iH1 = d_ih1(i);
    const int iH2 = d_ih2(i);
    const KK_FLOAT fx = qfactor * static_cast<KK_FLOAT>(ekx);
    const KK_FLOAT fy = qfactor * static_cast<KK_FLOAT>(eky);
    const KK_FLOAT fz = qfactor * static_cast<KK_FLOAT>(ekz);

    this->f(i, 0) += static_cast<KK_ACC_FLOAT>(fx * (1.0 - alpha_kk));
    this->f(i, 1) += static_cast<KK_ACC_FLOAT>(fy * (1.0 - alpha_kk));
    if (this->slabflag != 2)
      this->f(i, 2) += static_cast<KK_ACC_FLOAT>(fz * (1.0 - alpha_kk));

    Kokkos::atomic_add(&this->f(iH1, 0), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fx));
    Kokkos::atomic_add(&this->f(iH1, 1), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fy));
    if (this->slabflag != 2)
      Kokkos::atomic_add(&this->f(iH1, 2), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fz));

    Kokkos::atomic_add(&this->f(iH2, 0), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fx));
    Kokkos::atomic_add(&this->f(iH2, 1), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fy));
    if (this->slabflag != 2)
      Kokkos::atomic_add(&this->f(iH2, 2), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fz));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_peratom()
{
  int nlocal = this->atomKK->nlocal;

  this->shiftone_kk = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_fieldforce_peratom>(0, nlocal),
                       *this);
  this->copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPM_fieldforce_peratom,
                                                                    const int &i) const
{
  int l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;
  FFT_SCALAR u_pa, v0, v1, v2, v3, v4, v5;

  nx = this->d_part2grid(i, 0);
  ny = this->d_part2grid(i, 1);
  nz = this->d_part2grid(i, 2);
  dx = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nx) + this->shiftone_kk -
                               (d_xM(i, 0) - this->boxlo_kk[0]) * this->delxinv_kk);
  dy = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(ny) + this->shiftone_kk -
                               (d_xM(i, 1) - this->boxlo_kk[1]) * this->delyinv_kk);
  dz = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nz) + this->shiftone_kk -
                               (d_xM(i, 2) - this->boxlo_kk[2]) * this->delzinv_kk);

  nz -= this->nzlo_out;
  ny -= this->nylo_out;
  nx -= this->nxlo_out;

  this->compute_rho1d(i, dx, dy, dz);

  u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
  for (n = this->nlower; n <= this->nupper; n++) {
    mz = n + nz;
    z0 = this->d_rho1d(i, n + this->order / 2, 2);
    for (m = this->nlower; m <= this->nupper; m++) {
      my = m + ny;
      y0 = z0 * this->d_rho1d(i, m + this->order / 2, 1);
      for (l = this->nlower; l <= this->nupper; l++) {
        mx = l + nx;
        x0 = y0 * this->d_rho1d(i, l + this->order / 2, 0);
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

  const int itype = d_type(i);
  const int iH1 = d_ih1(i);
  const int iH2 = d_ih2(i);

  if (this->eflag_atom) {
    if (itype != this->typeO) {
      this->d_eatom[i] += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(u_pa);
    } else {
      this->d_eatom[i] +=
          static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(u_pa) *
          static_cast<KK_ACC_FLOAT>(1.0 - alpha_kk);
      this->d_eatom[iH1] +=
          static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(u_pa) *
          static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5);
      this->d_eatom[iH2] +=
          static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(u_pa) *
          static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5);
    }
  }
  if (this->vflag_atom) {
    if (itype != this->typeO) {
      this->d_vatom(i, 0) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v0);
      this->d_vatom(i, 1) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v1);
      this->d_vatom(i, 2) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v2);
      this->d_vatom(i, 3) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v3);
      this->d_vatom(i, 4) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v4);
      this->d_vatom(i, 5) += static_cast<KK_ACC_FLOAT>(this->q[i]) * static_cast<KK_ACC_FLOAT>(v5);
    } else {
      this->d_vatom(i, 0) +=
          static_cast<KK_ACC_FLOAT>(v0) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(i, 1) +=
          static_cast<KK_ACC_FLOAT>(v1) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(i, 2) +=
          static_cast<KK_ACC_FLOAT>(v2) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(i, 3) +=
          static_cast<KK_ACC_FLOAT>(v3) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(i, 4) +=
          static_cast<KK_ACC_FLOAT>(v4) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(i, 5) +=
          static_cast<KK_ACC_FLOAT>(v5) * static_cast<KK_ACC_FLOAT>((1.0 - alpha_kk) * this->q[i]);
      this->d_vatom(iH1, 0) +=
          static_cast<KK_ACC_FLOAT>(v0) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH1, 1) +=
          static_cast<KK_ACC_FLOAT>(v1) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH1, 2) +=
          static_cast<KK_ACC_FLOAT>(v2) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH1, 3) +=
          static_cast<KK_ACC_FLOAT>(v3) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH1, 4) +=
          static_cast<KK_ACC_FLOAT>(v4) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH1, 5) +=
          static_cast<KK_ACC_FLOAT>(v5) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 0) +=
          static_cast<KK_ACC_FLOAT>(v0) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 1) +=
          static_cast<KK_ACC_FLOAT>(v1) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 2) +=
          static_cast<KK_ACC_FLOAT>(v2) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 3) +=
          static_cast<KK_ACC_FLOAT>(v3) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 4) +=
          static_cast<KK_ACC_FLOAT>(v4) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
      this->d_vatom(iH2, 5) +=
          static_cast<KK_ACC_FLOAT>(v5) * static_cast<KK_ACC_FLOAT>(alpha_kk * 0.5 * this->q[i]);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::slabcorr()
{
  this->zprd_slab = this->domain->zprd * this->slab_volfactor;
  int nlocal = this->atomKK->nlocal;

  double dipole = 0.0;
  this->copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr1>(0, nlocal), *this,
                          dipole);
  this->copymode = 0;

  MPI_Allreduce(&dipole, &this->dipole_all, 1, MPI_DOUBLE, MPI_SUM, this->world);

  this->dipole_r2 = 0.0;
  if (this->eflag_atom || fabs(this->qsum) > SMALL) {
    this->copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr2>(0, nlocal), *this,
                            this->dipole_r2);
    this->copymode = 0;

    double tmp;
    MPI_Allreduce(&this->dipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, this->world);
    this->dipole_r2 = tmp;
  }

  const double e_slabcorr =
      MY_2PI * (this->dipole_all * this->dipole_all - this->qsum * this->dipole_r2 -
                this->qsum * this->qsum * this->zprd_slab * this->zprd_slab / 12.0) /
      this->volume;
  this->qscale = this->force->qqrd2e * this->scale;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  if (this->eflag_global) this->energy += this->qscale * e_slabcorr;

  if (this->eflag_atom) {
    this->efact = this->qscale * MY_2PI / this->volume;
    this->copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr3>(0, nlocal), *this);
    this->copymode = 0;
  }

  this->ffact = this->qscale * (-4.0 * MY_PI / this->volume);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr4>(0, nlocal), *this);
  this->copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr1,
                                                                    const int &i,
                                                                    double &dipole) const
{
  const int itype = d_type(i);
  if (itype == this->typeO)
    dipole += static_cast<double>(this->q[i] * d_xM(i, 2));
  else
    dipole += static_cast<double>(this->q[i] * this->x(i, 2));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr2,
                                                                    const int &i,
                                                                    double &dipole_r2) const
{
  dipole_r2 += static_cast<double>(this->q[i] * this->x(i, 2) * this->x(i, 2));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr3,
                                                                    const int &i) const
{
  double z_i = static_cast<double>(this->x(i, 2));
  double q_i = static_cast<double>(this->q[i]);
  this->d_eatom[i] +=
      static_cast<KK_ACC_FLOAT>(this->efact * q_i *
                                (z_i * this->dipole_all - 0.5 * (this->dipole_r2 + this->qsum * z_i * z_i) -
                                 this->qsum * this->zprd_slab * this->zprd_slab / 12.0));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr4,
                                                                    const int &i) const
{
  double z_i = static_cast<double>(this->x(i, 2));
  double q_i = static_cast<double>(this->q[i]);
  const double fzi_corr = this->ffact * q_i * (this->dipole_all - this->qsum * z_i);
  const int itype = d_type(i);
  if (itype == this->typeO) {
    const int iH1 = d_ih1(i);
    const int iH2 = d_ih2(i);
    this->f(i, 2) += static_cast<KK_ACC_FLOAT>(fzi_corr * (1.0 - alpha_kk));
    Kokkos::atomic_add(&this->f(iH1, 2), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fzi_corr));
    Kokkos::atomic_add(&this->f(iH2, 2), static_cast<KK_ACC_FLOAT>(0.5 * alpha_kk * fzi_corr));
  } else
    this->f(i, 2) += static_cast<KK_ACC_FLOAT>(fzi_corr);
}

namespace LAMMPS_NS {
template class PPPMTIP4PKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PPPMTIP4PKokkos<LMPHostType>;
#endif
}
