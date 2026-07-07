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

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int OFFSET = 16384;
static constexpr double SMALL = 0.00001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PPPMTIP4PKokkos<DeviceType>::PPPMTIP4PKokkos(LAMMPS *lmp) : Base(lmp)
{
  this->tip4pflag = 1;
  k_tip4p_flag = DAT::tdual_int_scalar("PPPM:tip4p_flag");
  nmax_tip4p = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::init()
{
  // TIP4P PPPM requires newton on, b/c it computes forces on ghost atoms

  if (this->force->newton == 0)
    this->error->all(FLERR,"Kspace style pppm/tip4p/kk requires newton on");

  if (this->domain->triclinic)
    this->error->all(FLERR,"Kspace style pppm/tip4p/kk does not (yet) support triclinic boxes");

  Base::init();
}

/* ----------------------------------------------------------------------
   bind TIP4P device data and (re)compute the off-atom M (charge) sites
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::compute_newsites()
{
  this->atomKK->sync(this->execution_space, TYPE_MASK | TAG_MASK);
  d_type = this->atomKK->k_type.template view<DeviceType>();
  d_tag = this->atomKK->k_tag.template view<DeviceType>();
  this->atomKK->k_sametag.template sync<DeviceType>();
  d_sametag = this->atomKK->k_sametag.template view<DeviceType>();

  map_style = this->atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = this->atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = this->atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  alpha_kk = static_cast<KK_FLOAT>(this->alpha);
  typeO_kk = this->typeO;
  typeH_kk = this->typeH;

  if (this->atom->nmax > nmax_tip4p) {
    nmax_tip4p = this->atom->nmax;
    d_xM = typename AT::t_kkfloat_1d_3("pppm/tip4p/kk:xM", nmax_tip4p);
    d_iH1 = typename AT::t_int_1d("pppm/tip4p/kk:iH1", nmax_tip4p);
    d_iH2 = typename AT::t_int_1d("pppm/tip4p/kk:iH2", nmax_tip4p);
  }

  k_tip4p_flag.view_host()() = 0;
  k_tip4p_flag.modify_host();
  k_tip4p_flag.template sync<DeviceType>();

  this->nlocal = this->atomKK->nlocal;
  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_findM>(0,this->nlocal),*this);
  this->copymode = 0;

  k_tip4p_flag.template modify<DeviceType>();
  k_tip4p_flag.sync_host();
  if (k_tip4p_flag.view_host()())
    this->error->one(FLERR,"TIP4P hydrogen is missing");
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_findM, const int &i) const
{
  if (d_type(i) == typeO_kk) {
    int iH1 = AtomKokkos::map_kokkos<DeviceType>(d_tag(i)+1,map_style,k_map_array,k_map_hash);
    int iH2 = AtomKokkos::map_kokkos<DeviceType>(d_tag(i)+2,map_style,k_map_array,k_map_hash);
    if (iH1 < 0 || iH2 < 0) {
      k_tip4p_flag.template view<DeviceType>()() = 1;
      d_iH1(i) = i; d_iH2(i) = i;
      d_xM(i,0) = x(i,0); d_xM(i,1) = x(i,1); d_xM(i,2) = x(i,2);
      return;
    }
    iH1 = closest_image(i,iH1);
    iH2 = closest_image(i,iH2);
    d_iH1(i) = iH1;
    d_iH2(i) = iH2;
    // xM = xO + alpha * 0.5 * ((xH1 - xO) + (xH2 - xO))
    d_xM(i,0) = x(i,0) + alpha_kk*(KK_FLOAT)0.5*(x(iH1,0) + x(iH2,0) - (KK_FLOAT)2.0*x(i,0));
    d_xM(i,1) = x(i,1) + alpha_kk*(KK_FLOAT)0.5*(x(iH1,1) + x(iH2,1) - (KK_FLOAT)2.0*x(i,1));
    d_xM(i,2) = x(i,2) + alpha_kk*(KK_FLOAT)0.5*(x(iH1,2) + x(iH2,2) - (KK_FLOAT)2.0*x(i,2));
  } else {
    d_xM(i,0) = x(i,0); d_xM(i,1) = x(i,1); d_xM(i,2) = x(i,2);
  }
}

/* ----------------------------------------------------------------------
   find grid points for all my particles, using the M site for O atoms
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::particle_map()
{
  this->shift_kk = static_cast<KK_FLOAT>(this->shift);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  if (!std::isfinite(this->boxlo[0]) || !std::isfinite(this->boxlo[1]) || !std::isfinite(this->boxlo[2]))
    this->error->one(FLERR,"Non-numeric box dimensions - simulation unstable" + utils::errorurl(6));

  // compute the off-atom M sites (also binds the TIP4P device data)

  compute_newsites();

  this->k_flag.view_host()() = 0;
  this->k_flag.modify_host();
  this->k_flag.template sync<DeviceType>();

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_particle_map>(0,this->nlocal),*this);
  this->copymode = 0;

  this->k_flag.template modify<DeviceType>();
  this->k_flag.sync_host();
  if (this->k_flag.view_host()())
    this->error->one(FLERR, Error::NOLASTLINE, "Out of range atoms - cannot compute PPPM" + utils::errorurl(4));
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_particle_map, const int &i) const
{
  const int nx = static_cast<int> ((d_xM(i,0)-boxlo_kk[0])*delxinv_kk+shift_kk) - OFFSET;
  const int ny = static_cast<int> ((d_xM(i,1)-boxlo_kk[1])*delyinv_kk+shift_kk) - OFFSET;
  const int nz = static_cast<int> ((d_xM(i,2)-boxlo_kk[2])*delzinv_kk+shift_kk) - OFFSET;

  d_part2grid(i,0) = nx;
  d_part2grid(i,1) = ny;
  d_part2grid(i,2) = nz;

  if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
      ny+nlower < nylo_out || ny+nupper > nyhi_out ||
      nz+nlower < nzlo_out || nz+nupper > nzhi_out)
    k_flag.template view<DeviceType>()() = 1;
}

/* ----------------------------------------------------------------------
   create discretized "density" on grid; O charge sits on its M site
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::make_rho()
{
  // clear 3d density array

  numz_out = nzhi_out-nzlo_out + 1;
  numy_out = nyhi_out-nylo_out + 1;
  numx_out = nxhi_out-nxlo_out + 1;
  const int inum_out = numz_out*numy_out*numx_out;

  this->shiftone_kk = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
  this->delvolinv_kk = static_cast<KK_FLOAT>(this->delvolinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_make_rho_zero>(0,inum_out),*this);
  this->copymode = 0;

  this->nlocal = this->atomKK->nlocal;

#ifdef LMP_KOKKOS_GPU
  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_make_rho_atomic>(0,this->nlocal),*this);
  this->copymode = 0;
#else
  ix = nxhi_out-nxlo_out + 1;
  iy = nyhi_out-nylo_out + 1;

  this->copymode = 1;
  Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho> config(this->lmp->kokkos->nthreads,1);
  Kokkos::parallel_for(config,*this);
  this->copymode = 0;
#endif
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_make_rho_atomic, const int &i) const
{
  // The density_brick array is atomic for Half/Thread neighbor style
  Kokkos::View<FFT_SCALAR***,Kokkos::LayoutRight,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_density_brick = d_density_brick;

  int nx = d_part2grid(i,0);
  int ny = d_part2grid(i,1);
  int nz = d_part2grid(i,2);
  const FFT_SCALAR dx = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nx)+shiftone_kk - (d_xM(i,0)-boxlo_kk[0])*delxinv_kk);
  const FFT_SCALAR dy = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(ny)+shiftone_kk - (d_xM(i,1)-boxlo_kk[1])*delyinv_kk);
  const FFT_SCALAR dz = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nz)+shiftone_kk - (d_xM(i,2)-boxlo_kk[2])*delzinv_kk);

  nz -= nzlo_out;
  ny -= nylo_out;
  nx -= nxlo_out;

  compute_rho1d(i,dx,dy,dz);

  const FFT_SCALAR z0 = static_cast<FFT_SCALAR>(delvolinv_kk * q[i]);
  for (int n = nlower; n <= nupper; n++) {
    const int mz = n+nz;
    const FFT_SCALAR y0 = z0*d_rho1d(i,n+order/2,2);
    for (int m = nlower; m <= nupper; m++) {
      const int my = m+ny;
      const FFT_SCALAR x0 = y0*d_rho1d(i,m+order/2,1);
      for (int l = nlower; l <= nupper; l++) {
        const int mx = l+nx;
        a_density_brick(mz,my,mx) += x0*d_rho1d(i,l+order/2,0);
      }
    }
  }
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator() (TagPPPMTIP4P_make_rho, typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type dev) const {
  // determine range of grid points handled by this thread
  int tid = dev.league_rank();
  const int nthreads = dev.league_size();
  const int idelta = 1 + ngrid/nthreads;
  int ifrom = tid*idelta;
  int ito = ((ifrom + idelta) > ngrid) ? ngrid : ifrom + idelta;

  for (int i = 0; i < nlocal; i++) {

    int nx = d_part2grid(i,0);
    int ny = d_part2grid(i,1);
    int nz = d_part2grid(i,2);

    // pre-screen whether this atom will ever come within
    // reach of the data segement this thread is updating.
    if ( ((nz+nlower-nzlo_out)*ix*iy >= ito)
         || ((nz+nupper-nzlo_out+1)*ix*iy < ifrom) ) continue;

    const FFT_SCALAR dx = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nx)+shiftone_kk - (d_xM(i,0)-boxlo_kk[0])*delxinv_kk);
    const FFT_SCALAR dy = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(ny)+shiftone_kk - (d_xM(i,1)-boxlo_kk[1])*delyinv_kk);
    const FFT_SCALAR dz = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nz)+shiftone_kk - (d_xM(i,2)-boxlo_kk[2])*delzinv_kk);

    nz -= nzlo_out;
    ny -= nylo_out;
    nx -= nxlo_out;

    compute_rho1d(i,dx,dy,dz);

    const FFT_SCALAR z0 = static_cast<FFT_SCALAR>(delvolinv_kk * q[i]);
    for (int n = nlower; n <= nupper; n++) {
      const int mz = n+nz;
      const int in = mz*ix*iy;
      const FFT_SCALAR y0 = z0*d_rho1d(i,n+order/2,2);
      for (int m = nlower; m <= nupper; m++) {
        const int my = m+ny;
        const int im = in+my*ix;
        const FFT_SCALAR x0 = y0*d_rho1d(i,m+order/2,1);
        for (int l = nlower; l <= nupper; l++) {
          const int mx = l+nx;
          const int il = im+mx;
          if (il >= ito) break;
          if (il < ifrom) continue;
          d_density_brick(mz,my,mx) += x0*d_rho1d(i,l+order/2,0);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate electric field and force; redistribute O force onto O + 2H
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_ik()
{
  this->nlocal = this->atomKK->nlocal;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_fieldforce_ik>(0,this->nlocal),*this);
  this->copymode = 0;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_fieldforce_ik, const int &i) const
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  nx = d_part2grid(i,0);
  ny = d_part2grid(i,1);
  nz = d_part2grid(i,2);

  nz -= nzlo_out;
  ny -= nylo_out;
  nx -= nxlo_out;

  ekx = eky = ekz = (FFT_SCALAR)0.0;
  for (n = nlower; n <= nupper; n++) {
    mz = n+nz;
    z0 = d_rho1d(i,n+order/2,2);
    for (m = nlower; m <= nupper; m++) {
      my = m+ny;
      y0 = z0*d_rho1d(i,m+order/2,1);
      for (l = nlower; l <= nupper; l++) {
        mx = l+nx;
        x0 = y0*d_rho1d(i,l+order/2,0);
        ekx -= x0*d_vdx_brick(mz,my,mx);
        eky -= x0*d_vdy_brick(mz,my,mx);
        ekz -= x0*d_vdz_brick(mz,my,mx);
      }
    }
  }

  // convert E-field to force, redistribute force on O onto O + 2 H atoms

  const KK_FLOAT qfactor = qscale_kk * q[i];
  if (d_type(i) != typeO_kk) {
    Kokkos::atomic_add(&f(i,0), static_cast<KK_ACC_FLOAT>(qfactor*static_cast<KK_FLOAT>(ekx)));
    Kokkos::atomic_add(&f(i,1), static_cast<KK_ACC_FLOAT>(qfactor*static_cast<KK_FLOAT>(eky)));
    if (slabflag != 2)
      Kokkos::atomic_add(&f(i,2), static_cast<KK_ACC_FLOAT>(qfactor*static_cast<KK_FLOAT>(ekz)));
  } else {
    const KK_FLOAT fx = qfactor*static_cast<KK_FLOAT>(ekx);
    const KK_FLOAT fy = qfactor*static_cast<KK_FLOAT>(eky);
    const KK_FLOAT fz = qfactor*static_cast<KK_FLOAT>(ekz);
    const int iH1 = d_iH1(i);
    const int iH2 = d_iH2(i);
    Kokkos::atomic_add(&f(i,0), static_cast<KK_ACC_FLOAT>(fx*((KK_FLOAT)1.0-alpha_kk)));
    Kokkos::atomic_add(&f(i,1), static_cast<KK_ACC_FLOAT>(fy*((KK_FLOAT)1.0-alpha_kk)));
    Kokkos::atomic_add(&f(iH1,0), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fx));
    Kokkos::atomic_add(&f(iH1,1), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fy));
    Kokkos::atomic_add(&f(iH2,0), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fx));
    Kokkos::atomic_add(&f(iH2,1), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fy));
    if (slabflag != 2) {
      Kokkos::atomic_add(&f(i,2), static_cast<KK_ACC_FLOAT>(fz*((KK_FLOAT)1.0-alpha_kk)));
      Kokkos::atomic_add(&f(iH1,2), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fz));
      Kokkos::atomic_add(&f(iH2,2), static_cast<KK_ACC_FLOAT>((KK_FLOAT)0.5*alpha_kk*fz));
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate per-atom energy/virial; redistribute O contribution onto O + 2H
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::fieldforce_peratom()
{
  this->nlocal = this->atomKK->nlocal;

  this->shiftone_kk = static_cast<KK_FLOAT>(this->shiftone);
  this->delxinv_kk = static_cast<KK_FLOAT>(this->delxinv);
  this->delyinv_kk = static_cast<KK_FLOAT>(this->delyinv);
  this->delzinv_kk = static_cast<KK_FLOAT>(this->delzinv);
  this->boxlo_kk[0] = static_cast<KK_FLOAT>(this->boxlo[0]);
  this->boxlo_kk[1] = static_cast<KK_FLOAT>(this->boxlo[1]);
  this->boxlo_kk[2] = static_cast<KK_FLOAT>(this->boxlo[2]);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_fieldforce_peratom>(0,this->nlocal),*this);
  this->copymode = 0;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_fieldforce_peratom, const int &i) const
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

  nx = d_part2grid(i,0);
  ny = d_part2grid(i,1);
  nz = d_part2grid(i,2);
  dx = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nx)+shiftone_kk - (d_xM(i,0)-boxlo_kk[0])*delxinv_kk);
  dy = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(ny)+shiftone_kk - (d_xM(i,1)-boxlo_kk[1])*delyinv_kk);
  dz = static_cast<FFT_SCALAR>(static_cast<KK_FLOAT>(nz)+shiftone_kk - (d_xM(i,2)-boxlo_kk[2])*delzinv_kk);

  nz -= nzlo_out;
  ny -= nylo_out;
  nx -= nxlo_out;

  compute_rho1d(i,dx,dy,dz);

  u = v0 = v1 = v2 = v3 = v4 = v5 = (FFT_SCALAR)0.0;
  for (n = nlower; n <= nupper; n++) {
    mz = n+nz;
    z0 = d_rho1d(i,n+order/2,2);
    for (m = nlower; m <= nupper; m++) {
      my = m+ny;
      y0 = z0*d_rho1d(i,m+order/2,1);
      for (l = nlower; l <= nupper; l++) {
        mx = l+nx;
        x0 = y0*d_rho1d(i,l+order/2,0);
        if (eflag_atom) u += x0*d_u_brick(mz,my,mx);
        if (vflag_atom) {
          v0 += x0*d_v0_brick(mz,my,mx);
          v1 += x0*d_v1_brick(mz,my,mx);
          v2 += x0*d_v2_brick(mz,my,mx);
          v3 += x0*d_v3_brick(mz,my,mx);
          v4 += x0*d_v4_brick(mz,my,mx);
          v5 += x0*d_v5_brick(mz,my,mx);
        }
      }
    }
  }

  const KK_FLOAT qi = q[i];
  if (d_type(i) != typeO_kk) {
    if (eflag_atom) Kokkos::atomic_add(&d_eatom[i], static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(u)));
    if (vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v0)));
      Kokkos::atomic_add(&d_vatom(i,1), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v1)));
      Kokkos::atomic_add(&d_vatom(i,2), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v2)));
      Kokkos::atomic_add(&d_vatom(i,3), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v3)));
      Kokkos::atomic_add(&d_vatom(i,4), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v4)));
      Kokkos::atomic_add(&d_vatom(i,5), static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(v5)));
    }
  } else {
    const int iH1 = d_iH1(i);
    const int iH2 = d_iH2(i);
    const KK_FLOAT a = alpha_kk;
    if (eflag_atom) {
      Kokkos::atomic_add(&d_eatom[i],   static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(u)*((KK_FLOAT)1.0-a)));
      Kokkos::atomic_add(&d_eatom[iH1], static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(u)*a*(KK_FLOAT)0.5));
      Kokkos::atomic_add(&d_eatom[iH2], static_cast<KK_ACC_FLOAT>(qi*static_cast<KK_FLOAT>(u)*a*(KK_FLOAT)0.5));
    }
    if (vflag_atom) {
      const FFT_SCALAR vv[6] = {v0,v1,v2,v3,v4,v5};
      for (int k = 0; k < 6; k++) {
        const KK_FLOAT vk = qi*static_cast<KK_FLOAT>(vv[k]);
        Kokkos::atomic_add(&d_vatom(i,k),   static_cast<KK_ACC_FLOAT>(vk*((KK_FLOAT)1.0-a)));
        Kokkos::atomic_add(&d_vatom(iH1,k), static_cast<KK_ACC_FLOAT>(vk*a*(KK_FLOAT)0.5));
        Kokkos::atomic_add(&d_vatom(iH2,k), static_cast<KK_ACC_FLOAT>(vk*a*(KK_FLOAT)0.5));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   slab correction with all terms evaluated at the M site;
   per-atom energy and force redistributed onto O + 2H
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMTIP4PKokkos<DeviceType>::slabcorr()
{
  this->zprd_slab = this->domain->zprd*this->slab_volfactor;
  this->nlocal = this->atomKK->nlocal;

  // M-site dipole moment (uses xM for O atoms)

  double dipole = 0.0;
  this->copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr1>(0,this->nlocal),*this,dipole);
  this->copymode = 0;

  MPI_Allreduce(&dipole,&this->dipole_all,1,MPI_DOUBLE,MPI_SUM,this->world);

  // dipole_r2 (also evaluated at the M site for O atoms)

  this->dipole_r2 = 0.0;
  if (this->eflag_atom || fabs(this->qsum) > SMALL) {
    double dipole_r2 = 0.0;
    this->copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr2>(0,this->nlocal),*this,dipole_r2);
    this->copymode = 0;

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,this->world);
    this->dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(this->dipole_all*this->dipole_all -
    this->qsum*this->dipole_r2 - this->qsum*this->qsum*this->zprd_slab*this->zprd_slab/12.0)/this->volume;
  this->qscale = this->qqrd2e * this->scale;
  this->qscale_kk = static_cast<KK_FLOAT>(this->qscale);

  if (this->eflag_global) this->energy += this->qscale * e_slabcorr;

  // per-atom energy, evaluated at the M site and redistributed onto O + 2H

  if (this->eflag_atom) {
    this->efact = this->qscale * MY_2PI/this->volume;
    this->copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr3>(0,this->nlocal),*this);
    this->copymode = 0;
  }

  // add on force corrections, redistributing O force onto O + 2H

  this->ffact = this->qscale * (-4.0*MY_PI/this->volume);

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPMTIP4P_slabcorr4>(0,this->nlocal),*this);
  this->copymode = 0;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr1, const int &i, double &dipole) const
{
  dipole += static_cast<double>(q[i])*static_cast<double>(d_xM(i,2));
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr2, const int &i, double &dipole_r2) const
{
  const double zM = static_cast<double>(d_xM(i,2));
  dipole_r2 += static_cast<double>(q[i])*zM*zM;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr3, const int &i) const
{
  const double zM = static_cast<double>(d_xM(i,2));
  const double e_pa = efact * static_cast<double>(q[i])*(zM*dipole_all -
    0.5*(dipole_r2 + qsum*zM*zM) - qsum*zprd_slab*zprd_slab/12.0);
  if (d_type(i) != typeO_kk) {
    Kokkos::atomic_add(&d_eatom[i], static_cast<KK_ACC_FLOAT>(e_pa));
  } else {
    const double a = static_cast<double>(alpha_kk);
    Kokkos::atomic_add(&d_eatom[i],        static_cast<KK_ACC_FLOAT>(e_pa*(1.0-a)));
    Kokkos::atomic_add(&d_eatom[d_iH1(i)], static_cast<KK_ACC_FLOAT>(0.5*a*e_pa));
    Kokkos::atomic_add(&d_eatom[d_iH2(i)], static_cast<KK_ACC_FLOAT>(0.5*a*e_pa));
  }
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PPPMTIP4PKokkos<DeviceType>::operator()(TagPPPMTIP4P_slabcorr4, const int &i) const
{
  // d_xM holds the M site for O atoms and the atom position otherwise
  const double z_i = static_cast<double>(d_xM(i,2));
  const double q_i = static_cast<double>(q[i]);
  const double fzi = ffact * q_i*(dipole_all - qsum*z_i);
  if (d_type(i) != typeO_kk) {
    Kokkos::atomic_add(&f(i,2), static_cast<KK_ACC_FLOAT>(fzi));
  } else {
    const int iH1 = d_iH1(i);
    const int iH2 = d_iH2(i);
    const double a = static_cast<double>(alpha_kk);
    Kokkos::atomic_add(&f(i,2),   static_cast<KK_ACC_FLOAT>(fzi*(1.0-a)));
    Kokkos::atomic_add(&f(iH1,2), static_cast<KK_ACC_FLOAT>(0.5*a*fzi));
    Kokkos::atomic_add(&f(iH2,2), static_cast<KK_ACC_FLOAT>(0.5*a*fzi));
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PPPMTIP4PKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PPPMTIP4PKokkos<LMPHostType>;
#endif
}
