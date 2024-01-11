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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pppm_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fft3d_kokkos.h"
#include "force.h"
#include "grid3d_kokkos.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "pair.h"
#include "remap_kokkos.h"
#include "kokkos_few.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecialKokkos;

#define MAXORDER 7
#define OFFSET 16384
#define LARGE 10000.0
#define SMALL 0.00001
#define EPS_HOC 1.0e-7

enum{REVERSE_RHO};
enum{FORWARD_IK,FORWARD_IK_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PPPMKokkos<DeviceType>::PPPMKokkos(LAMMPS *lmp) : PPPM(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK;
  datamask_modify = F_MASK;

  group_group_enable = 0;
  triclinic_support = 1;

  peratom_allocate_flag = 0;

  // define acons coefficients for estimation of kspace errors
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

  acons = typename Kokkos::DualView<F_FLOAT[8][7],Kokkos::LayoutRight,DeviceType>::t_host("pppm:acons");
  acons(1,0) = 2.0 / 3.0;
  acons(2,0) = 1.0 / 50.0;
  acons(2,1) = 5.0 / 294.0;
  acons(3,0) = 1.0 / 588.0;
  acons(3,1) = 7.0 / 1440.0;
  acons(3,2) = 21.0 / 3872.0;
  acons(4,0) = 1.0 / 4320.0;
  acons(4,1) = 3.0 / 1936.0;
  acons(4,2) = 7601.0 / 2271360.0;
  acons(4,3) = 143.0 / 28800.0;
  acons(5,0) = 1.0 / 23232.0;
  acons(5,1) = 7601.0 / 13628160.0;
  acons(5,2) = 143.0 / 69120.0;
  acons(5,3) = 517231.0 / 106536960.0;
  acons(5,4) = 106640677.0 / 11737571328.0;
  acons(6,0) = 691.0 / 68140800.0;
  acons(6,1) = 13.0 / 57600.0;
  acons(6,2) = 47021.0 / 35512320.0;
  acons(6,3) = 9694607.0 / 2095994880.0;
  acons(6,4) = 733191589.0 / 59609088000.0;
  acons(6,5) = 326190917.0 / 11700633600.0;
  acons(7,0) = 1.0 / 345600.0;
  acons(7,1) = 3617.0 / 35512320.0;
  acons(7,2) = 745739.0 / 838397952.0;
  acons(7,3) = 56399353.0 / 12773376000.0;
  acons(7,4) = 25091609.0 / 1560084480.0;
  acons(7,5) = 1755948832039.0 / 36229939200000.0;
  acons(7,6) = 4887769399.0 / 37838389248.0;

  k_flag = DAT::tdual_int_scalar("PPPM:flag");

  // same name but different than base class

  gc = nullptr;
  fft1 = nullptr;
  fft2 = nullptr;
  remap = nullptr;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

template<class DeviceType>
PPPMKokkos<DeviceType>::~PPPMKokkos()
{
  if (copymode) return;

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::init()
{
  if (me == 0) utils::logmesg(lmp,"PPPM Kokkos initialization ...\n");

  // error check

  if (differentiation_flag == 1)
    error->all(FLERR,"Cannot (yet) use PPPM Kokkos with 'kspace_modify diff ad'");

  triclinic_check();

  if (triclinic != domain->triclinic)
    error->all(FLERR,"Must redefine kspace_style after changing to triclinic box");

  if (domain->triclinic && slabflag)
    error->all(FLERR,"Cannot (yet) use PPPM with triclinic box and slab correction");
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use PPPM with 2d simulation");

  if (!atom->q_flag)
    error->all(FLERR,"Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with PPPM");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPM");
  }

  if (order < 2 || order > MAXORDER)
    error->all(FLERR,"PPPM order cannot be < 2 or > {}",MAXORDER);

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  pair_check();

  int itmp = 0;
  auto p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == nullptr)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // if kspace is TIP4P, extract TIP4P params from pair style
  // bond/angle are not yet init(), so ensure equilibrium request is valid

  qdist = 0.0;

  if (tip4pflag)
      error->all(FLERR,"Cannot (yet) use PPPM Kokkos TIP4P");

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  int iteration = 0;

  while (order >= minorder) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing PPPM order b/c stencil extends "
                     "beyond nearest neighbor processor");

    set_grid_global();
    set_grid_local();
    if (overlap_allowed) break;

    gc = new Grid3dKokkos<DeviceType>(lmp,world,nx_pppm,ny_pppm,nz_pppm);
    gc->set_distance(0.5*neighbor->skin + qdist);
    gc->set_stencil_atom(-nlower,nupper);
    gc->set_shift_atom(shiftatom_lo,shiftatom_hi);
    gc->set_zfactor(slab_volfactor);

    gc->setup_grid(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

    int tmp1,tmp2;
    gc->setup_comm(tmp1,tmp2);
    if (gc->ghost_adjacent()) break;
    delete gc;

    order--;
    iteration++;
  }

  if (order < minorder) error->all(FLERR,"PPPM order < minimum allowed order");
  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR,"PPPM grid stencil extends beyond nearest neighbor processor");
  if (gc) delete gc;

  // adjust g_ewald

  if (!gewaldflag) adjust_gewald();

  // calculate the final accuracy

  double estimated_accuracy = final_accuracy();

  // allocate K-space dependent memory
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  compute_rho_coeff();

  k_rho_coeff.template modify<LMPHostType>();
  k_rho_coeff.template sync<DeviceType>();

  // print stats

  int ngrid_max,nfft_both_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    std::string mesg = fmt::format("  G vector (1/distance) = {:.8g}\n",g_ewald);
    mesg += fmt::format("  grid = {} {} {}\n",nx_pppm,ny_pppm,nz_pppm);
    mesg += fmt::format("  stencil order = {}\n",order);
    mesg += fmt::format("  estimated absolute RMS force accuracy = {:.8g}\n",
                       estimated_accuracy);
    mesg += fmt::format("  estimated relative force accuracy = {:.8g}\n",
                       estimated_accuracy/two_charge_force);
    mesg += "  using " LMP_FFT_PREC " precision " LMP_FFT_LIB "\n";
    mesg += fmt::format("  3d grid and FFT values/proc = {} {}\n",
                       ngrid_max,nfft_both_max);
    utils::logmesg(lmp,mesg);
  }
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::setup()
{
  if (triclinic) {
    setup_triclinic();
    return;
  }

  // perform some checks to avoid illegal boundaries with read_data

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with PPPM");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPM");
  }

  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  unitkx = (MY_2PI/xprd);
  unitky = (MY_2PI/yprd);
  unitkz = (MY_2PI/zprd_slab);

  // d_fkx,d_fky,d_fkz for my FFT grid pts

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup1>(nxlo_fft,nxhi_fft+1),*this);
  copymode = 0;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup2>(nylo_fft,nyhi_fft+1),*this);
  copymode = 0;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup3>(nzlo_fft,nzhi_fft+1),*this);
  copymode = 0;

  // merge three outer loops into one for better threading

  numz_fft = nzhi_fft-nzlo_fft + 1;
  numy_fft = nyhi_fft-nylo_fft + 1;
  numx_fft = nxhi_fft-nxlo_fft + 1;
  const int inum_fft = numz_fft*numy_fft*numx_fft;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup4>(0,inum_fft),*this);
  copymode = 0;

  compute_gf_ik();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup1, const int &i) const
{
  double per = i - nx_pppm*(2*i/nx_pppm);
  d_fkx[i-nxlo_fft] = unitkx*per;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup2, const int &i) const
{
  double per = i - ny_pppm*(2*i/ny_pppm);
  d_fky[i-nylo_fft] = unitky*per;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup3, const int &i) const
{
  double per = i - nz_pppm*(2*i/nz_pppm);
  d_fkz[i-nzlo_fft] = unitkz*per;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup4, const int &n) const
{
  const int k = n/(numy_fft*numx_fft);
  const int j = (n - k*numy_fft*numx_fft) / numx_fft;
  const int i = n - k*numy_fft*numx_fft - j*numx_fft;
  const double sqk = d_fkx[i]*d_fkx[i] + d_fky[j]*d_fky[j] + d_fkz[k]*d_fkz[k];
  if (sqk == 0.0) {
    d_vg(n,0) = 0.0;
    d_vg(n,1) = 0.0;
    d_vg(n,2) = 0.0;
    d_vg(n,3) = 0.0;
    d_vg(n,4) = 0.0;
    d_vg(n,5) = 0.0;
  } else {
    const double vterm = -2.0 * (1.0/sqk + 0.25/(g_ewald*g_ewald));
    d_vg(n,0) = 1.0 + vterm*d_fkx[i]*d_fkx[i];
    d_vg(n,1) = 1.0 + vterm*d_fky[j]*d_fky[j];
    d_vg(n,2) = 1.0 + vterm*d_fkz[k]*d_fkz[k];
    d_vg(n,3) = vterm*d_fkx[i]*d_fky[j];
    d_vg(n,4) = vterm*d_fkx[i]*d_fkz[k];
    d_vg(n,5) = vterm*d_fky[j]*d_fkz[k];
  }
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
   for a triclinic system
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::setup_triclinic()
{
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  prd = domain->prd;
  // Update simulation box parameters
  h = Few<double, 6>(domain->h);
  h_inv = Few<double, 6>(domain->h_inv);

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  // use lamda (0-1) coordinates

  delxinv = nx_pppm;
  delyinv = ny_pppm;
  delzinv = nz_pppm;
  delvolinv = delxinv*delyinv*delzinv/volume;

  numz_fft = nzhi_fft-nzlo_fft + 1;
  numy_fft = nyhi_fft-nylo_fft + 1;
  numx_fft = nxhi_fft-nxlo_fft + 1;
  const int inum_fft = numz_fft*numy_fft*numx_fft;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup_triclinic1>(0,inum_fft),*this);
  copymode = 0;

  // virial coefficients

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_setup_triclinic2>(0,nfft),*this);
  copymode = 0;

  compute_gf_ik_triclinic();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup_triclinic1, const int &n) const
{
  int k = n/(numy_fft*numx_fft);
  int j = (n - k*numy_fft*numx_fft) / numx_fft;
  int i = n - k*numy_fft*numx_fft - j*numx_fft;
  k += nzlo_fft;
  j += nylo_fft;
  i += nxlo_fft;

  double per_k = k - nz_pppm*(2*k/nz_pppm);
  double per_j = j - ny_pppm*(2*j/ny_pppm);
  double per_i = i - nx_pppm*(2*i/nx_pppm);

  double unitk_lamda[3];
  unitk_lamda[0] = 2.0*MY_PI*per_i;
  unitk_lamda[1] = 2.0*MY_PI*per_j;
  unitk_lamda[2] = 2.0*MY_PI*per_k;
  x2lamdaT_kokkos(&unitk_lamda[0],&unitk_lamda[0]);
  d_fkx[n] = unitk_lamda[0];
  d_fky[n] = unitk_lamda[1];
  d_fkz[n] = unitk_lamda[2];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_setup_triclinic2, const int &n) const
{
  const double sqk = d_fkx[n]*d_fkx[n] + d_fky[n]*d_fky[n] + d_fkz[n]*d_fkz[n];
  if (sqk == 0.0) {
    d_vg(n,0) = 0.0;
    d_vg(n,1) = 0.0;
    d_vg(n,2) = 0.0;
    d_vg(n,3) = 0.0;
    d_vg(n,4) = 0.0;
    d_vg(n,5) = 0.0;
  } else {
    const double vterm = -2.0 * (1.0/sqk + 0.25/(g_ewald*g_ewald));
    d_vg(n,0) = 1.0 + vterm*d_fkx[n]*d_fkx[n];
    d_vg(n,1) = 1.0 + vterm*d_fky[n]*d_fky[n];
    d_vg(n,2) = 1.0 + vterm*d_fkz[n]*d_fkz[n];
    d_vg(n,3) = vterm*d_fkx[n]*d_fky[n];
    d_vg(n,4) = vterm*d_fkx[n]*d_fkz[n];
    d_vg(n,5) = vterm*d_fky[n]*d_fkz[n];
  }
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::compute(int eflag, int vflag)
{
  int i;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  if (evflag_atom && !peratom_allocate_flag)
    allocate_peratom();

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();

  //nlocal = atomKK->nlocal;
  //nall = atomKK->nlocal + atomKK->nghost;
  //newton_pair = force->newton_pair;

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) {
    boxlo[0] = domain->boxlo[0];
    boxlo[1] = domain->boxlo[1];
    boxlo[2] = domain->boxlo[2];
  } else {
    boxlo[0] = domain->boxlo_lamda[0];
    boxlo[1] = domain->boxlo_lamda[1];
    boxlo[2] = domain->boxlo_lamda[2];
    domain->x2lamda(atomKK->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    nmax = atomKK->nmax;
    d_part2grid = typename AT::t_int_1d_3("pppm:part2grid",nmax);
    d_rho1d = typename FFT_AT::t_FFT_SCALAR_2d_3("pppm:rho1d",nmax,order/2+order/2+1);
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   k_gc_buf1,k_gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                   k_gc_buf1,k_gc_buf2,MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom)
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                     k_gc_buf1,k_gc_buf2,MPI_FFT_SCALAR);

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // ntotal accounts for TIP4P tallying eatom/vatom for ghost atoms

  if (evflag_atom) {
    int nlocal = atomKK->nlocal;
    int ntotal = nlocal;

    if (eflag_atom) {
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_self1>(0,nlocal),*this);
      copymode = 0;
    }

    if (vflag_atom) {
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_self2>(0,ntotal),*this);
      copymode = 0;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_self1, const int &i) const
{
  d_eatom[i] *= 0.5;
  d_eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
    (g_ewald*g_ewald*volume);
  d_eatom[i] *= qscale;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_self2, const int &i) const
{
  for (int j = 0; j < 6; j++) d_vatom(i,j) *= 0.5*qscale;
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::allocate()
{
  // create ghost grid object for rho and electric field communication
  // returns local owned and ghost grid bounds
  // setup communication patterns and buffers

  gc = new Grid3dKokkos<DeviceType>(lmp,world,nx_pppm,ny_pppm,nz_pppm);
  gc->set_distance(0.5*neighbor->skin + qdist);
  gc->set_stencil_atom(-nlower,nupper);
  gc->set_shift_atom(shiftatom_lo,shiftatom_hi);
  gc->set_zfactor(slab_volfactor);

  gc->setup_grid(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                 nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

  gc->setup_comm(ngc_buf1,ngc_buf2);

  npergrid = 3;

  k_gc_buf1 = FFT_DAT::tdual_FFT_SCALAR_1d("pppm:gc_buf1",npergrid*ngc_buf1);
  k_gc_buf2 = FFT_DAT::tdual_FFT_SCALAR_1d("pppm:gc_buf2",npergrid*ngc_buf2);

  // tally local grid sizes
  // ngrid = count of owned+ghost grid cells on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  //              same as count of owned grid cells
  // nfft = FFT points in x-pencil FFT decomposition on this proc
  // nfft_both = greater of nfft and nfft_brick

  ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);

  nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
    (nzhi_in-nzlo_in+1);

  nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);

  nfft_both = MAX(nfft,nfft_brick);

  // allocate distributed grid data

  d_density_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:density_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);

  memoryKK->create_kokkos(k_density_fft,density_fft,nfft_both,"pppm:d_density_fft");
  d_density_fft = k_density_fft.view<DeviceType>();

  d_greensfn = typename AT::t_float_1d("pppm:greensfn",nfft_both);
  memoryKK->create_kokkos(k_work1,work1,2*nfft_both,"pppm:work1");
  memoryKK->create_kokkos(k_work2,work2,2*nfft_both,"pppm:work2");
  d_work1 = k_work1.view<DeviceType>();
  d_work2 = k_work2.view<DeviceType>();
  d_vg = typename AT::t_virial_array("pppm:vg",nfft_both);

  if (triclinic == 0) {
    d_fkx = typename AT::t_float_1d("pppm:d_fkx",nxhi_fft-nxlo_fft+1);
    d_fky = typename AT::t_float_1d("pppm:d_fky",nyhi_fft-nylo_fft+1);
    d_fkz = typename AT::t_float_1d("pppm:d_fkz",nzhi_fft-nzlo_fft+1);
  } else {
    d_fkx = typename AT::t_float_1d("pppm:d_fkx",nfft_both);
    d_fky = typename AT::t_float_1d("pppm:d_fky",nfft_both);
    d_fkz = typename AT::t_float_1d("pppm:d_fkz",nfft_both);
  }

  d_vdx_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_vdx_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_vdy_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_vdy_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_vdz_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_vdz_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);

  // summation coeffs

  order_allocated = order;
  k_gf_b = typename DAT::tdual_float_1d("pppm:gf_b",order);
  d_gf_b = k_gf_b.view<DeviceType>();
  d_rho1d = typename FFT_AT::t_FFT_SCALAR_2d_3("pppm:rho1d",nmax,order/2+order/2+1);
  k_rho_coeff = FFT_DAT::tdual_FFT_SCALAR_2d("pppm:rho_coeff",order,order/2-(1-order)/2+1);
  d_rho_coeff = k_rho_coeff.view<DeviceType>();
  h_rho_coeff = k_rho_coeff.h_view;

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decomposition
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int collective_flag = 0; // not yet supported in Kokkos version
  int gpu_aware_flag = lmp->kokkos->gpu_aware_flag;
  int tmp;

  fft1 = new FFT3dKokkos<DeviceType>(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                         nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                         nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                         0,0,&tmp,collective_flag,gpu_aware_flag);

  fft2 = new FFT3dKokkos<DeviceType>(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                         nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                         nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                         0,0,&tmp,collective_flag,gpu_aware_flag);
  remap = new RemapKokkos<DeviceType>(lmp,world,
                          nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                          nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                          1,0,0,FFT_PRECISION,collective_flag,gpu_aware_flag);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::deallocate()
{
  delete gc;
  gc = nullptr;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);

  memoryKK->destroy_kokkos(d_density_fft,density_fft);
  memoryKK->destroy_kokkos(d_greensfn,greensfn);
  memoryKK->destroy_kokkos(d_work1,work1);
  memoryKK->destroy_kokkos(d_work2,work2);

  delete fft1;
  fft1 = nullptr;
  delete fft2;
  fft1 = nullptr;
  delete remap;
  remap = nullptr;
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::allocate_peratom()
{
  peratom_allocate_flag = 1;

  d_u_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:u_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);

  d_v0_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v0_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_v1_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v1_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_v2_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v2_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_v3_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v3_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_v4_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v4_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);
  d_v5_brick = typename FFT_AT::t_FFT_SCALAR_3d("pppm:d_v5_brick",nzhi_out-nzlo_out+1,nyhi_out-nylo_out+1,nxhi_out-nxlo_out+1);


  // use same GC ghost grid object for peratom grid communication
  // but need to reallocate a larger gc_buf1 and gc_buf2

  npergrid = 7;

  k_gc_buf1 = FFT_DAT::tdual_FFT_SCALAR_1d("pppm:gc_buf1",npergrid*ngc_buf1);
  k_gc_buf2 = FFT_DAT::tdual_FFT_SCALAR_1d("pppm:gc_buf2",npergrid*ngc_buf2);
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::deallocate_peratom()
{
  peratom_allocate_flag = 0;
}

/* ----------------------------------------------------------------------
   estimate kspace force error for ik method
------------------------------------------------------------------------- */

template<class DeviceType>
double PPPMKokkos<DeviceType>::estimate_ik_error(double h, double prd, bigint natoms)
{
  double sum = 0.0;
  for (int m = 0; m < order; m++)
    sum += acons(order,m) * pow(h*g_ewald,2.0*m);
  double value = q2 * pow(h*g_ewald,(double)order) *
    sqrt(g_ewald*prd*sqrt(MY_2PI)*sum/natoms) / (prd*prd);

  return value;
}

/* ----------------------------------------------------------------------
   set params which determine which owned and ghost cells this proc owns
   Grid3d uses these params to partition grid
   also partition FFT grid
     n xyz lo/hi fft = FFT columns that I own (all of x dim, 2d decomp in yz)
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::set_grid_local()
{
  PPPM::set_grid_local();

  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n)
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::compute_gf_denom()
{
  int k,l,m;

  for (l = 1; l < order; l++) k_gf_b.h_view[l] = 0.0;
  k_gf_b.h_view[0] = 1.0;

  for (m = 1; m < order; m++) {
    for (l = m; l > 0; l--)
      k_gf_b.h_view[l] = 4.0 * (k_gf_b.h_view[l]*(l-m)*(l-m-0.5)-k_gf_b.h_view[l-1]*(l-m-1)*(l-m-1));
    k_gf_b.h_view[0] = 4.0 * (k_gf_b.h_view[0]*(l-m)*(l-m-0.5));
  }

  bigint ifact = 1;
  for (k = 1; k < 2*order; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < order; l++) k_gf_b.h_view[l] *= gaminv;

  k_gf_b.template modify<LMPHostType>();
  k_gf_b.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::compute_gf_ik()
{
  const double * const prd = domain->prd;

  xprd = prd[0];
  yprd = prd[1];
  const double zprd = prd[2];
  zprd_slab = zprd*slab_volfactor;
  unitkx = (MY_2PI/xprd);
  unitky = (MY_2PI/yprd);
  unitkz = (MY_2PI/zprd_slab);

  nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                          pow(-log(EPS_HOC),0.25));
  nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                          pow(-log(EPS_HOC),0.25));
  nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                          pow(-log(EPS_HOC),0.25));
  twoorder = 2*order;

  // merge three outer loops into one for better threading

  numz_fft = nzhi_fft-nzlo_fft + 1;
  numy_fft = nyhi_fft-nylo_fft + 1;
  numx_fft = nxhi_fft-nxlo_fft + 1;
  const int inum_fft = numx_fft*numy_fft*numz_fft;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_compute_gf_ik>(0,inum_fft),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_compute_gf_ik, const int &n) const
{
  int m = n/(numy_fft*numx_fft);
  int l = (n - m*numy_fft*numx_fft) / numx_fft;
  int k = n - m*numy_fft*numx_fft - l*numx_fft;
  m += nzlo_fft;
  l += nylo_fft;
  k += nxlo_fft;

  const int mper = m - nz_pppm*(2*m/nz_pppm);
  const double snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));

  const int lper = l - ny_pppm*(2*l/ny_pppm);
  const double sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));

  const int kper = k - nx_pppm*(2*k/nx_pppm);
  const double snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));

  const double sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

  if (sqk != 0.0) {
    const double numerator = 12.5663706/sqk;
    const double denominator = gf_denom(snx,sny,snz);
    double sum1 = 0.0;

    for (int nx = -nbx; nx <= nbx; nx++) {
      const double qx = unitkx*(kper+nx_pppm*nx);
      const double sx = exp(-0.25*square(qx/g_ewald));
      const double argx = 0.5*qx*xprd/nx_pppm;
      const double wx = powsinxx(argx,twoorder);

      for (int ny = -nby; ny <= nby; ny++) {
        const double qy = unitky*(lper+ny_pppm*ny);
        const double sy = exp(-0.25*square(qy/g_ewald));
        const double argy = 0.5*qy*yprd/ny_pppm;
        const double wy = powsinxx(argy,twoorder);

        for (int nz = -nbz; nz <= nbz; nz++) {
          const double qz = unitkz*(mper+nz_pppm*nz);
          const double sz = exp(-0.25*square(qz/g_ewald));
          const double argz = 0.5*qz*zprd_slab/nz_pppm;
          const double wz = powsinxx(argz,twoorder);

          const double dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
          const double dot2 = qx*qx+qy*qy+qz*qz;
          sum1 += (dot1/dot2) * sx*sy*sz * wx*wy*wz;
        }
      }
    }
    d_greensfn[n] = numerator*sum1/denominator;
  } else d_greensfn[n] = 0.0;
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
   for a triclinic system
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::compute_gf_ik_triclinic()
{
  double tmp[3];
  tmp[0] = (g_ewald/(MY_PI*nx_pppm)) * pow(-log(EPS_HOC),0.25);
  tmp[1] = (g_ewald/(MY_PI*ny_pppm)) * pow(-log(EPS_HOC),0.25);
  tmp[2] = (g_ewald/(MY_PI*nz_pppm)) * pow(-log(EPS_HOC),0.25);
  lamda2xT(&tmp[0],&tmp[0]);
  nbx = static_cast<int> (tmp[0]);
  nby = static_cast<int> (tmp[1]);
  nbz = static_cast<int> (tmp[2]);

  // Update the local copy of the domain box tilt
  h = Few<double, 6>(domain->h);
  h_inv = Few<double, 6>(domain->h_inv);

  twoorder = 2*order;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_compute_gf_ik_triclinic>(nzlo_fft,nzhi_fft+1),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_compute_gf_ik_triclinic, const int &m) const
{
  int n = (m - nzlo_fft)*(nyhi_fft+1 - nylo_fft)*(nxhi_fft+1 - nxlo_fft);

  const int mper = m - nz_pppm*(2*m/nz_pppm);
  const double snz = square(sin(MY_PI*mper/nz_pppm));

  for (int l = nylo_fft; l <= nyhi_fft; l++) {
    const int lper = l - ny_pppm*(2*l/ny_pppm);
    const double sny = square(sin(MY_PI*lper/ny_pppm));

    for (int k = nxlo_fft; k <= nxhi_fft; k++) {
      const int kper = k - nx_pppm*(2*k/nx_pppm);
      const double snx = square(sin(MY_PI*kper/nx_pppm));

      double unitk_lamda[3];
      unitk_lamda[0] = 2.0*MY_PI*kper;
      unitk_lamda[1] = 2.0*MY_PI*lper;
      unitk_lamda[2] = 2.0*MY_PI*mper;
      x2lamdaT_kokkos(&unitk_lamda[0],&unitk_lamda[0]);

      const double sqk = square(unitk_lamda[0]) + square(unitk_lamda[1]) + square(unitk_lamda[2]);

      if (sqk != 0.0) {
        const double numerator = 12.5663706/sqk;
        const double denominator = gf_denom(snx,sny,snz);
        double sum1 = 0.0;

        for (int nx = -nbx; nx <= nbx; nx++) {
          const double argx = MY_PI*kper/nx_pppm + MY_PI*nx;
          const double wx = powsinxx(argx,twoorder);

          for (int ny = -nby; ny <= nby; ny++) {
            const double argy = MY_PI*lper/ny_pppm + MY_PI*ny;
            const double wy = powsinxx(argy,twoorder);

            for (int nz = -nbz; nz <= nbz; nz++) {
              const double argz = MY_PI*mper/nz_pppm + MY_PI*nz;
              const double wz = powsinxx(argz,twoorder);

              double b[3];
              b[0] = 2.0*MY_PI*nx_pppm*nx;
              b[1] = 2.0*MY_PI*ny_pppm*ny;
              b[2] = 2.0*MY_PI*nz_pppm*nz;
              x2lamdaT_kokkos(&b[0],&b[0]);

              const double qx = unitk_lamda[0]+b[0];
              const double sx = exp(-0.25*square(qx/g_ewald));

              const double qy = unitk_lamda[1]+b[1];
              const double sy = exp(-0.25*square(qy/g_ewald));

              const double qz = unitk_lamda[2]+b[2];
              const double sz = exp(-0.25*square(qz/g_ewald));

              const double dot1 = unitk_lamda[0]*qx + unitk_lamda[1]*qy + unitk_lamda[2]*qz;
              const double dot2 = qx*qx+qy*qy+qz*qz;
              sum1 += (dot1/dot2) * sx*sy*sz * wx*wy*wz;
            }
          }
        }
        d_greensfn[n++] = numerator*sum1/denominator;
      } else d_greensfn[n++] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::particle_map()
{
  int nlocal = atomKK->nlocal;

  k_flag.h_view() = 0;
  k_flag.template modify<LMPHostType>();
  k_flag.template sync<DeviceType>();

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_particle_map>(0,nlocal),*this);
  copymode = 0;

  k_flag.template modify<DeviceType>();
  k_flag.template sync<LMPHostType>();
  if (k_flag.h_view()) error->one(FLERR,"Out of range atoms - cannot compute PPPM");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_particle_map, const int &i) const
{
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // current particle coord can be outside global and local box
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  const int nx = static_cast<int> ((x(i,0)-boxlo[0])*delxinv+shift) - OFFSET;
  const int ny = static_cast<int> ((x(i,1)-boxlo[1])*delyinv+shift) - OFFSET;
  const int nz = static_cast<int> ((x(i,2)-boxlo[2])*delzinv+shift) - OFFSET;

  d_part2grid(i,0) = nx;
  d_part2grid(i,1) = ny;
  d_part2grid(i,2) = nz;

  // check that entire stencil around nx,ny,nz will fit in my 3d brick

  if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
      ny+nlower < nylo_out || ny+nupper > nyhi_out ||
      nz+nlower < nzlo_out || nz+nupper > nzhi_out)
    k_flag.view<DeviceType>()() = 1;
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::make_rho()
{
  // clear 3d density array

  numz_out = nzhi_out-nzlo_out + 1;
  numy_out = nyhi_out-nylo_out + 1;
  numx_out = nxhi_out-nxlo_out + 1;
  const int inum_out = numz_out*numy_out*numx_out;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_make_rho_zero>(0,inum_out),*this);
  copymode = 0;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global indices of moving stencil pt

  nlocal = atomKK->nlocal;

#ifdef LMP_KOKKOS_GPU
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_make_rho_atomic>(0,nlocal),*this);
  copymode = 0;
#else
  ix = nxhi_out-nxlo_out + 1;
  iy = nyhi_out-nylo_out + 1;

  copymode = 1;
  Kokkos::TeamPolicy<DeviceType, TagPPPM_make_rho> config(lmp->kokkos->nthreads,1);
  Kokkos::parallel_for(config,*this);
  copymode = 0;
#endif
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_make_rho_zero, const int &ii) const
{
  int iz = ii/(numy_out*numx_out);
  int iy = (ii - iz*numy_out*numx_out) / numx_out;
  int ix = ii - iz*numy_out*numx_out - iy*numx_out;
  d_density_brick(iz,iy,ix) = 0.0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_make_rho_atomic, const int &i) const
{
  // The density_brick array is atomic for Half/Thread neighbor style
  Kokkos::View<FFT_SCALAR***,Kokkos::LayoutRight,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_density_brick = d_density_brick;

  int nx = d_part2grid(i,0);
  int ny = d_part2grid(i,1);
  int nz = d_part2grid(i,2);
  const FFT_SCALAR dx = nx+shiftone - (x(i,0)-boxlo[0])*delxinv;
  const FFT_SCALAR dy = ny+shiftone - (x(i,1)-boxlo[1])*delyinv;
  const FFT_SCALAR dz = nz+shiftone - (x(i,2)-boxlo[2])*delzinv;


  nz -= nzlo_out;
  ny -= nylo_out;
  nx -= nxlo_out;

  compute_rho1d(i,dx,dy,dz);

  const FFT_SCALAR z0 = delvolinv * q[i];
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
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator() (TagPPPM_make_rho, typename Kokkos::TeamPolicy<DeviceType, TagPPPM_make_rho>::member_type dev) const {
  // adapted from OPENMP/pppm.cpp:

  // determine range of grid points handled by this thread
  int tid = dev.league_rank();
  // each thread works on a fixed chunk of grid points
  const int nthreads = dev.league_size();
  const int idelta = 1 + ngrid/nthreads;
  int ifrom = tid*idelta;
  int ito = ((ifrom + idelta) > ngrid) ? ngrid : ifrom + idelta;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt

  // loop over all local atoms for all threads
  for (int i = 0; i < nlocal; i++) {

    int nx = d_part2grid(i,0);
    int ny = d_part2grid(i,1);
    int nz = d_part2grid(i,2);

    // pre-screen whether this atom will ever come within
    // reach of the data segement this thread is updating.
    if ( ((nz+nlower-nzlo_out)*ix*iy >= ito)
         || ((nz+nupper-nzlo_out+1)*ix*iy < ifrom) ) continue;

    const FFT_SCALAR dx = nx+shiftone - (x(i,0)-boxlo[0])*delxinv;
    const FFT_SCALAR dy = ny+shiftone - (x(i,1)-boxlo[1])*delyinv;
    const FFT_SCALAR dz = nz+shiftone - (x(i,2)-boxlo[2])*delzinv;

    nz -= nzlo_out;
    ny -= nylo_out;
    nx -= nxlo_out;

    compute_rho1d(i,dx,dy,dz);

    const FFT_SCALAR z0 = delvolinv * q[i];
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
          // make sure each thread only updates
          // their elements of the density grid
          if (il >= ito) break;
          if (il < ifrom) continue;
          d_density_brick(mz,my,mx) += x0*d_rho1d(i,l+order/2,0);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::brick2fft()
{
  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  numz_inout = (nzhi_in-nzlo_out)-(nzlo_in-nzlo_out) + 1;
  numy_inout = (nyhi_in-nylo_out)-(nylo_in-nylo_out) + 1;
  numx_inout = (nxhi_in-nxlo_out)-(nxlo_in-nxlo_out) + 1;
  const int inum_inout = numz_inout*numy_inout*numx_inout;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_brick2fft>(0,inum_inout),*this);
  copymode = 0;

  remap->perform(d_density_fft,d_density_fft,d_work1);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_brick2fft, const int &ii) const
{
  const int n = ii;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_density_fft[n] = d_density_brick(k,j,i);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::poisson_ik()
{
  int j;

  // transform charge density (r -> k)

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik1>(0,nfft),*this);
  copymode = 0;

  fft1->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::FORWARD);

  // global energy and virial contribution

  scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    EV_FLOAT ev;
    if (vflag_global) {
      copymode = 1;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik2>(0,nfft),*this,ev);
      copymode = 0;
      for (j = 0; j < 6; j++) virial[j] += ev.v[j];
      energy += ev.ecoul;
    } else {
      copymode = 1;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik3>(0,nfft),*this,ev);
      copymode = 0;
      energy += ev.ecoul;
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik4>(0,nfft),*this);
  copymode = 0;

  // extra FFTs for per-atom energy/virial

  if (evflag_atom) poisson_peratom();

  // triclinic system

  if (triclinic) {
    poisson_ik_triclinic();
    return;
  }

  // compute gradients of V(r) in each of 3 dims by transforming ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // merge three outer loops into one for better threading

  numz_fft = nzhi_fft-nzlo_fft + 1;
  numy_fft = nyhi_fft-nylo_fft + 1;
  numx_fft = nxhi_fft-nxlo_fft + 1;
  const int inum_fft = numz_fft*numy_fft*numx_fft;

  numz_inout = (nzhi_in-nzlo_out)-(nzlo_in-nzlo_out) + 1;
  numy_inout = (nyhi_in-nylo_out)-(nylo_in-nylo_out) + 1;
  numx_inout = (nxhi_in-nxlo_out)-(nxlo_in-nxlo_out) + 1;
  const int inum_inout = numz_inout*numy_inout*numx_inout;

  // x direction gradient

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik5>(0,inum_fft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik6>(0,inum_inout),*this);
  copymode = 0;


  // y direction gradient

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik7>(0,inum_fft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik8>(0,inum_inout),*this);
  copymode = 0;

  // z direction gradient

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik9>(0,inum_fft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik10>(0,inum_inout),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik1, const int &i) const
{
  d_work1[2*i] = d_density_fft[i];
  d_work1[2*i+1] = ZEROF;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik2, const int &i, EV_FLOAT& ev) const
{
  const double eng = s2 * d_greensfn[i] * (d_work1[2*i]*d_work1[2*i] + d_work1[2*i+1]*d_work1[2*i+1]);
  for (int j = 0; j < 6; j++) ev.v[j] += eng*d_vg(i,j);
  if (eflag_global) ev.ecoul += eng;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik3, const int &i, EV_FLOAT& ev) const
{
  ev.ecoul +=
    s2 * d_greensfn[i] * (d_work1[2*i]*d_work1[2*i] + d_work1[2*i+1]*d_work1[2*i+1]);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik4, const int &i) const
{
  d_work1[2*i] *= scaleinv * d_greensfn[i];
  d_work1[2*i+1] *= scaleinv * d_greensfn[i];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik5, const int &ii) const
{
  const int n = ii*2;
  const int k = ii/(numy_fft*numx_fft);
  const int j = (ii - k*numy_fft*numx_fft) / numx_fft;
  const int i = ii - k*numy_fft*numx_fft - j*numx_fft;
  d_work2[n] = -d_fkx[i]*d_work1[n+1];
  d_work2[n+1] = d_fkx[i]*d_work1[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik6, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdx_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik7, const int &ii) const
{
  const int n = ii*2;
  const int k = ii/(numy_fft*numx_fft);
  const int j = (ii - k*numy_fft*numx_fft) / numx_fft;
  d_work2[n] = -d_fky[j]*d_work1[n+1];
  d_work2[n+1] = d_fky[j]*d_work1[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik8, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdy_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik9, const int &ii) const
{
  const int n = ii*2;
  const int k = ii/(numy_fft*numx_fft);
  d_work2[n] = -d_fkz[k]*d_work1[n+1];
  d_work2[n+1] = d_fkz[k]*d_work1[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik10, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdz_brick(k,j,i) = d_work2[n];
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik for a triclinic system
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::poisson_ik_triclinic()
{
  // compute gradients of V(r) in each of 3 dims by transforming ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  // merge three outer loops into one for better threading

  numz_fft = nzhi_fft-nzlo_fft + 1;
  numy_fft = nyhi_fft-nylo_fft + 1;
  numx_fft = nxhi_fft-nxlo_fft + 1;

  numz_inout = (nzhi_in-nzlo_out)-(nzlo_in-nzlo_out) + 1;
  numy_inout = (nyhi_in-nylo_out)-(nylo_in-nylo_out) + 1;
  numx_inout = (nxhi_in-nxlo_out)-(nxlo_in-nxlo_out) + 1;
  const int inum_inout = numz_inout*numy_inout*numx_inout;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic1>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic2>(0,inum_inout),*this);
  copymode = 0;

  // y direction gradient

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic3>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic4>(0,inum_inout),*this);
  copymode = 0;

  // z direction gradient

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic5>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_ik_triclinic6>(0,inum_inout),*this);
  copymode = 0;

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic1, const int &ii) const
{
  d_work2[2*ii] = -d_fkx[ii]*d_work1[2*ii+1];
  d_work2[2*ii+1] = d_fkx[ii]*d_work1[2*ii];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic2, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdx_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic3, const int &ii) const
{
//  int n = (k - (nzlo_in-nzlo_out))*((nyhi_in-nylo_out) - (nylo_in-nylo_out) + 1)*((nxhi_in-nxlo_out) - (nxlo_in-nxlo_out) + 1)*2;
  d_work2[2*ii] = -d_fky[ii]*d_work1[2*ii+1];
  d_work2[2*ii+1] = d_fky[ii]*d_work1[2*ii];

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic4, const int &ii) const
{
//  int n = (k - (nzlo_in-nzlo_out))*((nyhi_in-nylo_out) - (nylo_in-nylo_out) + 1)*((nxhi_in-nxlo_out) - (nxlo_in-nxlo_out) + 1)*2;
//

  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdy_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic5, const int &ii) const
{
//  int n = (k - (nzlo_in-nzlo_out))*((nyhi_in-nylo_out) - (nylo_in-nylo_out) + 1)*((nxhi_in-nxlo_out) - (nxlo_in-nxlo_out) + 1)*2;
//
  d_work2[2*ii] = -d_fkz[ii]*d_work1[2*ii+1];
  d_work2[2*ii+1] = d_fkz[ii]*d_work1[2*ii];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_ik_triclinic6, const int &ii) const
{
//  int n = (k - (nzlo_in-nzlo_out))*((nyhi_in-nylo_out) - (nylo_in-nylo_out) + 1)*((nxhi_in-nxlo_out) - (nxlo_in-nxlo_out) + 1)*2;
//
//  for (j = nylo_in-nylo_out; j <= nyhi_in-nylo_out; j++)
//    for (i = nxlo_in-nxlo_out; i <= nxhi_in-nxlo_out; i++) {
//      d_vdz_brick(k,j,i) = d_work2[n];
//      n += 2;
//    }
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_vdz_brick(k,j,i) = d_work2[n];
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for per-atom energy/virial
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::poisson_peratom()
{
  // merge three outer loops into one for better threading

  numz_inout = (nzhi_in-nzlo_out)-(nzlo_in-nzlo_out) + 1;
  numy_inout = (nyhi_in-nylo_out)-(nylo_in-nylo_out) + 1;
  numx_inout = (nxhi_in-nxlo_out)-(nxlo_in-nxlo_out) + 1;
  const int inum_inout = numz_inout*numy_inout*numx_inout;

  // energy

  if (eflag_atom) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom1>(0,nfft),*this);
    copymode = 0;

    fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom2>(0,inum_inout),*this);
    copymode = 0;

  }

  // 6 components of virial in v0 thru v5

  if (!vflag_atom) return;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom3>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom4>(0,inum_inout),*this);
  copymode = 0;


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom5>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom6>(0,inum_inout),*this);
  copymode = 0;


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom7>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom8>(0,inum_inout),*this);
  copymode = 0;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom9>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom10>(0,inum_inout),*this);
  copymode = 0;


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom11>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom12>(0,inum_inout),*this);
  copymode = 0;


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom13>(0,nfft),*this);
  copymode = 0;

  fft2->compute(d_work2,d_work2,FFT3dKokkos<DeviceType>::BACKWARD);


  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_poisson_peratom14>(0,inum_inout),*this);
  copymode = 0;

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom1, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n];
  d_work2[n+1] = d_work1[n+1];
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom2, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_u_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom3, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,0);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,0);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom4, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v0_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom5, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,1);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,1);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom6, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v1_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom7, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,2);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,2);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom8, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v2_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom9, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,3);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,3);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom10, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v3_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom11, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,4);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,4);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom12, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v4_brick(k,j,i) = d_work2[n];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom13, const int &i) const
{
  int n = 2*i;

  d_work2[n] = d_work1[n]*d_vg(i,5);
  d_work2[n+1] = d_work1[n+1]*d_vg(i,5);
  n += 2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_poisson_peratom14, const int &ii) const
{
  const int n = ii*2;
  int k = ii/(numy_inout*numx_inout);
  int j = (ii - k*numy_inout*numx_inout) / numx_inout;
  int i = ii - k*numy_inout*numx_inout - j*numx_inout;
  k += nzlo_in-nzlo_out;
  j += nylo_in-nylo_out;
  i += nxlo_in-nxlo_out;
  d_v5_brick(k,j,i) = d_work2[n];
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::fieldforce()
{
  fieldforce_ik();
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::fieldforce_ik()
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  int nlocal = atomKK->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_fieldforce_ik>(0,nlocal),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_fieldforce_ik, const int &i) const
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

  ekx = eky = ekz = ZEROF;
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

  // convert E-field to force

  const double qfactor = qqrd2e * scale * q[i];
  f(i,0) += qfactor*ekx;
  f(i,1) += qfactor*eky;
  if (slabflag != 2) f(i,2) += qfactor*ekz;
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::fieldforce_peratom()
{
  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int nlocal = atomKK->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_fieldforce_peratom>(0,nlocal),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_fieldforce_peratom, const int &i) const
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

  nx = d_part2grid(i,0);
  ny = d_part2grid(i,1);
  nz = d_part2grid(i,2);
  dx = nx+shiftone - (x(i,0)-boxlo[0])*delxinv;
  dy = ny+shiftone - (x(i,1)-boxlo[1])*delyinv;
  dz = nz+shiftone - (x(i,2)-boxlo[2])*delzinv;

  nz -= nzlo_out;
  ny -= nylo_out;
  nx -= nxlo_out;

  compute_rho1d(i,dx,dy,dz);

  u = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
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

  if (eflag_atom) d_eatom[i] += q[i]*u;
  if (vflag_atom) {
    d_vatom(i,0) += q[i]*v0;
    d_vatom(i,1) += q[i]*v1;
    d_vatom(i,2) += q[i]*v2;
    d_vatom(i,3) += q[i]*v3;
    d_vatom(i,4) += q[i]*v4;
    d_vatom(i,5) += q[i]*v5;
  }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::pack_forward_grid_kokkos(int flag, FFT_DAT::tdual_FFT_SCALAR_1d &k_buf, int nlist, DAT::tdual_int_2d &k_list, int index)
{
  typename AT::t_int_2d_um d_list = k_list.view<DeviceType>();
  d_list_index = Kokkos::subview(d_list,index,Kokkos::ALL());
  d_buf = k_buf.view<DeviceType>();

  nx = (nxhi_out-nxlo_out+1);
  ny = (nyhi_out-nylo_out+1);

  if (flag == FORWARD_IK) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_pack_forward1>(0,nlist),*this);
    copymode = 0;
  } else if (flag == FORWARD_IK_PERATOM) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_pack_forward2>(0,nlist),*this);
    copymode = 0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_pack_forward1, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  d_buf[3*i] = d_vdx_brick(iz,iy,ix);
  d_buf[3*i+1] = d_vdy_brick(iz,iy,ix);
  d_buf[3*i+2] = d_vdz_brick(iz,iy,ix);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_pack_forward2, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  if (eflag_atom) d_buf[7*i] = d_u_brick(iz,iy,ix);
  if (vflag_atom) {
    d_buf[7*i+1] = d_v0_brick(iz,iy,ix);
    d_buf[7*i+2] = d_v1_brick(iz,iy,ix);
    d_buf[7*i+3] = d_v2_brick(iz,iy,ix);
    d_buf[7*i+4] = d_v3_brick(iz,iy,ix);
    d_buf[7*i+5] = d_v4_brick(iz,iy,ix);
    d_buf[7*i+6] = d_v5_brick(iz,iy,ix);
  }
}
/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::unpack_forward_grid_kokkos(int flag, FFT_DAT::tdual_FFT_SCALAR_1d &k_buf, int offset, int nlist, DAT::tdual_int_2d &k_list, int index)
{
  typename AT::t_int_2d_um d_list = k_list.view<DeviceType>();
  d_list_index = Kokkos::subview(d_list,index,Kokkos::ALL());
  d_buf = k_buf.view<DeviceType>();
  unpack_offset = offset;

  nx = (nxhi_out-nxlo_out+1);
  ny = (nyhi_out-nylo_out+1);

  if (flag == FORWARD_IK) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_unpack_forward1>(0,nlist),*this);
    copymode = 0;
  } else if (flag == FORWARD_IK_PERATOM) {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_unpack_forward2>(0,nlist),*this);
    copymode = 0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_unpack_forward1, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  d_vdx_brick(iz,iy,ix) = d_buf[3*i   + unpack_offset];
  d_vdy_brick(iz,iy,ix) = d_buf[3*i+1 + unpack_offset];
  d_vdz_brick(iz,iy,ix) = d_buf[3*i+2 + unpack_offset];
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_unpack_forward2, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  if (eflag_atom) d_u_brick(iz,iy,ix) = d_buf[7*i];
  if (vflag_atom) {
    d_v0_brick(iz,iy,ix) = d_buf[7*i+1 + unpack_offset];
    d_v1_brick(iz,iy,ix) = d_buf[7*i+2 + unpack_offset];
    d_v2_brick(iz,iy,ix) = d_buf[7*i+3 + unpack_offset];
    d_v3_brick(iz,iy,ix) = d_buf[7*i+4 + unpack_offset];
    d_v4_brick(iz,iy,ix) = d_buf[7*i+5 + unpack_offset];
    d_v5_brick(iz,iy,ix) = d_buf[7*i+6 + unpack_offset];
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::pack_reverse_grid_kokkos(int /*flag*/, FFT_DAT::tdual_FFT_SCALAR_1d &k_buf, int nlist, DAT::tdual_int_2d &k_list, int index)
{
  typename AT::t_int_2d_um d_list = k_list.view<DeviceType>();
  d_list_index = Kokkos::subview(d_list,index,Kokkos::ALL());
  d_buf = k_buf.view<DeviceType>();

  nx = (nxhi_out-nxlo_out+1);
  ny = (nyhi_out-nylo_out+1);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_pack_reverse>(0,nlist),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_pack_reverse, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  d_buf[i] = d_density_brick(iz,iy,ix);
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::unpack_reverse_grid_kokkos(int /*flag*/, FFT_DAT::tdual_FFT_SCALAR_1d &k_buf, int offset, int nlist, DAT::tdual_int_2d &k_list, int index)
{
  typename AT::t_int_2d_um d_list = k_list.view<DeviceType>();
  d_list_index = Kokkos::subview(d_list,index,Kokkos::ALL());
  d_buf = k_buf.view<DeviceType>();
  unpack_offset = offset;

  nx = (nxhi_out-nxlo_out+1);
  ny = (nyhi_out-nylo_out+1);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_unpack_reverse>(0,nlist),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_unpack_reverse, const int &i) const
{
  const double dlist = (double) d_list_index[i];
  const int iz = (int) (dlist/(nx*ny));
  const int iy = (int) ((dlist - iz*nx*ny)/nx);
  const int ix = d_list_index[i] - iz*nx*ny - iy*nx;
  d_density_brick(iz,iy,ix) += d_buf[i + unpack_offset];
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::compute_rho1d(const int i, const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                         const FFT_SCALAR &dz) const
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-1; l >= 0; l--) {
      r1 = d_rho_coeff(l,k-(1-order)/2) + r1*dx;
      r2 = d_rho_coeff(l,k-(1-order)/2) + r2*dy;
      r3 = d_rho_coeff(l,k-(1-order)/2) + r3*dz;
    }
    d_rho1d(i,k+order/2,0) = r1;
    d_rho1d(i,k+order/2,1) = r2;
    d_rho1d(i,k+order/2,2) = r3;
  }
}

/* ----------------------------------------------------------------------
   generate coeffients for the weight function of order n

              (n-1)
  Wn(x) =     Sum    wn(k,x) , Sum is over every other integer
           k=-(n-1)
  For k=-(n-1),-(n-1)+2, ....., (n-1)-2,n-1
      k is odd integers if n is even and even integers if n is odd
              ---
             | n-1
             | Sum a(l,j)*(x-k/2)**l   if abs(x-k/2) < 1/2
  wn(k,x) = <  l=0
             |
             |  0                       otherwise
              ---
  a coeffients are packed into the array rho_coeff to eliminate zeros
  rho_coeff(l,((k+mod(n+1,2))/2) = a(l,k)
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::compute_rho_coeff()
{
  int j,k,l,m;
  FFT_SCALAR s;
  FFT_SCALAR **a = new FFT_SCALAR *[order];
  for (int i = 0; i < order; ++i)
    a[i] = new FFT_SCALAR[2*order+1];

  for (k = 0; k <= 2*order; k++)
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;

  a[0][order] = 1.0;
  for (j = 1; j < order; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
        a[l+1][k+order] = (a[l][k+1+order]-a[l][k-1+order]) / (l+1);
#ifdef FFT_SINGLE
        s += powf(0.5,(float) l+1) *
          (a[l][k-1+order] + powf(-1.0,(float) l) * a[l][k+1+order]) / (l+1);
#else
        s += pow(0.5,(double) l+1) *
          (a[l][k-1+order] + pow(-1.0,(double) l) * a[l][k+1+order]) / (l+1);
#endif
      }
      a[0][k+order] = s;
    }
  }

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      h_rho_coeff(l,m-(1-order)/2) = a[l][k+order];
    m++;
  }
  for (int i = 0; i < order; ++i)
    delete[] a[i];
  delete[] a;
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

template<class DeviceType>
void PPPMKokkos<DeviceType>::slabcorr()
{
  // compute local contribution to global dipole moment

  zprd_slab = domain->zprd*slab_volfactor;
  int nlocal = atomKK->nlocal;

  double dipole = 0.0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPM_slabcorr1>(0,nlocal),*this,dipole);
  copymode = 0;

  // sum local contributions to get global dipole moment

  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPPPM_slabcorr2>(0,nlocal),*this,dipole_r2);
    copymode = 0;

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    efact = qscale * MY_2PI/volume;
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_slabcorr3>(0,nlocal),*this);
    copymode = 0;
  }

  // add on force corrections

  ffact = qscale * (-4.0*MY_PI/volume);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_slabcorr4>(0,nlocal),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_slabcorr1, const int &i, double &dipole) const
{
  dipole += q[i]*x(i,2);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_slabcorr2, const int &i, double &dipole_r2) const
{
  dipole_r2 += q[i]*x(i,2)*x(i,2);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_slabcorr3, const int &i) const
{
  d_eatom[i] += efact * q[i]*(x(i,2)*dipole_all - 0.5*(dipole_r2 +
    qsum*x(i,2)*x(i,2)) - qsum*zprd_slab*zprd_slab/12.0);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_slabcorr4, const int &i) const
{
  f(i,2) += ffact * q[i]*(dipole_all - qsum*x(i,2));
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

template<class DeviceType>
int PPPMKokkos<DeviceType>::timing_1d(int n, double &time1d)
{
  double time1,time2;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_timing_zero>(0,2*nfft_both),*this);
  copymode = 0;

  MPI_Barrier(world);
  time1 = platform::walltime();

  for (int i = 0; i < n; i++) {
    fft1->timing1d(d_work1,nfft_both,FFT3dKokkos<DeviceType>::FORWARD);
    fft2->timing1d(d_work1,nfft_both,FFT3dKokkos<DeviceType>::BACKWARD);
    fft2->timing1d(d_work1,nfft_both,FFT3dKokkos<DeviceType>::BACKWARD);
    fft2->timing1d(d_work1,nfft_both,FFT3dKokkos<DeviceType>::BACKWARD);
  }

  MPI_Barrier(world);
  time2 = platform::walltime();
  time1d = time2 - time1;

  return 4;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PPPMKokkos<DeviceType>::operator()(TagPPPM_timing_zero, const int &i) const
{
  d_work1[i] = ZEROF;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

template<class DeviceType>
int PPPMKokkos<DeviceType>::timing_3d(int n, double &time3d)
{
  double time1,time2;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPPPM_timing_zero>(0,2*nfft_both),*this);
  copymode = 0;

  MPI_Barrier(world);
  time1 = platform::walltime();

  for (int i = 0; i < n; i++) {
    fft1->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::FORWARD);
    fft2->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::BACKWARD);
    fft2->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::BACKWARD);
    fft2->compute(d_work1,d_work1,FFT3dKokkos<DeviceType>::BACKWARD);
  }

  MPI_Barrier(world);
  time2 = platform::walltime();
  time3d = time2 - time1;

  return 4;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double PPPMKokkos<DeviceType>::memory_usage()
{
  double bytes = (double)nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  bytes += (double)4 * nbrick * sizeof(FFT_SCALAR);
  if (triclinic) bytes += (double)3 * nfft_both * sizeof(double);
  bytes += (double)6 * nfft_both * sizeof(double);
  bytes += (double)nfft_both * sizeof(double);
  bytes += (double)nfft_both*5 * sizeof(FFT_SCALAR);

  if (peratom_allocate_flag)
    bytes += (double)6 * nbrick * sizeof(FFT_SCALAR);

  // two Grid3d bufs

  bytes += (double)(ngc_buf1 + ngc_buf2) * npergrid * sizeof(FFT_SCALAR);

  return bytes;
}

namespace LAMMPS_NS {
template class PPPMKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PPPMKokkos<LMPHostType>;
#endif
}

