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
   Contributing authors: Jiuyang Liang, Libin Lu, Shidong Jiang (Flatiron)
     analytic diff (2 FFT) option is set to default and recommended
------------------------------------------------------------------------- */

#include "esp.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "math_pswf.h"
#include "remap_wrap.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr int MAXORDER = 32;
static constexpr int OFFSET = 16384;
static constexpr double SMALL = 0.00001;
static constexpr double EPS_HOC = 1.0e-7;
static constexpr FFT_SCALAR ZEROF = 0.0;

/* ---------------------------------------------------------------------- */

ESP::ESP(LAMMPS *lmp) : KSpace(lmp),
  factors(nullptr), density_brick(nullptr), vdx_brick(nullptr), vdy_brick(nullptr), vdz_brick(nullptr),
  u_brick(nullptr), v0_brick(nullptr), v1_brick(nullptr), v2_brick(nullptr), v3_brick(nullptr),
  v4_brick(nullptr), v5_brick(nullptr), greensfn(nullptr), greensfn2(nullptr), vg(nullptr), vg2(nullptr), fkx(nullptr), fky(nullptr),
  fkz(nullptr), density_fft(nullptr), work1(nullptr), work2(nullptr), rho1d(nullptr),
  rho_coeff(nullptr), drho1d(nullptr), drho_coeff(nullptr),
  sf_precoeff1(nullptr), sf_precoeff2(nullptr), sf_precoeff3(nullptr),
  sf_precoeff4(nullptr), sf_precoeff5(nullptr), sf_precoeff6(nullptr),
  fft1(nullptr), fft2(nullptr), remap(nullptr), gc(nullptr),
  gc_buf1(nullptr), gc_buf2(nullptr), density_A_brick(nullptr), density_B_brick(nullptr), density_A_fft(nullptr),
  density_B_fft(nullptr), part2grid(nullptr), boxlo(nullptr)
{
  peratom_allocate_flag = 0;
  group_allocate_flag = 0;

  espflag = 1;
  group_group_enable = 1;
  triclinic = domain->triclinic;

  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nfft_both = 0;
  nxhi_in = nxlo_in = nxhi_out = nxlo_out = 0;
  nyhi_in = nylo_in = nyhi_out = nylo_out = 0;
  nzhi_in = nzlo_in = nzhi_out = nzlo_out = 0;

  density_brick = vdx_brick = vdy_brick = vdz_brick = nullptr;
  density_fft = nullptr;
  u_brick = nullptr;
  v0_brick = v1_brick = v2_brick = v3_brick = v4_brick = v5_brick = nullptr;
  greensfn = nullptr;
  greensfn2 = nullptr;

  work1 = work2 = nullptr;
  vg = nullptr;
  vg2 = nullptr;
  fkx = fky = fkz = nullptr;

  sf_precoeff1 = sf_precoeff2 = sf_precoeff3 =
    sf_precoeff4 = sf_precoeff5 = sf_precoeff6 = nullptr;

  density_A_brick = density_B_brick = nullptr;
  density_A_fft = density_B_fft = nullptr;

  rho1d = rho_coeff = drho1d = drho_coeff = nullptr;

  fft1 = fft2 = nullptr;
  remap = nullptr;
  gc = nullptr;
  gc_buf1 = gc_buf2 = nullptr;

  nmax = 0;
  part2grid = nullptr;

  differentiation_flag = 1; // default to analytic differentiation
}

/* ---------------------------------------------------------------------- */

void ESP::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Illegal kspace_style {} command", force->kspace_style);

  accuracy_relative = fabs(utils::numeric(FLERR, arg[0], false, lmp)); // splitting accuracy

  if (narg == 1) {
    // Use the single kspace_style tolerance as an overall force-accuracy target.
    // A tighter default spread budget keeps the table-based stencil order accurate.
    spreading_accuracy = 0.25 * accuracy_relative;
  } else {
    spreading_accuracy = fabs(utils::numeric(FLERR, arg[1], false, lmp));
  }

  if (accuracy_relative <= 0.0 || spreading_accuracy <= 0.0 ||
      accuracy_relative > 1.0 || spreading_accuracy > 1.0)
    error->all(FLERR, "Invalid relative accuracy {:g} or spreading accuracy {:g} for kspace_style {}",
               accuracy_relative, spreading_accuracy, force->kspace_style);

  // Estimate order for ESP method
  order = estimate_order(accuracy_relative);
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

ESP::~ESP()
{
  if (copymode) return;

  delete[] factors;
  ESP::deallocate();
  if (peratom_allocate_flag) ESP::deallocate_peratom();
  if (group_allocate_flag) ESP::deallocate_groups();
  memory->destroy(part2grid);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void ESP::init()
{
  if (me == 0) utils::logmesg(lmp,"ESP initialization ...\n");

  // error check

  triclinic_check();

  if (triclinic != domain->triclinic)
    error->all(FLERR,"Must redefine kspace_style after changing to triclinic box");

  if (domain->triclinic && differentiation_flag == 1)
    error->all(FLERR,"Cannot (yet) use ESP with triclinic box and kspace_modify diff ad");
  if (domain->triclinic && slabflag)
    error->all(FLERR,"Cannot (yet) use ESP with triclinic box and slab correction");
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use ESP with 2d simulation");

  if (!atom->q_flag)
    error->all(FLERR,"Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with ESP");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab ESP");
  }

  if (order < 2 || order > MAXORDER)
    error->all(FLERR,"ESP order cannot be < 2 or > {}",MAXORDER);

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  pair_check();

  int itmp = 0;
  auto *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  if (p_cutoff == nullptr)
    error->all(FLERR, "KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  qdist = 0.0;

  if (tip4pflag) {
    if (me == 0) utils::logmesg(lmp,"  extracting TIP4P info from pair style\n");

    auto *p_qdist = (double *) force->pair->extract("qdist",itmp);
    auto *p_typeO = (int *) force->pair->extract("typeO",itmp);
    auto *p_typeH = (int *) force->pair->extract("typeH",itmp);
    auto *p_typeA = (int *) force->pair->extract("typeA",itmp);
    auto *p_typeB = (int *) force->pair->extract("typeB",itmp);
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all(FLERR,"Pair style is incompatible with TIP4P KSpace style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == nullptr || force->bond == nullptr ||
        force->angle->setflag == nullptr || force->bond->setflag == nullptr)
      error->all(FLERR,"Bond and angle potentials must be defined for TIP4P");
    if (typeA < 1 || typeA > atom->nangletypes ||
        force->angle->setflag[typeA] == 0)
      error->all(FLERR,"Bad TIP4P angle type for ESP/TIP4P");
    if (typeB < 1 || typeB > atom->nbondtypes ||
        force->bond->setflag[typeB] == 0)
      error->all(FLERR,"Bad TIP4P bond type for ESP/TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5*theta) * blen);
  }

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
  if (group_allocate_flag) deallocate_groups();

  // setup FFT grid resolution
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  gc = nullptr;
  int iteration = 0;

  MathPSWF::prolc180(accuracy_relative, select_c);
  MathPSWF::prolc180(spreading_accuracy, spreading_select_c);

  // Here debug
  while (order >= minorder) {
    if (iteration && me == 0)
      error->warning(FLERR,"Reducing ESP order b/c stencil extends "
                     "beyond nearest neighbor processor");

    set_grid_global();
    set_grid_local();
    if (overlap_allowed) break;

    gc = new Grid3d(lmp,world,nx_pppm,ny_pppm,nz_pppm);
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

    if (me == 0)
      error->warning(FLERR,"Order decreased to {} to fit stencil on nearest neighbor processors. May not achieve the desired accuracy.", order);
  }

  if (order < minorder) error->all(FLERR,"ESP order < minimum allowed order");
  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR,"ESP grid stencil extends beyond nearest neighbor processor");
  if (gc) delete gc;

  // allocate K-space dependent memory
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();

  // pre-compute tables for splitting/spreading function

  build_table(spreading_accuracy);
  if (differentiation_flag == 1) compute_sf_precoeff();

  // print stats

  int ngrid_max,nfft_both_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    std::string mesg = fmt::format("  Spreading parameter c = {:.8g}\n",spreading_select_c);
    mesg += fmt::format("  Splitting parameter c = {:.8g}\n",select_c);
    mesg += fmt::format("  grid = {} {} {}\n",nx_pppm,ny_pppm,nz_pppm);
    mesg += fmt::format("  stencil order = {}\n",order);
    mesg += fmt::format("  estimated relative splitting force accuracy = {:.8g}\n",
                       accuracy_relative);
    mesg += fmt::format("  estimated relative spreading accuracy = {:.8g}\n",
                       spreading_accuracy);
    mesg += "  using " LMP_FFT_PREC " precision " LMP_FFT_LIB "\n";
    mesg += fmt::format("  3d grid and FFT values/proc = {} {}\n",
                       ngrid_max,nfft_both_max);
    utils::logmesg(lmp,mesg);
  }

}

/* ----------------------------------------------------------------------
   adjust ESP coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void ESP::setup()
{
  if (triclinic) {
    setup_triclinic();
    return;
  }

  // perform some checks to avoid illegal boundaries with read_data

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with ESP");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab ESP");
  }

  int i,j,k,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab ESP
  // z dimension for 3d ESP is zprd since slab_volfactor = 1.0

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

  double unitkx = (MY_2PI/xprd);
  double unitky = (MY_2PI/yprd);
  double unitkz = (MY_2PI/zprd_slab);

  // fkx,fky,fkz for my FFT grid pts

  double per;

  for (i = nxlo_fft; i <= nxhi_fft; i++) {
    per = i - nx_pppm*(2*i/nx_pppm);
    fkx[i] = unitkx*per;
  }

  for (i = nylo_fft; i <= nyhi_fft; i++) {
    per = i - ny_pppm*(2*i/ny_pppm);
    fky[i] = unitky*per;
  }

  for (i = nzlo_fft; i <= nzhi_fft; i++) {
    per = i - nz_pppm*(2*i/nz_pppm);
    fkz[i] = unitkz*per;
  }

  // virial coefficients

  double sqk,vterm;

    n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        sqk = fkx[i]*fkx[i] + fky[j]*fky[j] + fkz[k]*fkz[k];
        if (sqk == 0.0) {
          vg[n][0] = 0.0;
          vg[n][1] = 0.0;
          vg[n][2] = 0.0;
          vg[n][3] = 0.0;
          vg[n][4] = 0.0;
          vg[n][5] = 0.0;
          vg2[n][0] = 0.0;
          vg2[n][1] = 0.0;
          vg2[n][2] = 0.0;
          vg2[n][3] = 0.0;
          vg2[n][4] = 0.0;
          vg2[n][5] = 0.0;
        } else {
          vterm = -2.0 / sqk;
          vg[n][0] = 1.0 + vterm*fkx[i]*fkx[i];
          vg[n][1] = 1.0 + vterm*fky[j]*fky[j];
          vg[n][2] = 1.0 + vterm*fkz[k]*fkz[k];
          vg[n][3] = vterm*fkx[i]*fky[j];
          vg[n][4] = vterm*fkx[i]*fkz[k];
          vg[n][5] = vterm*fky[j]*fkz[k];
          vg2[n][0] = fkx[i]*fkx[i];
          vg2[n][1] = fky[j]*fky[j];
          vg2[n][2] = fkz[k]*fkz[k];
          vg2[n][3] = fkx[i]*fky[j];
          vg2[n][4] = fkx[i]*fkz[k];
          vg2[n][5] = fky[j]*fkz[k];
        }
        n++;
      }
    }
  }

  if (differentiation_flag == 1) compute_gf_ad();
  else compute_gf_ik();
}

/* ----------------------------------------------------------------------
   adjust ESP coeffs, called initially and whenever volume has changed
   for a triclinic system
------------------------------------------------------------------------- */

void ESP::setup_triclinic()
{
  int i,j,k,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab ESP
  // z dimension for 3d ESP is zprd since slab_volfactor = 1.0

  prd = domain->prd;

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

  // fkx,fky,fkz for my FFT grid pts

  double per_i,per_j,per_k;

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    per_k = k - nz_pppm*(2*k/nz_pppm);
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      per_j = j - ny_pppm*(2*j/ny_pppm);
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        per_i = i - nx_pppm*(2*i/nx_pppm);

        double unitk_lamda[3];
        unitk_lamda[0] = 2.0*MY_PI*per_i;
        unitk_lamda[1] = 2.0*MY_PI*per_j;
        unitk_lamda[2] = 2.0*MY_PI*per_k;
        x2lamdaT(&unitk_lamda[0],&unitk_lamda[0]);
        fkx[n] = unitk_lamda[0];
        fky[n] = unitk_lamda[1];
        fkz[n] = unitk_lamda[2];
        n++;
      }
    }
  }

  // virial coefficients

  double sqk,vterm;

  for (n = 0; n < nfft; n++) {
    sqk = fkx[n]*fkx[n] + fky[n]*fky[n] + fkz[n]*fkz[n];
    if (sqk == 0.0) {
      vg[n][0] = 0.0;
      vg[n][1] = 0.0;
      vg[n][2] = 0.0;
      vg[n][3] = 0.0;
      vg[n][4] = 0.0;
      vg[n][5] = 0.0;
      vg2[n][0] = 0.0;
      vg2[n][1] = 0.0;
      vg2[n][2] = 0.0;
      vg2[n][3] = 0.0;
      vg2[n][4] = 0.0;
      vg2[n][5] = 0.0;
    } else {
      vterm = -2.0 / sqk;
      vg[n][0] = 1.0 + vterm*fkx[n]*fkx[n];
      vg[n][1] = 1.0 + vterm*fky[n]*fky[n];
      vg[n][2] = 1.0 + vterm*fkz[n]*fkz[n];
      vg[n][3] = vterm*fkx[n]*fky[n];
      vg[n][4] = vterm*fkx[n]*fkz[n];
      vg[n][5] = vterm*fky[n]*fkz[n];
      vg2[n][0] = fkx[n]*fkx[n];
      vg2[n][1] = fky[n]*fky[n];
      vg2[n][2] = fkz[n]*fkz[n];
      vg2[n][3] = fkx[n]*fky[n];
      vg2[n][4] = fkx[n]*fkz[n];
      vg2[n][5] = fky[n]*fkz[n];
    }
  }

  compute_gf_ik_triclinic();
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void ESP::reset_grid()
{
  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  if (group_allocate_flag) deallocate_groups();

  // reset splitting/spreading c
  MathPSWF::prolc180(accuracy_relative, select_c);
  MathPSWF::prolc180(spreading_accuracy, spreading_select_c);

  // reset portion of global grid that each proc owns

  set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();

  ////
  build_table(spreading_accuracy);

  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR,"ESP grid stencil extends beyond nearest neighbor processor");

  // pre-compute 1d charge distribution coefficients

  if (differentiation_flag == 1) compute_sf_precoeff();

  ////compute_rho_coeff();

  // pre-compute volume-dependent coeffs for portion of grid I now own

  setup();
}

/* ----------------------------------------------------------------------
   compute the ESP long-range force, energy, virial
------------------------------------------------------------------------- */

void ESP::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"esp:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM,6,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  }

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    double self_coeff = 0.00;

    for(int i=1;i<num_of_energy_poly;i++)
    {
      self_coeff += 2 * (i+0.00) * energy_poly_coeff[i] * (i%2==1?1.0:-1.0);
    }

    energy *= 0.5*volume;
    //energy -= (-energy_poly_coeff[1]/cutoff) * qsqsum / 2.0;
    energy -= (-self_coeff/cutoff) * qsqsum / 2.0;
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
    double *q = atom->q;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      double self_coeff = 0.00;
      for(int i=1;i<num_of_energy_poly;i++)
      {
        self_coeff += 2 * (i+0.00) * energy_poly_coeff[i] * (i%2==1?1.0:-1.0);
      }
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        //eatom[i] -= (-energy_poly_coeff[1]/cutoff) * q[i] * q[i] / 2.0;
        eatom[i] -= (-self_coeff/cutoff) * q[i] * q[i] / 2.0;
        eatom[i] *= qscale;
      }
      for (i = nlocal; i < ntotal; i++) eatom[i] *= 0.5*qscale;
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void ESP::allocate()
{
  // create ghost grid object for rho and electric field communication
  // returns local owned and ghost grid bounds
  // setup communication patterns and buffers

  gc = new Grid3d(lmp,world,nx_pppm,ny_pppm,nz_pppm);
  gc->set_distance(0.5*neighbor->skin + qdist);
  gc->set_stencil_atom(-nlower,nupper);
  gc->set_shift_atom(shiftatom_lo,shiftatom_hi);
  gc->set_zfactor(slab_volfactor);

  gc->setup_grid(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                 nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

  gc->setup_comm(ngc_buf1,ngc_buf2);

  if (differentiation_flag) npergrid = 1;
  else npergrid = 3;

  memory->create(gc_buf1,npergrid*ngc_buf1,"esp:gc_buf1");
  memory->create(gc_buf2,npergrid*ngc_buf2,"esp:gc_buf2");

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

  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:density_brick");

  memory->create(density_fft,nfft_both,"esp:density_fft");
  memory->create(greensfn,nfft_both,"esp:greensfn");
  memory->create(greensfn2,nfft_both,"esp:greensfn2");
  memory->create(work1,2*nfft_both,"esp:work1");
  memory->create(work2,2*nfft_both,"esp:work2");
  memory->create(vg,nfft_both,6,"esp:vg");
  memory->create(vg2,nfft_both,6,"esp:vg2");

  if (triclinic == 0) {
    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"esp:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"esp:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"esp:fkz");
  } else {
    memory->create(fkx,nfft_both,"esp:fkx");
    memory->create(fky,nfft_both,"esp:fky");
    memory->create(fkz,nfft_both,"esp:fkz");
  }

  if (differentiation_flag == 1) {
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:u_brick");

    memory->create(sf_precoeff1,nfft_both,"esp:sf_precoeff1");
    memory->create(sf_precoeff2,nfft_both,"esp:sf_precoeff2");
    memory->create(sf_precoeff3,nfft_both,"esp:sf_precoeff3");
    memory->create(sf_precoeff4,nfft_both,"esp:sf_precoeff4");
    memory->create(sf_precoeff5,nfft_both,"esp:sf_precoeff5");
    memory->create(sf_precoeff6,nfft_both,"esp:sf_precoeff6");

  } else {
    memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"esp:vdx_brick");
    memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"esp:vdy_brick");
    memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"esp:vdz_brick");
  }

  // summation coeffs

  order_allocated = order;
  memory->create2d_offset(rho1d,3,-order/2,order/2,"esp:rho1d");
  memory->create2d_offset(drho1d,3,-order/2,order/2,"esp:drho1d");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decomposition
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp,collective_flag,nonblocking_flag);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp,collective_flag,nonblocking_flag);

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,FFT_PRECISION,collective_flag,nonblocking_flag);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void ESP::deallocate()
{
  delete gc;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);

  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);

  if (differentiation_flag == 1) {
    memory->destroy3d_offset(u_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy(sf_precoeff1);
    memory->destroy(sf_precoeff2);
    memory->destroy(sf_precoeff3);
    memory->destroy(sf_precoeff4);
    memory->destroy(sf_precoeff5);
    memory->destroy(sf_precoeff6);
  } else {
    memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);
  }

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(greensfn2);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);
  memory->destroy(vg2);

  if (triclinic == 0) {
    memory->destroy1d_offset(fkx,nxlo_fft);
    memory->destroy1d_offset(fky,nylo_fft);
    memory->destroy1d_offset(fkz,nzlo_fft);
  } else {
    memory->destroy(fkx);
    memory->destroy(fky);
    memory->destroy(fkz);
  }

  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(drho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff, 0, (1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff, 0, (1-order_allocated)/2);

  if(force_poly_coeff)
    memory->destroy(force_poly_coeff);

  if(energy_poly_coeff)
    memory->destroy(energy_poly_coeff);

  if(fourier_split_poly_coeff)
    memory->destroy(fourier_split_poly_coeff);

  if(fourier_spread_poly_coeff)
    memory->destroy(fourier_spread_poly_coeff);

  delete fft1;
  delete fft2;
  delete remap;
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void ESP::allocate_peratom()
{
  peratom_allocate_flag = 1;

  if (differentiation_flag != 1)
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"esp:u_brick");

  memory->create3d_offset(v0_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v0_brick");

  memory->create3d_offset(v1_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v1_brick");
  memory->create3d_offset(v2_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v2_brick");
  memory->create3d_offset(v3_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v3_brick");
  memory->create3d_offset(v4_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v4_brick");
  memory->create3d_offset(v5_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:v5_brick");

  // use same GC ghost grid object for peratom grid communication
  // but need to reallocate a larger gc_buf1 and gc_buf2

  if (differentiation_flag) npergrid = 6;
  else npergrid = 7;

  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
  memory->create(gc_buf1,npergrid*ngc_buf1,"esp:gc_buf1");
  memory->create(gc_buf2,npergrid*ngc_buf2,"esp:gc_buf2");
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void ESP::deallocate_peratom()
{
  peratom_allocate_flag = 0;

  memory->destroy3d_offset(v0_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v1_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v2_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v3_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v4_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(v5_brick,nzlo_out,nylo_out,nxlo_out);

  if (differentiation_flag != 1)
    memory->destroy3d_offset(u_brick,nzlo_out,nylo_out,nxlo_out);
}

/* ----------------------------------------------------------------------
   set global size of ESP grid = nx,ny,nz_pppm
   used for charge accumulation, FFTs, and electric field interpolation
------------------------------------------------------------------------- */

void ESP::set_grid_global()
{
  // use xprd,yprd,zprd (even if triclinic, and then scale later)
  // adjust z dimension for 2d slab ESP
  // 3d ESP just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {

    if (differentiation_flag == 1) {

      h_x = h_y = h_z = MY_PI * cutoff / select_c; // set initial h

      nx_pppm = static_cast<int>(ceil(xprd/h_x));
      ny_pppm = static_cast<int>(ceil(yprd/h_y));
      nz_pppm = static_cast<int>(ceil(zprd_slab/h_z));

      if (nx_pppm <= 1) nx_pppm = 2;
      if (ny_pppm <= 1) ny_pppm = 2;
      if (nz_pppm <= 1) nz_pppm = 2;

    } else {

      h_x = h_y = h_z = MY_PI * cutoff / select_c;

      nx_pppm = static_cast<int>(xprd/h_x) + 1;
      ny_pppm = static_cast<int>(yprd/h_y) + 1;
      nz_pppm = static_cast<int>(zprd_slab/h_z) + 1;
    }

    // scale grid for triclinic skew

    if (triclinic) {
      double tmp[3];
      tmp[0] = nx_pppm/xprd;
      tmp[1] = ny_pppm/yprd;
      tmp[2] = nz_pppm/zprd;
      lamda2xT(&tmp[0],&tmp[0]);
      nx_pppm = static_cast<int>(tmp[0]) + 1;
      ny_pppm = static_cast<int>(tmp[1]) + 1;
      nz_pppm = static_cast<int>(tmp[2]) + 1;
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;

  if (triclinic == 0) {
    h_x = xprd/nx_pppm;
    h_y = yprd/ny_pppm;
    h_z = zprd_slab/nz_pppm;
  } else {
    double tmp[3];
    tmp[0] = nx_pppm;
    tmp[1] = ny_pppm;
    tmp[2] = nz_pppm;
    x2lamdaT(&tmp[0],&tmp[0]);
    h_x = 1.0/tmp[0];
    h_y = 1.0/tmp[1];
    h_z = 1.0/tmp[2];
  }

  if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
    error->all(FLERR,"ESP grid is too large");
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int ESP::factorable(int n)
{
  int i;

  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        break;
      }
    }
    if (i == nfactors) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   set params which determine which owned and ghost cells this proc owns
   Grid3d uses these params to partition grid
   also partition FFT grid
     n xyz lo/hi fft = FFT columns that I own (all of x dim, 2d decomp in yz)
------------------------------------------------------------------------- */

void ESP::set_grid_local()
{
  // shift values for particle <-> grid mapping depend on stencil order
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1
  // used in particle_map() and make_rho() and fieldforce()

  if (order % 2) shift = OFFSET + 0.5;
  else shift = OFFSET;

  if (order % 2) shiftone = 0.0;
  else shiftone = 0.5;

  // nlower/nupper = stencil size for mapping particles to grid

  nlower = -(order-1)/2;
  nupper = order/2;

  // shiftatom lo/hi are passed to Grid3d to determine ghost cell extents
  // shiftatom_lo = min shift on lo side
  // shiftatom_hi = max shift on hi side

  if (order % 2) {
    shiftatom_lo = 0.5;
    shiftatom_hi = 0.5;
  }
  else if (order % 2 == 0) {
    shiftatom_lo = 0.0;
    shiftatom_hi = 0.0;
  }

  // x-pencil decomposition of FFT mesh
  // global indices range from 0 to N-1
  // each proc owns entire x-dimension, clumps of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  //   NOTE: commented out lines support this
  //     need to ensure fft3d.cpp and remap.cpp support 2D planes
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own in x-pencil decomposition

  int npey_fft,npez_fft;
  procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx_pppm - 1;
  nylo_fft = me_y*ny_pppm/npey_fft;
  nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
  nzlo_fft = me_z*nz_pppm/npez_fft;
  nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void ESP::compute_gf_ik()
{
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz;
  double wx,wy,wz,qx,qy,qz;
  double sum1,dot1,dot2;
  double numerator,denominator;
  double sqk;

  int k,l,m,n,nx,ny,nz,kper,lper,mper;

  const int nbx = static_cast<int>(static_cast<int>(select_c * xprd / (MY_2PI * cutoff * nx_pppm)) * pow(-log(EPS_HOC),0.25));
  const int nby = static_cast<int>(static_cast<int>(select_c * yprd / (MY_2PI * cutoff * ny_pppm)) * pow(-log(EPS_HOC),0.25));
  const int nbz = static_cast<int>(static_cast<int>(select_c * zprd / (MY_2PI * cutoff * nz_pppm)) * pow(-log(EPS_HOC),0.25));

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = unitkz*mper;

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = unitky*lper;

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = unitkx*kper;

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          numerator = 1.0/sqk;
          denominator = gf_denom_psw(snx,sny,snz,xprd/nx_pppm,yprd/ny_pppm,zprd_slab/nz_pppm);

          sum1 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);

            double ph_2_kx_c = order * (xprd/nx_pppm) / 2.0 * fabs(qx) / spreading_select_c;
            wx = 0.00;
            if(ph_2_kx_c <= 1.00){
                ph_2_kx_c = 2.0 * ph_2_kx_c - 1.0;
                double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
                double Fourier_spreading_r = 1.0;
                for(int i=1; i<fourier_spreading_order; i++){
                    Fourier_spreading_r *= ph_2_kx_c;
                    Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
                }
                wx = order / 2.0 * Fourier_spreading_appx;
                wx = wx * wx;
            }

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);

              double ph_2_ky_c = order * (yprd/ny_pppm) / 2.0 * fabs(qy) / spreading_select_c;
              wy = 0.00;
              if(ph_2_ky_c <= 1.00){
                ph_2_ky_c = 2.0 * ph_2_ky_c - 1.0;
                double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
                double Fourier_spreading_r = 1.0;
                for(int i=1; i<fourier_spreading_order; i++){
                    Fourier_spreading_r *= ph_2_ky_c;
                    Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
                }
                wy = order / 2.0 * Fourier_spreading_appx;
                wy = wy * wy;
              }

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);

                double ph_2_kz_c = order * (zprd/nz_pppm) / 2.0 * fabs(qz) / spreading_select_c;
                wz = 0.00;
                if(ph_2_kz_c <= 1.00){
                  ph_2_kz_c = 2.0 * ph_2_kz_c - 1.0;
                  double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
                  double Fourier_spreading_r = 1.0;
                  for(int i=1; i<fourier_spreading_order; i++){
                    Fourier_spreading_r *= ph_2_kz_c;
                    Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
                  }
                  wz = order / 2.0 * Fourier_spreading_appx;
                  wz = wz * wz;
                }

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                double k_rc_c = sqrt(qx*qx+qy*qy+qz*qz) * cutoff / select_c;
                dot2 = 0.00;
                if(k_rc_c <= 1.00){
                  k_rc_c = 2.0 * k_rc_c - 1.0;
                  double Fourier_poly_appx = fourier_split_poly_coeff[0];
                  double Fourier_poly_r = 1.0;
                  for(int i=1; i<num_of_Fourier_poly; i++){
                    Fourier_poly_r *= k_rc_c;
                    Fourier_poly_appx += fourier_split_poly_coeff[i] * Fourier_poly_r;
                  }
                  dot2 = 2 * 3.14159265 * Fourier_poly_appx / (qx*qx+qy*qy+qz*qz);
                }
                else
                  dot2 = 0.0;

                sum1 += dot1 * dot2 * wx * wy * wz;
              }
            }
          }

          if(denominator==0.00){
             greensfn[n++] = 0.00;
          }
          else{
             greensfn[n++] = numerator*sum1/denominator;
          }

        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
   for a triclinic system
------------------------------------------------------------------------- */

void ESP::compute_gf_ik_triclinic()
{
  double wx,wy,wz,qx,qy,qz;
  double sum1,dot1,dot2;
  double numerator,denominator;
  double sqk;

  int k,l,m,n,nx,ny,nz,kper,lper,mper;

  const double * const prd = domain->prd;
  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd * slab_volfactor;

  double tmp[3];
  tmp[0] = static_cast<int> (select_c / (MY_2PI * cutoff * nx_pppm)) * pow(-log(EPS_HOC),0.25);
  tmp[1] = static_cast<int> (select_c / (MY_2PI * cutoff * ny_pppm)) * pow(-log(EPS_HOC),0.25);
  tmp[2] = static_cast<int> (select_c / (MY_2PI * cutoff * nz_pppm)) * pow(-log(EPS_HOC),0.25);

  lamda2xT(&tmp[0], &tmp[0]);
  const int nbx = static_cast<int> (tmp[0]);
  const int nby = static_cast<int> (tmp[1]);
  const int nbz = static_cast<int> (tmp[2]);

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        double unitk_lamda[3];
        unitk_lamda[0] = 2.0*MY_PI*kper;
        unitk_lamda[1] = 2.0*MY_PI*lper;
        unitk_lamda[2] = 2.0*MY_PI*mper;
        x2lamdaT(&unitk_lamda[0],&unitk_lamda[0]);

        sqk = square(unitk_lamda[0]) + square(unitk_lamda[1]) + square(unitk_lamda[2]);

        if (sqk != 0.0) {
          numerator = 1.0 / sqk;
          denominator = gf_denom_psw(unitk_lamda[0], unitk_lamda[1], unitk_lamda[2],
                                     xprd / nx_pppm, yprd / ny_pppm, zprd_slab / nz_pppm);
          sum1 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            //qx = unitk_lamda[0] + 2.0 * MY_PI * nx_pppm * nx;

            double ph_2_kx_c = order * fabs(MY_PI*kper/nx_pppm + MY_PI*nx) / spreading_select_c;
            wx = 0.00;
            if (ph_2_kx_c <= 1.00) {
              ph_2_kx_c = 2.0 * ph_2_kx_c - 1.0;
              double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
              double Fourier_spreading_r = 1.0;
              for (int i = 1; i < fourier_spreading_order; i++) {
                Fourier_spreading_r *= ph_2_kx_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
              }
              wx = order / 2.0 * Fourier_spreading_appx;
              wx = wx * wx;
            }

            for (ny = -nby; ny <= nby; ny++) {
              //qy = unitk_lamda[1] + 2.0 * MY_PI * ny_pppm * ny;

              double ph_2_ky_c = order * fabs(MY_PI*lper/ny_pppm + MY_PI*ny) / spreading_select_c;
              wy = 0.00;
              if (ph_2_ky_c <= 1.00) {
                ph_2_ky_c = 2.0 * ph_2_ky_c - 1.0;
                double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
                double Fourier_spreading_r = 1.0;
                for (int i = 1; i < fourier_spreading_order; i++) {
                  Fourier_spreading_r *= ph_2_ky_c;
                  Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
                }
                wy = order / 2.0 * Fourier_spreading_appx;
                wy = wy * wy;
              }

              for (nz = -nbz; nz <= nbz; nz++) {
                //qz = unitk_lamda[2] + 2.0 * MY_PI * nz_pppm * nz;

                double ph_2_kz_c = order * fabs(MY_PI*mper/nz_pppm + MY_PI*nz) / spreading_select_c;
                wz = 0.00;
                if (ph_2_kz_c <= 1.00) {
                  ph_2_kz_c = 2.0 * ph_2_kz_c - 1.0;
                  double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
                  double Fourier_spreading_r = 1.0;
                  for (int i = 1; i < fourier_spreading_order; i++) {
                    Fourier_spreading_r *= ph_2_kz_c;
                    Fourier_spreading_appx += fourier_spread_poly_coeff[i] * Fourier_spreading_r;
                  }
                  wz = order / 2.0 * Fourier_spreading_appx;
                  wz = wz * wz;
                }

                double b[3];
                b[0] = 2.0*MY_PI*nx_pppm*nx;
                b[1] = 2.0*MY_PI*ny_pppm*ny;
                b[2] = 2.0*MY_PI*nz_pppm*nz;
                x2lamdaT(&b[0],&b[0]);

                qx = unitk_lamda[0]+b[0];
                qy = unitk_lamda[1]+b[1];
                qz = unitk_lamda[2]+b[2];

                dot1 = unitk_lamda[0]*qx + unitk_lamda[1]*qy + unitk_lamda[2]*qz;
                double k_rc_c = sqrt(qx * qx + qy * qy + qz * qz) * cutoff / select_c;
                dot2 = 0.00;
                if (k_rc_c <= 1.00) {
                  k_rc_c = 2.0 * k_rc_c - 1.0;
                  double Fourier_poly_appx = fourier_split_poly_coeff[0];
                  double Fourier_poly_r = 1.0;
                  for (int i = 1; i < num_of_Fourier_poly; i++) {
                    Fourier_poly_r *= k_rc_c;
                    Fourier_poly_appx += fourier_split_poly_coeff[i] * Fourier_poly_r;
                  }
                  dot2 = 2 * MY_PI * Fourier_poly_appx / (qx * qx + qy * qy + qz * qz);
                } else
                  dot2 = 0.0;

                sum1 += dot1 * dot2 * wx * wy * wz;
              }
            }
          }

          if(denominator==0.00){
             greensfn[n++] = 0.00;
          }
          else{
             greensfn[n++] = numerator*sum1/denominator;
          }

        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
------------------------------------------------------------------------- */

void ESP::compute_gf_ad()
{
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];

  const double zprd_slab = zprd * slab_volfactor;

  const double unitkx = (MY_2PI / xprd);
  const double unitky = (MY_2PI / yprd);
  const double unitkz = (MY_2PI / zprd_slab);

  const double sqk_cut2 = (select_c / cutoff) * (select_c / cutoff);

  const double scale_x = (order * MY_PI) / ((double)nx_pppm * spreading_select_c);
  const double scale_y = (order * MY_PI) / ((double)ny_pppm * spreading_select_c);
  const double scale_z = (order * MY_PI) / ((double)nz_pppm * spreading_select_c * slab_volfactor);

  for (int i = 0; i < 6; i++) sf_coeff[i] = 0.0;

  const int nx = nxhi_fft - nxlo_fft + 1;
  const int ny = nyhi_fft - nylo_fft + 1;
  const int nz = nzhi_fft - nzlo_fft + 1;

  // --- precompute per-dimension tables
  std::vector<int>    kper(nx), lper(ny), mper(nz);
  std::vector<double> qx(nx), qy(ny), qz(nz);
  std::vector<double> qx2(nx), qy2(ny), qz2(nz);

  std::vector<double> wx_denom(nx), wy_denom(ny), wz_denom(nz);

  // k dimension
  for (int ii = 0; ii < nx; ++ii) {
    const int k = nxlo_fft + ii;
    const int kp = k - nx_pppm * (2 * k / nx_pppm);
    kper[ii] = kp;

    const double q = unitkx * (double)kp;
    qx[ii]  = q;
    qx2[ii] = q * q;

    wx_denom[ii] = spreading_weight2_from_abs_index(std::abs(kp), scale_x);
  }

  // l dimension
  for (int ll = 0; ll < ny; ++ll) {
    const int l = nylo_fft + ll;
    const int lp = l - ny_pppm * (2 * l / ny_pppm);
    lper[ll] = lp;

    const double q = unitky * (double)lp;
    qy[ll]  = q;
    qy2[ll] = q * q;

    wy_denom[ll] = spreading_weight2_from_abs_index(std::abs(lp), scale_y);
  }

  // m dimension
  for (int mm = 0; mm < nz; ++mm) {
    const int m = nzlo_fft + mm;
    const int mp = m - nz_pppm * (2 * m / nz_pppm);
    mper[mm] = mp;

    const double q = unitkz * (double)mp; // slab kz
    qz[mm]  = q;
    qz2[mm] = q * q;

    wz_denom[mm] = spreading_weight2_from_abs_index(std::abs(mp), scale_z);
  }

  // local accumulators
  double c0=0.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0, c5=0.0;

  // pointers
  double *g1 = greensfn;
  double *g2 = greensfn2;

  const double *p1 = sf_precoeff1;
  const double *p2 = sf_precoeff2;
  const double *p3 = sf_precoeff3;
  const double *p4 = sf_precoeff4;
  const double *p5 = sf_precoeff5;
  const double *p6 = sf_precoeff6;

  int n = 0;
  for (int mm = 0; mm < nz; ++mm) {
    const double wz = wz_denom[mm];

    for (int ll = 0; ll < ny; ++ll) {
      const double wy = wy_denom[ll];

      const double wyz = wy * wz;

      for (int ii = 0; ii < nx; ++ii, ++n) {
        const double wx = wx_denom[ii];

        // base weight product
        const double wxyz = wx * wyz;

        // sqk
        const double sqk = qx2[ii] + qy2[ll] + qz2[mm];
        g1[n] = 0.0;
        g2[n] = 0.0;
        if (sqk == 0.0) continue;

        // denom = (sumx*sumy*sumz)^2  (sumx=denomx[ii], etc.)
        const double denom = wxyz;
        if (denom == 0.0) continue;

        // outside cutoff => dot2=0, greensfn stays 0
        if (sqk > sqk_cut2) continue;

        // k_rc_c_old in [0,1], then map to x in [-1,1]
        const double krc_old = std::sqrt(sqk) * cutoff / select_c;
        const double x = 2.0 * krc_old - 1.0;

        double poly, dpoly_dx;
        poly_and_deriv_horner(x, fourier_split_poly_coeff, num_of_Fourier_poly, poly, dpoly_dx);

        const double poly_virial = 2.0 * krc_old * dpoly_dx;

        const double invden = 1.0 / denom;

        g1[n] = (MY_2PI * poly) * invden / sqk;

        const double sqk2 = sqk * sqk;
        g2[n] = (MY_2PI * poly_virial) * invden / sqk2;

        const double gg = g1[n];
        c0 += p1[n] * gg;  c1 += p2[n] * gg;
        c2 += p3[n] * gg;  c3 += p4[n] * gg;
        c4 += p5[n] * gg;  c5 += p6[n] * gg;
      }
    }
  }

  sf_coeff[0] = c0; sf_coeff[1] = c1; sf_coeff[2] = c2;
  sf_coeff[3] = c3; sf_coeff[4] = c4; sf_coeff[5] = c5;

  // self-force correction
  double prex = MY_PI / volume;
  double prey = MY_PI / volume;
  double prez = MY_PI / volume;

  prex *= nx_pppm / xprd;
  prey *= ny_pppm / yprd;
  prez *= nz_pppm / zprd_slab;

  sf_coeff[0] *= prex;
  sf_coeff[1] *= prex * 2.0;
  sf_coeff[2] *= prey;
  sf_coeff[3] *= prey * 2.0;
  sf_coeff[4] *= prez;
  sf_coeff[5] *= prez * 2.0;

  // MPI reduce
  double tmp[6];
  MPI_Allreduce(sf_coeff, tmp, 6, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; i++) sf_coeff[i] = tmp[i];
}

/* ----------------------------------------------------------------------
   compute self force coefficients for ad-differentiation scheme
------------------------------------------------------------------------- */

void ESP::compute_sf_precoeff()
{
  int i,k,l,m,n;
  int nx,ny,nz,kper,lper,mper;
  double wx0[5],wy0[5],wz0[5],wx1[5],wy1[5],wz1[5],wx2[5],wy2[5],wz2[5];
  double qx0,qy0,qz0,qx1,qy1,qz1,qx2,qy2,qz2;
  double u0,u1,u2,u3,u4,u5,u6;
  double sum1,sum2,sum3,sum4,sum5,sum6;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = 0.0;
        for (i = 0; i < 5; i++) {

          qx0 = MY_2PI*(kper+nx_pppm*(i-2));
          qx1 = MY_2PI*(kper+nx_pppm*(i-1));
          qx2 = MY_2PI*(kper+nx_pppm*(i  ));

          double ph_2_kx_c = (0.5 * order / nx_pppm) * (std::abs(qx0) / spreading_select_c);
          wx0[i] = 0.00;
          if(ph_2_kx_c <= 1.00){
            ph_2_kx_c = 2.0 * ph_2_kx_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kx_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wx0[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_kx_c = (0.5 * order / nx_pppm) * (std::abs(qx1) / spreading_select_c);
          wx1[i] = 0.00;
          if(ph_2_kx_c <= 1.00){
            ph_2_kx_c = 2.0 * ph_2_kx_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kx_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wx1[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_kx_c = (0.5 * order / nx_pppm) * (std::abs(qx2) / spreading_select_c);
          wx2[i] = 0.00;
          if(ph_2_kx_c <= 1.00){
            ph_2_kx_c = 2.0 * ph_2_kx_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kx_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wx2[i] = order / 2.0 * Fourier_spreading_appx;
          }

          qy0 = MY_2PI*(lper+ny_pppm*(i-2));
          qy1 = MY_2PI*(lper+ny_pppm*(i-1));
          qy2 = MY_2PI*(lper+ny_pppm*(i  ));

          double ph_2_ky_c = (0.5 * order / ny_pppm) * (std::abs(qy0) / spreading_select_c);
          wy0[i] = 0.00;
          if(ph_2_ky_c <= 1.00){
            ph_2_ky_c = 2.0 * ph_2_ky_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_ky_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wy0[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_ky_c = (0.5 * order / ny_pppm) * (std::abs(qy1) / spreading_select_c);
          wy1[i] = 0.00;
          if(ph_2_ky_c <= 1.00){
            ph_2_ky_c = 2.0 * ph_2_ky_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_ky_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wy1[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_ky_c = (0.5 * order / ny_pppm) * (std::abs(qy2) / spreading_select_c);
          wy2[i] = 0.00;
          if(ph_2_ky_c <= 1.00){
            ph_2_ky_c = 2.0 * ph_2_ky_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_ky_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wy2[i] = order / 2.0 * Fourier_spreading_appx;
          }

          qz0 = MY_2PI*(mper+nz_pppm*(i-2));
          qz1 = MY_2PI*(mper+nz_pppm*(i-1));
          qz2 = MY_2PI*(mper+nz_pppm*(i  ));

          double ph_2_kz_c = (0.5 * order / nz_pppm) * (std::abs(qz0) / spreading_select_c);
          wz0[i] = 0.00;
          if(ph_2_kz_c <= 1.00){
            ph_2_kz_c = 2.0 * ph_2_kz_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kz_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wz0[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_kz_c = (0.5 * order / nz_pppm) * (std::abs(qz1) / spreading_select_c);

          wz1[i] = 0.00;
          if(ph_2_kz_c <= 1.00){
            ph_2_kz_c = 2.0 * ph_2_kz_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kz_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wz1[i] = order / 2.0 * Fourier_spreading_appx;
          }

          ph_2_kz_c = (0.5 * order / nz_pppm) * (std::abs(qz2) / spreading_select_c);
          wz2[i] = 0.00;
          if(ph_2_kz_c <= 1.00){
            ph_2_kz_c = 2.0 * ph_2_kz_c - 1.0;
            double Fourier_spreading_appx = fourier_spread_poly_coeff[0];
            double Fourier_spreading_r = 1.0;
            for(int ii=1; ii<fourier_spreading_order; ii++){
                Fourier_spreading_r *= ph_2_kz_c;
                Fourier_spreading_appx += fourier_spread_poly_coeff[ii] * Fourier_spreading_r;
            }
            wz2[i] = order / 2.0 * Fourier_spreading_appx;
          }
        }

        for (nx = 0; nx < 5; nx++) {
          for (ny = 0; ny < 5; ny++) {
            for (nz = 0; nz < 5; nz++) {
              u0 = wx0[nx]*wy0[ny]*wz0[nz];
              u1 = wx1[nx]*wy0[ny]*wz0[nz];
              u2 = wx2[nx]*wy0[ny]*wz0[nz];
              u3 = wx0[nx]*wy1[ny]*wz0[nz];
              u4 = wx0[nx]*wy2[ny]*wz0[nz];
              u5 = wx0[nx]*wy0[ny]*wz1[nz];
              u6 = wx0[nx]*wy0[ny]*wz2[nz];

              sum1 += u0*u1;
              sum2 += u0*u2;
              sum3 += u0*u3;
              sum4 += u0*u4;
              sum5 += u0*u5;
              sum6 += u0*u6;
            }
          }
        }

        sf_precoeff1[n] = sum1;
        sf_precoeff2[n] = sum2;
        sf_precoeff3[n] = sum3;
        sf_precoeff4[n] = sum4;
        sf_precoeff5[n] = sum5;
        sf_precoeff6[n++] = sum6;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void ESP::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  for (int i = 0; i < nlocal; i++) {

    // order = even:
    //   (nx,ny,nz) = global index of grid pt to "lower left" of charge
    // order = odd:
    //   (nx,ny,nz) = global index of grid pt closest to charge due to shift
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute ESP");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void ESP::make_rho()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global indices of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void ESP::brick2fft()
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_fft[n++] = density_brick[iz][iy][ix];

  remap->perform(density_fft,density_fft,work1);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver
------------------------------------------------------------------------- */

void ESP::poisson()
{
  if (differentiation_flag == 1) poisson_ad();
  else poisson_ik();
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik
------------------------------------------------------------------------- */

void ESP::poisson_ik()
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }

  fft1->compute(work1,work1,FFT3d::FORWARD);

  // global energy and virial contribution

  bigint ngridtotal = (bigint) nx_pppm * ny_pppm * nz_pppm;
  double scaleinv = 1.0/ngridtotal;
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        for (j = 0; j < 6; j++) virial[j] += eng*vg[i][j];
        if (eflag_global) energy += eng;
        n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
        energy +=
          s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

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

  // x direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = -fkx[i]*work1[n+1];
        work2[n+1] = fkx[i]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdx_brick[k][j][i] = work2[n];
        n += 2;
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = -fky[j]*work1[n+1];
        work2[n+1] = fky[j]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdy_brick[k][j][i] = work2[n];
        n += 2;
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = -fkz[k]*work1[n+1];
        work2[n+1] = fkz[k]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdz_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik for a triclinic system
------------------------------------------------------------------------- */

void ESP::poisson_ik_triclinic()
{
  int i,j,k,n;

  // compute gradients of V(r) in each of 3 dims by transforming ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = -fkx[i]*work1[n+1];
    work2[n+1] = fkx[i]*work1[n];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdx_brick[k][j][i] = work2[n];
        n += 2;
      }

  // y direction gradient

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = -fky[i]*work1[n+1];
    work2[n+1] = fky[i]*work1[n];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdy_brick[k][j][i] = work2[n];
        n += 2;
      }

  // z direction gradient

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = -fkz[i]*work1[n+1];
    work2[n+1] = fkz[i]*work1[n];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdz_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ad
------------------------------------------------------------------------- */

void ESP::poisson_ad()
{
  int i,j,k,n;
  double eng,eng_virial;

  // transform charge density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }

  fft1->compute(work1,work1,FFT3d::FORWARD);

  // global energy and virial contribution

  bigint ngridtotal = (bigint) nx_pppm * ny_pppm * nz_pppm;
  double scaleinv = 1.0/ngridtotal;
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        eng_virial = s2 * greensfn2[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        for (j = 0; j < 6; j++)
        {
          virial[j] += eng*vg[i][j]; // The first term
          virial[j] += eng_virial*vg2[i][j]; // The second term
        }
        if (eflag_global) energy += eng;
        n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
        energy +=
          s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

  // extra FFTs for per-atom energy/virial

  if (vflag_atom) poisson_peratom();

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n];
    work2[n+1] = work1[n+1];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        u_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for per-atom energy/virial
------------------------------------------------------------------------- */

void ESP::poisson_peratom()
{
  int i,j,k,n;

  // energy

  if (eflag_atom && differentiation_flag != 1) {
    n = 0;
    for (i = 0; i < nfft; i++) {
      work2[n] = work1[n];
      work2[n+1] = work1[n+1];
      n += 2;
    }

    fft2->compute(work2,work2,FFT3d::BACKWARD);

    n = 0;
    for (k = nzlo_in; k <= nzhi_in; k++)
      for (j = nylo_in; j <= nyhi_in; j++)
        for (i = nxlo_in; i <= nxhi_in; i++) {
          u_brick[k][j][i] = work2[n];
          n += 2;
        }
  }

  // 6 components of virial in v0 thru v5

  if (!vflag_atom) return;

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][0];
    work2[n+1] = work1[n+1]*vg[i][0];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v0_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][1];
    work2[n+1] = work1[n+1]*vg[i][1];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v1_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][2];
    work2[n+1] = work1[n+1]*vg[i][2];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v2_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][3];
    work2[n+1] = work1[n+1]*vg[i][3];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v3_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][4];
    work2[n+1] = work1[n+1]*vg[i][4];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v4_brick[k][j][i] = work2[n];
        n += 2;
      }

  n = 0;
  for (i = 0; i < nfft; i++) {
    work2[n] = work1[n]*vg[i][5];
    work2[n+1] = work1[n+1]*vg[i][5];
    n += 2;
  }

  fft2->compute(work2,work2,FFT3d::BACKWARD);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        v5_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

void ESP::fieldforce()
{
  if (differentiation_flag == 1) fieldforce_ad();
  else fieldforce_ik();
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void ESP::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const double qfactor = qqrd2e * scale * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void ESP::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz;
  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  prd = domain->prd;
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);
    compute_drho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }

    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force and subtract self forces

    const double qfactor = qqrd2e * scale;

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*q[i]*q[i];

    f[i][0] += qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*q[i]*q[i];
    f[i][1] += qfactor*(eky*q[i] - sf);

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*q[i]*q[i];

    if (slabflag != 2)
    {
      f[i][2] += qfactor*(ekz*q[i] - sf);
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void ESP::fieldforce_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    u = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) u += x0*u_brick[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0_brick[mz][my][mx];
            v1 += x0*v1_brick[mz][my][mx];
            v2 += x0*v2_brick[mz][my][mx];
            v3 += x0*v3_brick[mz][my][mx];
            v4 += x0*v4_brick[mz][my][mx];
            v5 += x0*v5_brick[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) eatom[i] += q[i]*u;
    if (vflag_atom) {
      vatom[i][0] += q[i]*v0;
      vatom[i][1] += q[i]*v1;
      vatom[i][2] += q[i]*v2;
      vatom[i][3] += q[i]*v3;
      vatom[i][4] += q[i]*v4;
      vatom[i][5] += q[i]*v5;
    }
  }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void ESP::pack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  auto *buf = (FFT_SCALAR *) vbuf;

  int n = 0;

  if (flag == FORWARD_IK) {
    FFT_SCALAR *xsrc = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ysrc = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zsrc = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = xsrc[list[i]];
      buf[n++] = ysrc[list[i]];
      buf[n++] = zsrc[list[i]];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *src = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) buf[n++] = esrc[list[i]];
      if (vflag_atom) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src[list[i]];
      buf[n++] = v1src[list[i]];
      buf[n++] = v2src[list[i]];
      buf[n++] = v3src[list[i]];
      buf[n++] = v4src[list[i]];
      buf[n++] = v5src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void ESP::unpack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  auto *buf = (FFT_SCALAR *) vbuf;

  int n = 0;

  if (flag == FORWARD_IK) {
    FFT_SCALAR *xdest = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ydest = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zdest = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      xdest[list[i]] = buf[n++];
      ydest[list[i]] = buf[n++];
      zdest[list[i]] = buf[n++];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *dest = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) esrc[list[i]] = buf[n++];
      if (vflag_atom) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[n++];
      v1src[list[i]] = buf[n++];
      v2src[list[i]] = buf[n++];
      v3src[list[i]] = buf[n++];
      v4src[list[i]] = buf[n++];
      v5src[list[i]] = buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void ESP::pack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  auto *buf = (FFT_SCALAR *) vbuf;

  if (flag == REVERSE_RHO) {
    FFT_SCALAR *src = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void ESP::unpack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  auto *buf = (FFT_SCALAR *) vbuf;

  if (flag == REVERSE_RHO) {
    FFT_SCALAR *dest = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
------------------------------------------------------------------------- */

void ESP::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
{
  // loop thru all possible factorizations of nprocs
  // surf = surface area of largest proc sub-domain
  // innermost if test minimizes surface area and surface/volume ratio

  int bestsurf = 2 * (nx + ny);
  int bestboxx = 0;
  int bestboxy = 0;

  int boxx,boxy,surf,ipx,ipy;

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = nx/ipx;
      if (nx % ipx) boxx++;
      boxy = ny/ipy;
      if (ny % ipy) boxy++;
      surf = boxx + boxy;
      if (surf < bestsurf ||
          (surf == bestsurf && boxx*boxy > bestboxx*bestboxy)) {
        bestsurf = surf;
        bestboxx = boxx;
        bestboxy = boxy;
        *px = ipx;
        *py = ipy;
      }
    }
    ipx++;
  }
}

/* ----------------------------------------------------------------------
   estimate the spreading order of the ESP method
------------------------------------------------------------------------- */

int ESP::estimate_order(double accuracy)
{
    // p = -log10(acc)
    double p = -std::log10(accuracy);

    double p_round = std::round(p);
    const double eps = 0.2;

    int order;
    if (std::abs(p - p_round) < eps) {
        order = static_cast<int>(p_round) * 2 - 2;
    } else {
        order = static_cast<int>(std::ceil(p)) * 2 - 3;
    }

    if (order < 2) order = 2;
    if (order > 16) order = 16;
    return order;
}

/* —----------------------------------------------------------------------
   build table for ESP method
------------------------------------------------------------------------- */

void ESP::build_table(double spreading_accuracy)
{
  // force kernel
  std::vector<double> poly_coeff;
  MathPSWF::force_poly(0.1*accuracy_relative, select_c, poly_coeff);
  num_of_force_poly = poly_coeff.size();
  memory->create(force_poly_coeff, num_of_force_poly, "ESP:force_poly_coeff");
  for (int i = 0; i < num_of_force_poly; i++) force_poly_coeff[i] = poly_coeff[i];
  if (me == 0) utils::logmesg(lmp,"  Force poly size: {}\n",num_of_force_poly);

  // energy kernel
  poly_coeff.clear();
  MathPSWF::energy_poly(0.01*accuracy_relative, select_c, poly_coeff);
  num_of_energy_poly = poly_coeff.size();
  memory->create(energy_poly_coeff, num_of_energy_poly, "ESP:energy_poly_coeff");
  for (int i = 0; i < num_of_energy_poly; i++) energy_poly_coeff[i] = poly_coeff[i];
  if (me == 0) utils::logmesg(lmp,"  Energy poly size: {}\n",num_of_energy_poly);

  // Fourier space kernel
  poly_coeff.clear();
  MathPSWF::fourier_poly(0.1*accuracy_relative, select_c, Lambda_0, poly_coeff);
  num_of_Fourier_poly = poly_coeff.size();
  memory->create(fourier_split_poly_coeff, num_of_Fourier_poly, "ESP:fourier_split_poly_coeff");
  for (int i = 0; i < num_of_Fourier_poly; i++) fourier_split_poly_coeff[i] = poly_coeff[i];
  if (me == 0) utils::logmesg(lmp,"  Fourier poly size: {}\n",num_of_Fourier_poly);

  // spreading kernel in real space - need to be consistent with the spreading accuracy
  poly_coeff.clear();
  MathPSWF::spread_real_poly(order, 0.1*spreading_accuracy, spreading_select_c, poly_coeff);
  poly_order = poly_coeff.size() / order;

  memory->create2d_offset(rho_coeff,poly_order,(1-order)/2,order/2,"esp:rho_coeff");
  memory->create2d_offset(drho_coeff,poly_order,(1-order)/2,order/2,"esp:drho_coeff");

  for(int i=0; i<poly_order; i++){
    for(int j=0; j<order; j++){
        rho_coeff[i][j+(1-order)/2] = poly_coeff[i*order+j];
    }
  }
  for (int m = -(order-1)/2; m <= order/2; m += 1) {
    for (int l = 1; l < poly_order; l++)
        drho_coeff[l-1][m] = l*rho_coeff[l][m]; // Coefficients for l x^l-1 terms
    drho_coeff[poly_order-1][m] = 0.00;
  }
  if (me == 0) utils::logmesg(lmp,"  Spread real poly size: {}\n",poly_order);

  // spreading kernel in Fourier space
  poly_coeff.clear();
  MathPSWF::spread_fourier_poly(0.1*spreading_accuracy, spreading_select_c, spreading_Lambda_0, poly_coeff);
  fourier_spreading_order = poly_coeff.size();
  memory->create(fourier_spread_poly_coeff, fourier_spreading_order, "ESP:fourier_spread_poly_coeff");
  for (int i = 0; i < fourier_spreading_order; i++) fourier_spread_poly_coeff[i] = poly_coeff[i];
  if (me == 0) utils::logmesg(lmp,"  Fourier spreading poly size: {}\n",fourier_spreading_order);
}


/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void ESP::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                         const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  // k: order of spreading points
  // l: order of polynomials
  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = poly_order-1; l >=0; l--) {
      r1 = rho_coeff[l][k] + r1*dx;
      r2 = rho_coeff[l][k] + r2*dy;
      r3 = rho_coeff[l][k] + r3*dz;
    }
    rho1d[0][k] = r1;
    rho1d[1][k] = r2;
    rho1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into drho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void ESP::compute_drho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                          const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = poly_order-2; l >= 0; l--) {
      r1 = drho_coeff[l][k] + r1*dx;
      r2 = drho_coeff[l][k] + r2*dy;
      r3 = drho_coeff[l][k] + r3*dz;
    }
    drho1d[0][k] = r1;
    drho1d[1][k] = r2;
    drho1d[2][k] = r3;
  }
}

void ESP::compute_rho_coeff()
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;
  memory->create2d_offset(a,order,-order,order,"esp:a");

  for (k = -order; k <= order; k++)
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;

  a[0][0] = 1.0;
  for (j = 1; j < order; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
        a[l+1][k] = (a[l][k+1]-a[l][k-1]) / (l+1);
#ifdef FFT_SINGLE
        s += powf(0.5,(float) l+1) *
          (a[l][k-1] + powf(-1.0,(float) l) * a[l][k+1]) / (l+1);
#else
        s += pow(0.5,(double) l+1) *
          (a[l][k-1] + pow(-1.0,(double) l) * a[l][k+1]) / (l+1);
#endif
      }
      a[0][k] = s;
    }
  }

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      rho_coeff[l][m] = a[l][k]; // Coefficients for x^l terms
    for (l = 1; l < order; l++)
      drho_coeff[l-1][m] = l*a[l][k]; // Coefficients for l x^l-1 terms
    m++;
  }

  memory->destroy2d_offset(a,-order);
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void ESP::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd_slab = domain->zprd*slab_volfactor;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd_slab*zprd_slab/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += ffact * q[i]*(dipole_all - qsum*x[i][2]);
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

int ESP::timing_1d(int n, double &time1d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;

  MPI_Barrier(world);
  time1 = platform::walltime();

  for (int i = 0; i < n; i++) {
    fft1->timing1d(work1,nfft_both,FFT3d::FORWARD);
    fft2->timing1d(work1,nfft_both,FFT3d::BACKWARD);
    if (differentiation_flag != 1) {
      fft2->timing1d(work1,nfft_both,FFT3d::BACKWARD);
      fft2->timing1d(work1,nfft_both,FFT3d::BACKWARD);
    }
  }

  MPI_Barrier(world);
  time2 = platform::walltime();
  time1d = time2 - time1;

  if (differentiation_flag) return 2;
  return 4;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

int ESP::timing_3d(int n, double &time3d)
{
  double time1,time2;

  for (int i = 0; i < 2*nfft_both; i++) work1[i] = ZEROF;

  MPI_Barrier(world);
  time1 = platform::walltime();

  for (int i = 0; i < n; i++) {
    fft1->compute(work1,work1,FFT3d::FORWARD);
    fft2->compute(work1,work1,FFT3d::BACKWARD);
    if (differentiation_flag != 1) {
      fft2->compute(work1,work1,FFT3d::BACKWARD);
      fft2->compute(work1,work1,FFT3d::BACKWARD);
    }
  }

  MPI_Barrier(world);
  time2 = platform::walltime();
  time3d = time2 - time1;

  if (differentiation_flag) return 2;
  return 4;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double ESP::memory_usage()
{
  double bytes = (double)nmax*3 * sizeof(double);

  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  if (differentiation_flag == 1) {
    bytes += (double)2 * nbrick * sizeof(FFT_SCALAR);
  } else {
    bytes += (double)4 * nbrick * sizeof(FFT_SCALAR);
  }

  if (triclinic) bytes += (double)3 * nfft_both * sizeof(double);
  bytes += (double)6 * nfft_both * sizeof(double);
  bytes += (double)nfft_both * sizeof(double);
  bytes += (double)nfft_both*5 * sizeof(FFT_SCALAR);

  if (peratom_allocate_flag)
    bytes += (double)6 * nbrick * sizeof(FFT_SCALAR);

  if (group_allocate_flag) {
    bytes += (double)2 * nbrick * sizeof(FFT_SCALAR);
    bytes += (double)2 * nfft_both * sizeof(FFT_SCALAR);
  }

  // two Grid3d bufs

  bytes += (double)(ngc_buf1 + ngc_buf2) * npergrid * sizeof(FFT_SCALAR);

  return bytes;
}

/* ----------------------------------------------------------------------
   group-group interactions
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the ESP total long-range force and energy for groups A and B
 ------------------------------------------------------------------------- */

void ESP::compute_group_group(int groupbit_A, int groupbit_B, int AA_flag)
{
  if (slabflag && triclinic)
    error->all(FLERR,"Cannot (yet) use K-space slab "
               "correction with compute group/group for triclinic systems");

  if (differentiation_flag)
    error->all(FLERR,"Cannot (yet) use kspace_modify "
               "diff ad with compute group/group");

  if (!group_allocate_flag) allocate_groups();

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  e2group = 0.0; //energy
  f2group[0] = 0.0; //force in x-direction
  f2group[1] = 0.0; //force in y-direction
  f2group[2] = 0.0; //force in z-direction

  // map my particle charge onto my local 3d density grid

  make_rho_groups(groupbit_A,groupbit_B,AA_flag);

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  // temporarily store and switch pointers so we can
  //  use brick2fft() for groups A and B (without
  //  writing an additional function)

  FFT_SCALAR ***density_brick_real = density_brick;
  FFT_SCALAR *density_fft_real = density_fft;

  // group A

  density_brick = density_A_brick;
  density_fft = density_A_fft;

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // group B

  density_brick = density_B_brick;
  density_fft = density_B_fft;

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // switch back pointers

  density_brick = density_brick_real;
  density_fft = density_fft_real;

  // compute potential gradient on my FFT grid and
  //   portion of group-group energy/force on this proc's FFT grid

  poisson_groups(AA_flag);

  const double qscale = qqrd2e * scale;

  // total group A <--> group B energy
  // self and boundary correction terms are in compute_group_group.cpp

  double e2group_all;
  MPI_Allreduce(&e2group,&e2group_all,1,MPI_DOUBLE,MPI_SUM,world);
  e2group = e2group_all;

  e2group *= qscale*0.5*volume;

  // total group A <--> group B force

  double f2group_all[3];
  MPI_Allreduce(f2group,f2group_all,3,MPI_DOUBLE,MPI_SUM,world);

  f2group[0] = qscale*volume*f2group_all[0];
  f2group[1] = qscale*volume*f2group_all[1];
  if (slabflag != 2) f2group[2] = qscale*volume*f2group_all[2];

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);

  if (slabflag == 1)
    slabcorr_groups(groupbit_A, groupbit_B, AA_flag);
}

/* ----------------------------------------------------------------------
 allocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */

void ESP::allocate_groups()
{
  group_allocate_flag = 1;

  memory->create3d_offset(density_A_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:density_A_brick");
  memory->create3d_offset(density_B_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"esp:density_B_brick");
  memory->create(density_A_fft,nfft_both,"esp:density_A_fft");
  memory->create(density_B_fft,nfft_both,"esp:density_B_fft");
}

/* ----------------------------------------------------------------------
 deallocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */

void ESP::deallocate_groups()
{
  group_allocate_flag = 0;

  memory->destroy3d_offset(density_A_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(density_B_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(density_A_fft);
  memory->destroy(density_B_fft);
}

/* ----------------------------------------------------------------------
 create discretized "density" on section of global grid due to my particles
 density(x,y,z) = charge "density" at grid points of my 3d brick
 (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
 in global grid for group-group interactions
 ------------------------------------------------------------------------- */

void ESP::make_rho_groups(int groupbit_A, int groupbit_B, int AA_flag)
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density arrays

  memset(&(density_A_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  memset(&(density_B_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  for (int i = 0; i < nlocal; i++) {

    if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
      if (AA_flag) continue;

    if ((mask[i] & groupbit_A) || (mask[i] & groupbit_B)) {

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

      compute_rho1d(dx,dy,dz);

      z0 = delvolinv * q[i];
      for (n = nlower; n <= nupper; n++) {
        mz = n+nz;
        y0 = z0*rho1d[2][n];
        for (m = nlower; m <= nupper; m++) {
          my = m+ny;
          x0 = y0*rho1d[1][m];
          for (l = nlower; l <= nupper; l++) {
            mx = l+nx;

            // group A

            if (mask[i] & groupbit_A)
              density_A_brick[mz][my][mx] += x0*rho1d[0][l];

            // group B

            if (mask[i] & groupbit_B)
              density_B_brick[mz][my][mx] += x0*rho1d[0][l];
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for group-group interactions
 ------------------------------------------------------------------------- */

void ESP::poisson_groups(int AA_flag)
{
  int i,j,k,n;

  // reuse memory (already declared)

  FFT_SCALAR *work_A = work1;
  FFT_SCALAR *work_B = work2;

  // transform charge density (r -> k)

  // group A

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_A[n++] = density_A_fft[i];
    work_A[n++] = ZEROF;
  }

  fft1->compute(work_A,work_A,FFT3d::FORWARD);

  // group B

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_B[n++] = density_B_fft[i];
    work_B[n++] = ZEROF;
  }

  fft1->compute(work_B,work_B,FFT3d::FORWARD);

  // group-group energy and force contribution,
  //  keep everything in reciprocal space so
  //  no inverse FFTs needed

  bigint ngridtotal = (bigint) nx_pppm * ny_pppm * nz_pppm;
  double scaleinv = 1.0/ngridtotal;
  double s2 = scaleinv*scaleinv;

  // energy

  n = 0;
  for (i = 0; i < nfft; i++) {
    e2group += s2 * greensfn[i] *
      (work_A[n]*work_B[n] + work_A[n+1]*work_B[n+1]);
    n += 2;
  }

  if (AA_flag) return;

  // multiply by Green's function and s2
  //  (only for work_A so it is not squared below)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_A[n++] *= s2 * greensfn[i];
    work_A[n++] *= s2 * greensfn[i];
  }

  // triclinic system

  if (triclinic) {
    poisson_groups_triclinic();
    return;
  }

  double partial_group;

  // force, x direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
        f2group[0] += fkx[i] * partial_group;
        n += 2;
      }

  // force, y direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
        f2group[1] += fky[j] * partial_group;
        n += 2;
      }

  // force, z direction

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
        f2group[2] += fkz[k] * partial_group;
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for group-group interactions
   for a triclinic system
 ------------------------------------------------------------------------- */

void ESP::poisson_groups_triclinic()
{
  int i,n;

  // reuse memory (already declared)

  FFT_SCALAR *work_A = work1;
  FFT_SCALAR *work_B = work2;

  double partial_group;

  // force, x direction

  n = 0;
  for (i = 0; i < nfft; i++) {
    partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
    f2group[0] += fkx[i] * partial_group;
    n += 2;
  }

  // force, y direction

  n = 0;
  for (i = 0; i < nfft; i++) {
    partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
    f2group[1] += fky[i] * partial_group;
    n += 2;
  }

  // force, z direction

  n = 0;
  for (i = 0; i < nfft; i++) {
    partial_group = work_A[n]*work_B[n+1] - work_A[n+1]*work_B[n];
    f2group[2] += fkz[i] * partial_group;
    n += 2;
  }
}

/* ----------------------------------------------------------------------
   slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void ESP::slabcorr_groups(int groupbit_A, int groupbit_B, int AA_flag)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd_slab = domain->zprd*slab_volfactor;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double qsum_A = 0.0;
  double qsum_B = 0.0;
  double dipole_A = 0.0;
  double dipole_B = 0.0;
  double dipole_r2_A = 0.0;
  double dipole_r2_B = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
      if (AA_flag) continue;

    if (mask[i] & groupbit_A) {
      qsum_A += q[i];
      dipole_A += q[i]*x[i][2];
      dipole_r2_A += q[i]*x[i][2]*x[i][2];
    }

    if (mask[i] & groupbit_B) {
      qsum_B += q[i];
      dipole_B += q[i]*x[i][2];
      dipole_r2_B += q[i]*x[i][2]*x[i][2];
    }
  }

  // sum local contributions to get total charge and global dipole moment
  //  for each group

  double tmp;
  MPI_Allreduce(&qsum_A,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_A = tmp;

  MPI_Allreduce(&qsum_B,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_B = tmp;

  MPI_Allreduce(&dipole_A,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  dipole_A = tmp;

  MPI_Allreduce(&dipole_B,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  dipole_B = tmp;

  MPI_Allreduce(&dipole_r2_A,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  dipole_r2_A = tmp;

  MPI_Allreduce(&dipole_r2_B,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  dipole_r2_B = tmp;

  // compute corrections

  const double qscale = qqrd2e * scale;
  const double efact = qscale * MY_2PI/volume;

  e2group += efact * (dipole_A*dipole_B - 0.5*(qsum_A*dipole_r2_B +
    qsum_B*dipole_r2_A) - qsum_A*qsum_B*zprd_slab*zprd_slab/12.0);

  // add on force corrections

  const double ffact = qscale * (-4.0*MY_PI/volume);
  f2group[2] += ffact * (qsum_A*dipole_B - qsum_B*dipole_A);
}
