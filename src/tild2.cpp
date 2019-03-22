/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    tabedzki@seas.upenn.edu zvicars@seas.upenn.edu
-------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "tild2.h"
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "remap_wrap.h"
#include "error.h"
#include <iostream>
#include "fft3d_wrap.h"
#include "pppm.h"

using namespace LAMMPS_NS;
using namespace MathConst;

double *TILD::uG, TILD::a_squared;

#define SMALL 0.00001
#define OFFSET 16384

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

enum{REVERSE_RHO};

/* ------------------------------------------------------------ */

TILD::TILD(LAMMPS *lmp) : KSpace(lmp),
  factors(NULL), density_brick(NULL), vdx_brick(NULL), vdy_brick(NULL), vdz_brick(NULL),
  u_brick(NULL), v0_brick(NULL), v1_brick(NULL), v2_brick(NULL), v3_brick(NULL),
  v4_brick(NULL), v5_brick(NULL), greensfn(NULL), vg(NULL), fkx(NULL), fky(NULL),
  fkz(NULL), density_fft(NULL), work1(NULL), work2(NULL), gf_b(NULL), rho1d(NULL),
  rho_coeff(NULL), drho1d(NULL), drho_coeff(NULL), sf_precoeff1(NULL), sf_precoeff2(NULL),
  sf_precoeff3(NULL), sf_precoeff4(NULL), sf_precoeff5(NULL), sf_precoeff6(NULL),
  acons(NULL), density_A_brick(NULL), density_B_brick(NULL), density_A_fft(NULL),
  density_B_fft(NULL), fft1(NULL), fft2(NULL), remap(NULL), cg(NULL), cg_peratom(NULL),
  part2grid(NULL), boxlo(NULL)
 {
    if (screen) fprintf(screen,"TILD construction...\n");
    if (logfile) fprintf(logfile,"TILD construction...\n");
  peratom_allocate_flag = 0;
  group_allocate_flag = 0;

  pppmflag = 1;
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

  density_brick = vdx_brick = vdy_brick = vdz_brick = NULL;
  density_fft = NULL;
  u_brick = NULL;
  v0_brick = v1_brick = v2_brick = v3_brick = v4_brick = v5_brick = NULL;
  greensfn = NULL;
  work1 = work2 = NULL;
  vg = NULL;
  fkx = fky = fkz = NULL;

  sf_precoeff1 = sf_precoeff2 = sf_precoeff3 =
    sf_precoeff4 = sf_precoeff5 = sf_precoeff6 = NULL;

  density_A_brick = density_B_brick = NULL;
  density_A_fft = density_B_fft = NULL;

  gf_b = NULL;
  rho1d = rho_coeff = drho1d = drho_coeff = NULL;

  fft1 = fft2 = NULL;
  remap = NULL;
  cg = NULL;
  cg_peratom = NULL;

  nmax = 0;
  part2grid = NULL;

  // define acons coefficients for estimation of kspace errors
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

  memory->create(acons,8,7,"pppm:acons");
  acons[1][0] = 2.0 / 3.0;
  acons[2][0] = 1.0 / 50.0;
  acons[2][1] = 5.0 / 294.0;
  acons[3][0] = 1.0 / 588.0;
  acons[3][1] = 7.0 / 1440.0;
  acons[3][2] = 21.0 / 3872.0;
  acons[4][0] = 1.0 / 4320.0;
  acons[4][1] = 3.0 / 1936.0;
  acons[4][2] = 7601.0 / 2271360.0;
  acons[4][3] = 143.0 / 28800.0;
  acons[5][0] = 1.0 / 23232.0;
  acons[5][1] = 7601.0 / 13628160.0;
  acons[5][2] = 143.0 / 69120.0;
  acons[5][3] = 517231.0 / 106536960.0;
  acons[5][4] = 106640677.0 / 11737571328.0;
  acons[6][0] = 691.0 / 68140800.0;
  acons[6][1] = 13.0 / 57600.0;
  acons[6][2] = 47021.0 / 35512320.0;
  acons[6][3] = 9694607.0 / 2095994880.0;
  acons[6][4] = 733191589.0 / 59609088000.0;
  acons[6][5] = 326190917.0 / 11700633600.0;
  acons[7][0] = 1.0 / 345600.0;
  acons[7][1] = 3617.0 / 35512320.0;
  acons[7][2] = 745739.0 / 838397952.0;
  acons[7][3] = 56399353.0 / 12773376000.0;
  acons[7][4] = 25091609.0 / 1560084480.0;
  acons[7][5] = 1755948832039.0 / 36229939200000.0;
  acons[7][6] = 4887769399.0 / 37838389248.0;
  };

void TILD::settings(int narg, char **arg)
{
    if (narg < 1) error->all(FLERR,"Illegal kspace_style tild command");
    accuracy_relative = fabs(force->numeric(FLERR,arg[0]));
}

TILD::~TILD(){
  return;
}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }
}


void TILD::setup(){

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with PPPMDisp");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDisp");
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

 // compute fkx,fky,fkz for my FFT grid pts

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  init_gauss();

  return;
}

void TILD::setup_grid()
{
  // free all arrays previously allocated

  deallocate();
  deallocate_peratom();

  // reset portion of global grid that each proc owns

  set_fft_parameters(nx_pppm, ny_pppm, nz_pppm,
                      nxlo_fft, nylo_fft, nzlo_fft,
                      nxhi_fft, nyhi_fft, nzhi_fft,
                      nxlo_in, nylo_in, nzlo_in,
                      nxhi_in, nyhi_in, nzhi_in,
                      nxlo_out, nylo_out, nzlo_out,
                      nxhi_out, nyhi_out, nzhi_out,
                      nlower, nupper,
                      ngrid, nfft, nfft_both,
                      shift, shiftone, order);


  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

  if (function[0]) {
    cg->ghost_notify();
    if (overlap_allowed == 0 && cg->ghost_overlap())
      error->all(FLERR,"PPPM grid stencil extends "
                 "beyond nearest neighbor processor");
    cg->setup();
  }
  if (function[1] + function[2] + function[3]) {
    cg_6->ghost_notify();
    if (overlap_allowed == 0 && cg_6->ghost_overlap())
      error->all(FLERR,"PPPM grid stencil extends "
                 "beyond nearest neighbor processor");
    cg_6->setup();
  }

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  if (function[0]) {
    compute_gf_denom(gf_b, order);
    compute_rho_coeff(rho_coeff, drho_coeff, order);
    if (differentiation_flag == 1)
      compute_sf_precoeff(nx_pppm, ny_pppm, nz_pppm, order,
                          nxlo_fft, nylo_fft, nzlo_fft,
                          nxhi_fft, nyhi_fft, nzhi_fft,
                          sf_precoeff1, sf_precoeff2, sf_precoeff3,
                          sf_precoeff4, sf_precoeff5, sf_precoeff6);
  }
  if (function[1] + function[2] + function[3]) {
    compute_gf_denom(gf_b_6, order_6);
    compute_rho_coeff(rho_coeff_6, drho_coeff_6, order_6);
    if (differentiation_flag == 1)
      compute_sf_precoeff(nx_pppm_6, ny_pppm_6, nz_pppm_6, order_6,
                          nxlo_fft_6, nylo_fft_6, nzlo_fft_6,
                          nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                          sf_precoeff1_6, sf_precoeff2_6, sf_precoeff3_6,
                          sf_precoeff4_6, sf_precoeff5_6, sf_precoeff6_6);
  }

  // pre-compute volume-dependent coeffs

  setup();
}

void TILD::compute(int i1, int i2){
  
  int i; 
  // convert atoms from box to lambda coords

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    if (function[0]) {
      cg_peratom->ghost_notify();
      cg_peratom->setup();
    }
    if (function[1] + function[2] + function[3]) {
      cg_peratom_6->ghost_notify();
      cg_peratom_6->setup();
    }
    peratom_allocate_flag = 1;
  }

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }
  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {

    if (function[0]) memory->destroy(part2grid);
    if (function[1] + function[2] + function[3]) memory->destroy(part2grid_6);
    nmax = atom->nmax;
    if (function[0]) memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
    if (function[1] + function[2] + function[3])
      memory->create(part2grid_6,nmax,3,"pppm/disp:part2grid_6");
  }

  // find grid points for all my particles 
  // distribute particles' densities on the grid
  // communication between processors and remapping two fft
  // Convolution in k-space and backtransformation 
  // communication between processors
  // calculation of forces

    particle_map(delxinv_den, delyinv_den, delzinv_den, shift_den, part2grid_den, nupper_den, nlower_den,
                 nxlo_out_den, nylo_out_den, nzlo_out_den, nxhi_out_den, nyhi_out_den, nzhi_out_den);
    make_rho_g();

  return;
}


double TILD::memory_usage(){
  return 0;
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::allocate()
{
  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick");

  memory->create(density_fft,nfft_both,"pppm:density_fft");
  memory->create(greensfn,nfft_both,"pppm:greensfn");
  memory->create(work1,2*nfft_both,"pppm:work1");
  memory->create(work2,2*nfft_both,"pppm:work2");
  memory->create(vg,nfft_both,6,"pppm:vg");
  memory->create(uG,nfft,"pppm:uG");
  memory->create(grad_uG,domain->dimension,nfft_both,"pppm:grad_uG");
  memory->create(grad_uG_hat,domain->dimension,nfft_both,"pppm:grad_uG_hat");

  if (triclinic == 0) {
    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm:fkz");
  } else {
    memory->create(fkx,nfft_both,"pppm:fkx");
    memory->create(fky,nfft_both,"pppm:fky");
    memory->create(fkz,nfft_both,"pppm:fkz");
  }

  if (differentiation_flag == 1) {
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:u_brick");

    memory->create(sf_precoeff1,nfft_both,"pppm:sf_precoeff1");
    memory->create(sf_precoeff2,nfft_both,"pppm:sf_precoeff2");
    memory->create(sf_precoeff3,nfft_both,"pppm:sf_precoeff3");
    memory->create(sf_precoeff4,nfft_both,"pppm:sf_precoeff4");
    memory->create(sf_precoeff5,nfft_both,"pppm:sf_precoeff5");
    memory->create(sf_precoeff6,nfft_both,"pppm:sf_precoeff6");

  } else {
    memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdx_brick");
    memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdy_brick");
    memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdz_brick");
  }

  // summation coeffs

  order_allocated = order;
  if (!stagger_flag) memory->create(gf_b,order,"pppm:gf_b");
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(drho1d,3,-order/2,order/2,"pppm:drho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,
                          "pppm:drho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   0,0,&tmp,collective_flag);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   0,0,&tmp,collective_flag);

  remap = new Remap(lmp,world,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
                    1,0,0,FFT_PRECISION,collective_flag);

  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;

  if (differentiation_flag == 1)
    cg = new GridComm(lmp,world,1,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  else
    cg = new GridComm(lmp,world,3,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::deallocate()
{
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
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  if (triclinic == 0) {
    memory->destroy1d_offset(fkx,nxlo_fft);
    memory->destroy1d_offset(fky,nylo_fft);
    memory->destroy1d_offset(fkz,nzlo_fft);
  } else {
    memory->destroy(fkx);
    memory->destroy(fky);
    memory->destroy(fkz);
  }

  memory->destroy(gf_b);
  if (stagger_flag) gf_b = NULL;
  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(drho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);

  // memory->destroy(uG);
  memory->destroy(grad_uG);
  memory->destroy(grad_uG_hat);

  delete fft1;
  delete fft2;
  delete remap;
  delete cg;
}


/* ----------------------------------------------------------------------
   set global size of PPPM grid = nx,ny,nz_pppm
   used for charge accumulation, FFTs, and electric field interpolation
------------------------------------------------------------------------- */

void TILD::set_grid_global()
{
  // use xprd,yprd,zprd (even if triclinic, and then scale later)
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (accuracy <= 0.0)
      error->all(FLERR,"KSpace accuracy must be > 0");
    if (q2 == 0.0)
      error->all(FLERR,"Must use kspace_modify gewald for uncharged system");
    g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
    if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/cutoff;
    else g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {

    if (differentiation_flag == 1 || stagger_flag) {

      h = h_x = h_y = h_z = 4.0/g_ewald;
      int count = 0;
      while (1) {

        // set grid dimension
        nx_pppm = static_cast<int> (xprd/h_x);
        ny_pppm = static_cast<int> (yprd/h_y);
        nz_pppm = static_cast<int> (zprd_slab/h_z);

        if (nx_pppm <= 1) nx_pppm = 2;
        if (ny_pppm <= 1) ny_pppm = 2;
        if (nz_pppm <= 1) nz_pppm = 2;

        //set local grid dimension
        int npey_fft,npez_fft;
        if (nz_pppm >= nprocs) {
          npey_fft = 1;
          npez_fft = nprocs;
        } else procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);

        int me_y = me % npey_fft;
        int me_z = me / npey_fft;

        nxlo_fft = 0;
        nxhi_fft = nx_pppm - 1;
        nylo_fft = me_y*ny_pppm/npey_fft;
        nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
        nzlo_fft = me_z*nz_pppm/npez_fft;
        nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;

        double df_kspace = compute_df_kspace();

        count++;

        // break loop if the accuracy has been reached or
        // too many loops have been performed

        if (df_kspace <= accuracy) break;
        if (count > 500) error->all(FLERR, "Could not compute grid size");
        h *= 0.95;
        h_x = h_y = h_z = h;
      }

    } else {

      double err;
      h_x = h_y = h_z = 1.0/g_ewald;

      nx_pppm = static_cast<int> (xprd/h_x) + 1;
      ny_pppm = static_cast<int> (yprd/h_y) + 1;
      nz_pppm = static_cast<int> (zprd_slab/h_z) + 1;

      err = estimate_ik_error(h_x,xprd,natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_x,xprd,natoms);
        nx_pppm++;
        h_x = xprd/nx_pppm;
      }

      err = estimate_ik_error(h_y,yprd,natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_y,yprd,natoms);
        ny_pppm++;
        h_y = yprd/ny_pppm;
      }

      err = estimate_ik_error(h_z,zprd_slab,natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_z,zprd_slab,natoms);
        nz_pppm++;
        h_z = zprd_slab/nz_pppm;
      }
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
    error->all(FLERR,"PPPM grid is too large");
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::allocate_peratom()
{
  peratom_allocate_flag = 1;

  if (differentiation_flag != 1)
    memory->create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:u_brick");

  memory->create3d_offset(v0_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v0_brick");

  memory->create3d_offset(v1_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v1_brick");
  memory->create3d_offset(v2_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v2_brick");
  memory->create3d_offset(v3_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v3_brick");
  memory->create3d_offset(v4_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v4_brick");
  memory->create3d_offset(v5_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:v5_brick");

  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;

  if (differentiation_flag == 1)
    cg_peratom =
      new GridComm(lmp,world,6,1,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                   procneigh[0][0],procneigh[0][1],procneigh[1][0],
                   procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  else
    cg_peratom =
      new GridComm(lmp,world,7,1,
                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                   nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                   procneigh[0][0],procneigh[0][1],procneigh[1][0],
                   procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::deallocate_peratom()
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

  delete cg_peratom;
}

// Need to create functionality that would loop over all places owned by 
// this processor and use that to generate the position for the gaussian. and assign them as well.
void TILD::init_gauss(){
  
  // Represents the 0 to N-1 points that are captured on this grid.
  int nlo_in, nhi_in;
  int nrange = nhi_in - nlo_in +1;
  int n_loc=0;
  int Dim = domain->dimension;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  double mdr2;
  double V = domain->xprd * domain->yprd * domain->zprd ;

  // clear 3d density arrays

  memset(&(density_A_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  memset(&(density_B_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));


  double pref = V / ( pow(2.0* sqrt(MY_PI * gauss_a2), Dim ));


  // decomposition of FFT mesh
  // global indices range from 0 to N-1
  // proc owns entire x-dimension, clumps of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own

  int npey_fft,npez_fft;
  if (nz_pppm >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny_pppm,nz_pppm,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx_pppm - 1;
  nylo_fft = me_y*ny_pppm/npey_fft;
  nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1;
  nzlo_fft = me_z*nz_pppm/npez_fft;
  nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1;

  // PPPM grid pts owned by this proc, including ghosts

  ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);

  // FFT grids owned by this proc, without ghosts
  // nfft = FFT points in FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
    (nzhi_fft-nzlo_fft+1);

  double xloca, zloca, yloca;
  // double mdr2;

  double yper, xper, zper;
  int k;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  n = 0;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    zper = (m * (zprd - 0.5)) / npez_fft;
    for (l = nylo_fft; l <= nyhi_fft; l++) {
      yper = (l * (yprd - 0.5)) / npey_fft;
      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        xper = (l * (xprd - 0.5)) / npey_fft;
        mdr2=xper*xper + yper*yper + zper*zper;
        uG[n++] = exp(-mdr2 * 0.25 / a_squared) * pref;
      }
    }
  }

  // Do the field gradient of the uG

  // Do the FFT of the Gaussian function



}



void TILD::compute_group_group(int groupbit_A, int groupbit_B, int AA_flag)
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

  cg->reverse_comm(this,REVERSE_RHO);
  brick2fft();

  // group B

  density_brick = density_B_brick;
  density_fft = density_B_fft;

  cg->reverse_comm(this,REVERSE_RHO);
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

  return;
}

void TILD::field_groups(int AA_flag){
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

  fft1->compute(work_A,work_A,1);

  // group B

  n = 0;
  for (i = 0; i < nfft; i++) {
    work_B[n++] = density_B_fft[i];
    work_B[n++] = ZEROF;
  }

  fft1->compute(work_B,work_B,1);

  // group-group energy and force contribution,
  //  keep everything in reciprocal space so
  //  no inverse FFTs needed

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
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
    error->all(FLERR,"TILD doesn't support triclinic yet.");
    return;
  }

  double partial_group;


  fft1->compute(work_A, work_A,1);


  /************************************************************************
   * OLD RIGGLEMAN CODE
  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////
  //Set forces for local particles to 0
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
  //Set forces for ghost particles associated with current proc to 0
  for ( i=0 ; i<total_ghost ; i++ ) {
    id = ghost_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
 
  //Sets gradients to 0 
  for ( i=0 ; i<ML ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwA[j][i] = gradwB[j][i] = gradwC[j][i] = gradwC[j][i] = gradwLC[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////
  
  // A acting on B, C //
  fftw_fwd( rho[0] , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ ) ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      if ( chiAB != 0.0 )
        gradwB[j][i] += tmp[i] * chiAB / rho0 ;
      if ( chiAC != 0.0 )
        gradwC[j][i] += tmp[i] * chiAC / rho0 ;
    }
  }
  ***************************************************************************/

  

}


// void TILD::get_k_alias(int id, double k[domain->dimension]){

// }

void TILD::field_gradient(FFT_SCALAR *in, 
                          FFT_SCALAR *out, int flag,
                          int dir)
{
  int i; 
  int Dim = domain->dimension;
  double k2, kv[Dim];
  fft1->compute(in, in, 1);

  for ( i = 0; i < GRID_POINTS ; i++){
    TILD::get_k_alias;
  }

  fft2->compute(in, out, -1);

}


int TILD::factorable(int n)
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

double TILD::get_k_alias(int id, double k[]){
  double kmag = 0.0;
  int dim = domain->dimensions; 
  int i, id2, n[Dim] , j , has_nyquist = 0;
  for ( i=0 ; i<Dim ; i++ )
    if ( Nx[i] % 2 == 0 )
      has_nyquist = 1;

  id2 = unstack_stack( id ) ;
  unstack(id2, n);


  if ( Nx[0] % 2 == 0 && n[0] == Nx[0] / 2 )
    k[0] = 0.0 ;
  else if ( double(n[0]) < double(Nx[0]) / 2.)
   k[0] = 2*PI*double(n[0])/L[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/L[0];

  if (Dim>1) {
    if ( Nx[1] % 2 == 0 && n[1] == Nx[1] / 2 )
      k[1] = 0.0 ;
    else if ( double(n[1]) < double(Nx[1]) / 2.)
      k[1] = 2*PI*double(n[1])/L[1];
    else
      k[1] = 2*PI*double(n[1]-Nx[1])/L[1];
  }

  if (Dim==3) {
    if ( Nx[2] % 2 == 0 && n[2] == Nx[2] / 2 )
      k[2] = 0.0 ;
    else if ( double(n[2]) < double(Nx[2]) / 2.)
      k[2] = 2*PI*double(n[2])/L[2];
    else
      k[2] = 2*PI*double(n[2]-Nx[2])/L[2];
  }

  // Kills off the Nyquist modes
  if ( id2 != 0 && has_nyquist ) {
    for ( i=0 ; i<Dim ; i++ ) {
      if ( k[i] == 0.0 ) {
        for ( j=0 ; j<Dim ; j++ ) k[j] = 0.0 ;
        kmag = 0.0;
        break;
      }
    }
  }
  
  for (i=0; i<Dim; i++) kmag += k[i]*k[i];

  return kmag;

}
