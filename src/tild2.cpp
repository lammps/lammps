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
#include <complex>
#include "group.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001
#define OFFSET 16384
#define PI 3.141592653589793238462643383279
#define MAXORDER   7
#define I std::complex<double>(0.0, 1.0)

enum{REVERSE_RHO};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif


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
  deallocate();
  deallocate_groups();
  deallocate_peratom();

}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use PPPMDisp with 2d simulation");
  if (comm->style != 0)
    error->universe_all(FLERR,"PPPMDisp can only currently be used with "
                        "comm_style brick");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with PPPMDisp");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab PPPMDisp");
  }

  if (order > MAXORDER || order_6 > MAXORDER) {
    char str[128];
    sprintf(str,"TILD coulomb order cannot be greater than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated
  deallocate();
  deallocate_groups();
  deallocate_peratom();

  memory->create(param, group->ngroup, group->ngroup,"tild:param");

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

  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  init_gauss();

  return;
}

void TILD::setup_grid()
{
  // free all arrays previously allocated

  deallocate();
  deallocate_peratom();
  if (group_allocate_flag) deallocate_groups();

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


  set_grid_global();
  set_grid_local();

  
  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate peratom() or group(), will be allocated when needed


  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

    cg->ghost_notify();
    if (overlap_allowed == 0 && cg->ghost_overlap())
      error->all(FLERR,"PPPM grid stencil extends "
                 "beyond nearest neighbor processor");
    cg->setup();

  // pre-compute volume-dependent coeffs

  setup();
}

void TILD::compute(int eflag, int vflag){

  double density;
  output->thermo->evaluate_keyword("density",&density);
  
  int i; 

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // convert atoms from box to lamda coords

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
      cg_peratom->ghost_notify();
      cg_peratom->setup();
    peratom_allocate_flag = 1;
  }

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }
  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {

    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
  }

  particle_map(delxinv, delyinv, delzinv, shift, part2grid, nupper, nlower,
                nxlo_out, nylo_out, nzlo_out, nxhi_out, nyhi_out, nzhi_out);

  make_rho_none();

  cg->reverse_comm(this, REVERSE_RHO);

  brick2fft_none();
  // convert atoms from box to lambda coords

  // if (eflag || vflag) ev_setup(eflag,vflag);
  // else evflag = evflag_atom = eflag_global = vflag_global =
  //        eflag_atom = vflag_atom = 0;

  // if (evflag_atom && !peratom_allocate_flag) {
  //   allocate_peratom();
  //   if (function[0]) {
  //     cg_peratom->ghost_notify();
  //     cg_peratom->setup();
  //   }
  //   if (function[1] + function[2] + function[3]) {
  //     cg_peratom_6->ghost_notify();
  //     cg_peratom_6->setup();
  //   }
  //   peratom_allocate_flag = 1;
  // }

  // if (triclinic == 0) boxlo = domain->boxlo;
  // else {
  //   boxlo = domain->boxlo_lamda;
  //   domain->x2lamda(atom->nlocal);
  // }
  // // extend size of per-atom arrays if necessary

  // if (atom->nmax > nmax) {

  //   if (function[0]) memory->destroy(part2grid);
  //   if (function[1] + function[2] + function[3]) memory->destroy(part2grid_6);
  //   nmax = atom->nmax;
  //   if (function[0]) memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
  //   if (function[1] + function[2] + function[3])
  //     memory->create(part2grid_6,nmax,3,"pppm/disp:part2grid_6");
  // }

  // find grid points for all my particles 
  // distribute particles' densities on the grid
  // communication between processors and remapping two fft
  // Convolution in k-space and backtransformation 
  // communication between processors
  // calculation of forces

  //   make_rho_g();

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

  memory->destroy(uG);
  memory->destroy(temp);
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

        double qopt = compute_qopt();

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
  int Dim = domain->dimension;

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
  memory->create(gradWgroup, group->ngroup, Dim, nfft, "tild:gradWgroup");
  memory->create(ktmp, nfft, "tild:ktmp");
  memory->create(ktmp2, nfft, "tild:ktmp2");
  memory->create(tmp, nfft, "tild:tmp");

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


  // double pref = V / ( pow(2.0* sqrt(MY_PI * gauss_a2), Dim ));


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
        uG[n++] = exp(-mdr2 * 0.25 / a_squared) ;
      }
    }
  }

  // Do the field gradient of the uG
  field_gradient(uG, grad_uG, 0);

  // Do the FFT of the Gaussian function
  for (int i = 0; i < Dim; i++)
  fft1->compute(grad_uG[i], grad_uG_hat[i], 1);

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

  // make_rho_groups(groupbit_A,groupbit_B,AA_flag);
  
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

  // poisson_groups(AA_flag);

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

  // if (slabflag == 1)
  //   slabcorr_groups(groupbit_A, groupbit_B, AA_flag);

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

void TILD::field_gradient(FFT_SCALAR *in, 
                          FFT_SCALAR **out, int flag)
{
  int i; 
  int Dim = domain->dimension;
  double k2, kv[Dim];
  fft1->compute(in, in, 1);

  get_k_alias();

  for (i = 0; i < Dim; i++)
  fft2->compute(in, out[i], -1);

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

void TILD::get_k_alias(){
  int Dim = domain->dimension; 
  int k[Dim];
  int x, y, z;
  int n=0;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  // id2 = unstack_stack( id ) ;
  // unstack(id2, n);
  for (z = nzlo_fft; z <= nzhi_fft; z++) {
    for (y = nylo_fft; y <= nyhi_fft; y++) {
      for (x = nxlo_fft; x <= nxhi_fft; x++) {

          if (nx_pppm % 2 == 0 && x == nx_pppm / 2)
              k[0] = 0.0;
          else if (double(x) < double(nx_pppm) / 2.)
              k[0] = 2 * PI * double(x) / xprd;
          else
              k[0] = 2 * PI * double(y - nx_pppm) / xprd;

          if (ny_pppm % 2 == 0 && y == ny_pppm / 2)
              k[1] = 0.0;
          else if (double(y) < double(ny_pppm) / 2.)
              k[1] = 2 * PI * double(y) / yprd;
          else
              k[1] = 2 * PI * double(y - ny_pppm) / yprd;

          if (Dim == 3) {
              if (nz_pppm % 2 == 0 && z == nz_pppm / 2)
                  k[2] = 0.0;
              else if (double(z) < double(nz_pppm) / 2.)
                  k[2] = 2 * PI * double(z) / zprd;
              else
                  k[2] = 2 * PI * double(z - nz_pppm) / zprd;
          }

          grad_uG[0][n] = uG[n] * k[0];
          grad_uG[1][n] = uG[n] * k[1];
          if (Dim == 3)
            grad_uG[2][n] = uG[n] * k[2];
          n++;
      }
    }
  }



}

void TILD::particle_map(double delx, double dely, double delz,
                             double sft, int** p2g, int nup, int nlow,
                             int nxlo, int nylo, int nzlo,
                             int nxhi, int nyhi, int nzhi)
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delx+sft) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*dely+sft) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delz+sft) - OFFSET;

    p2g[i][0] = nx;
    p2g[i][1] = ny;
    p2g[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlow < nxlo || nx+nup > nxhi ||
        ny+nlow < nylo || ny+nup > nyhi ||
        nz+nlow < nzlo || nz+nup > nzhi)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute TILD");
}


void TILD::particle_map_c(double delx, double dely, double delz,
                               double sft, int** p2g, int nup, int nlow,
                               int nxlo, int nylo, int nzlo,
                               int nxhi, int nyhi, int nzhi)
{
  particle_map(delx, dely, delz, sft, p2g, nup, nlow,
               nxlo, nylo, nzlo, nxhi, nyhi, nzhi);
}

int TILD::modify_param(int narg, char** arg)
{
  int i;

    if (strcmp(arg[0], "tild/params") == 0) {

  if (domain->box_exist == 0)
            error->all(FLERR, "TILD command before simulation box is defined");
        if (narg < 3)
            error->all(FLERR, "Illegal kspace_modify tild command");

        if (strcmp(arg[1], "all") == 0) {
            if (narg != 3)
                error->all(FLERR, "Illegal kspace_modify tild command");

    kappa = atof(arg[1]);
            param[0][0] = atof(arg[1]);
        } else {
            if (narg != 4)
                error->all(FLERR, "Illegal kspace_modify tild command");
    int igroup1 = group->find(arg[1]);
    int igroup2 = group->find(arg[2]);
        
            if (igroup1 == -1) {
      error->all(FLERR, "group1 not found in kspace_modify tild command");
    }
            if (igroup2 == -1) {
      error->all(FLERR, "group2 not found in kspace_modify tild command");
  }
    if (igroup1 == 0 || igroup2 == 0)
      error->all(FLERR, "all group specified in 'group1 group2 param' format");

    param[igroup1][igroup2] = param[igroup2][igroup1] = atof(arg[3]);
}
  } else 
        error->all(FLERR, "Illegal kspace_modify tild command");
  return narg;
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition
   remap density from 3d brick decomposition to FFTdecomposition
   for coulomb interaction or dispersion interaction with geometric
   mixing
------------------------------------------------------------------------- */
void TILD::brick2fft(int nxlo_i, int nylo_i, int nzlo_i,
                         int nxhi_i, int nyhi_i, int nzhi_i,
                         FFT_SCALAR*** dbrick, FFT_SCALAR* dfft, FFT_SCALAR* work,
                         LAMMPS_NS::Remap* rmp)
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_i; iz <= nzhi_i; iz++)
    for (iy = nylo_i; iy <= nyhi_i; iy++)
      for (ix = nxlo_i; ix <= nxhi_i; ix++)
        dfft[n++] = dbrick[iz][iy][ix];

  rmp->perform(dfft,dfft,work);
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */
void TILD::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
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
void TILD::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
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
void TILD::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  if (flag == REVERSE_RHO) {
    FFT_SCALAR *src = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */
void TILD::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  if (flag == REVERSE_RHO) {
    FFT_SCALAR *dest = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  }
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
------------------------------------------------------------------------- */
void TILD::procs2grid2d(int nprocs, int nx, int ny, int *px, int *py)
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
   set the FFT parameters
------------------------------------------------------------------------- */
void TILD::set_fft_parameters(int& nx_p,int& ny_p,int& nz_p,
                                   int& nxlo_f,int& nylo_f,int& nzlo_f,
                                   int& nxhi_f,int& nyhi_f,int& nzhi_f,
                                   int& nxlo_i,int& nylo_i,int& nzlo_i,
                                   int& nxhi_i,int& nyhi_i,int& nzhi_i,
                                   int& nxlo_o,int& nylo_o,int& nzlo_o,
                                   int& nxhi_o,int& nyhi_o,int& nzhi_o,
                                   int& nlow, int& nupp,
                                   int& ng, int& nf, int& nfb,
                                   double& sft,double& sftone, int& ord)
{
  // global indices of PPPM grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that I own without ghost cells
  // for slab PPPM, assign z grid as if it were not extended

  nxlo_i = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_p);
  nxhi_i = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_p) - 1;

  nylo_i = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_p);
  nyhi_i = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_p) - 1;

  nzlo_i = static_cast<int>
      (comm->zsplit[comm->myloc[2]] * nz_p/slab_volfactor);
  nzhi_i = static_cast<int>
      (comm->zsplit[comm->myloc[2]+1] * nz_p/slab_volfactor) - 1;

  // nlow,nupp = stencil size for mapping particles to PPPM grid

  nlow = -(ord-1)/2;
  nupp = ord/2;

  // sft values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (ord % 2) sft = OFFSET + 0.5;
  else sft = OFFSET;
  if (ord % 2) sftone = 0.0;
  else sftone = 0.5;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 + qdist
  //   qdist = offset due to TIP4P fictitious charge
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab PPPM, assign z grid as if it were not extended

  double *prd,*sublo,*subhi;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double dist[3];
  double cuthalf = 0.5*neighbor->skin + qdist;
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = cuthalf;
  else {
    dist[0] = cuthalf/domain->prd[0];
    dist[1] = cuthalf/domain->prd[1];
    dist[2] = cuthalf/domain->prd[2];
  }

  int nlo,nhi;

  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_p/xprd + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_p/xprd + sft) - OFFSET;
  nxlo_o = nlo + nlow;
  nxhi_o = nhi + nupp;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_p/yprd + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_p/yprd + sft) - OFFSET;
  nylo_o = nlo + nlow;
  nyhi_o = nhi + nupp;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_p/zprd_slab + sft) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_p/zprd_slab + sft) - OFFSET;
  nzlo_o = nlo + nlow;
  nzhi_o = nhi + nupp;

  // for slab PPPM, change the grid boundary for processors at +z end
  //   to include the empty volume between periodically repeating slabs
  // for slab PPPM, want charge data communicated from -z proc to +z proc,
  //   but not vice versa, also want field data communicated from +z proc to
  //   -z proc, but not vice versa
  // this is accomplished by nzhi_i = nzhi_o on +z end (no ghost cells)

  if (slabflag && (comm->myloc[2] == comm->procgrid[2]-1)) {
    nzhi_i = nz_p - 1;
    nzhi_o = nz_p - 1;
  }

  // decomposition of FFT mesh
  // global indices range from 0 to N-1
  // proc owns entire x-dimension, clump of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own

  int npey_fft,npez_fft;
  if (nz_p >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else procs2grid2d(nprocs,ny_p,nz_p,&npey_fft,&npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_f = 0;
  nxhi_f = nx_p - 1;
  nylo_f = me_y*ny_p/npey_fft;
  nyhi_f = (me_y+1)*ny_p/npey_fft - 1;
  nzlo_f = me_z*nz_p/npez_fft;
  nzhi_f = (me_z+1)*nz_p/npez_fft - 1;

  // PPPM grid for this proc, including ghosts

  ng = (nxhi_o-nxlo_o+1) * (nyhi_o-nylo_o+1) *
    (nzhi_o-nzlo_o+1);

  // FFT arrays on this proc, without ghosts
  // nfft = FFT points in FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nf = (nxhi_f-nxlo_f+1) * (nyhi_f-nylo_f+1) *
    (nzhi_f-nzlo_f+1);
  int nfft_brick = (nxhi_i-nxlo_i+1) * (nyhi_i-nylo_i+1) *
    (nzhi_i-nzlo_i+1);
  nfb = MAX(nf,nfft_brick);

}

/* ----------------------------------------------------------------------
   estimate kspace force error for ik method
------------------------------------------------------------------------- */
double TILD::estimate_ik_error(double h, double prd, bigint natoms)
{
  double sum = 0.0;
  if (natoms == 0) return 0.0;
  for (int m = 0; m < order; m++)
    sum += acons[order][m] * pow(h*g_ewald,2.0*m);
  double value = q2 * pow(h*g_ewald,(double)order) *
    sqrt(g_ewald*prd*sqrt(MY_2PI)*sum/natoms) / (prd*prd);

  return value;
}

/* ----------------------------------------------------------------------
 allocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */
void TILD::allocate_groups()
{
  group_allocate_flag = 1;

  memory->create3d_offset(density_A_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_A_brick");
  memory->create3d_offset(density_B_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_B_brick");
  memory->create(density_A_fft,nfft_both,"pppm:density_A_fft");
  memory->create(density_B_fft,nfft_both,"pppm:density_B_fft");
}

/* ----------------------------------------------------------------------
 deallocate group-group memory that depends on # of K-vectors and order
 ------------------------------------------------------------------------- */
void TILD::deallocate_groups()
{
  group_allocate_flag = 0;

  memory->destroy3d_offset(density_A_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(density_B_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(density_A_fft);
  memory->destroy(density_B_fft);
}

/* ----------------------------------------------------------------------
 compute estimated kspace force error
------------------------------------------------------------------------- */
double TILD::compute_df_kspace()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd*slab_volfactor;
  bigint natoms = atom->natoms;
  double df_kspace = 0.0;
  // if (differentiation_flag == 1 || stagger_flag) {
  //   double qopt = compute_qopt();
  //   df_kspace = sqrt(qopt/natoms)*q2/(xprd*yprd*zprd_slab);
  // } else {
  //   double lprx = estimate_ik_error(h_x,xprd,natoms);
  //   double lpry = estimate_ik_error(h_y,yprd,natoms);
  //   double lprz = estimate_ik_error(h_z,zprd_slab,natoms);
  //   df_kspace = sqrt(lprx*lprx + lpry*lpry + lprz*lprz) / sqrt(3.0);
  // }
  return df_kspace;
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */
void TILD::brick2fft()
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
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- case when mixing rules don't apply
------------------------------------------------------------------------- */
void TILD::make_rho_none()
{
  int k,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  int ngroups = group->ngroup;
  // clear 3d density array
  // for (k = 0; k < nsplit_alloc; k++)
  //   memset(&(density_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
  //          ngrid_6*sizeof(FFT_SCALAR));
  for (k = 0; k < ngroups; k++)
    memset(&(density_brick_types[k][nzlo_out][nylo_out][nxlo_out]),0,
           ngrid*sizeof(FFT_SCALAR));


  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int type;
  double **x = atom->x;
  int *mask = atom->mask;
  int *mass= atom->mass;
  int nlocal = atom->nlocal;
  memory->create(groupbits, ngroups, "tild:groupbits");
  for (k = 0; k < ngroups; k++) groupbits[k] = group->bitmask[k];

  for (int i = 0; i < nlocal; i++) {

    // do the following for all 4 grids
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;
    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);
    type = atom->type[i];
    z0 = delvolinv * mass[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          w = x0*rho1d[0][l];
          for (k = 0; k < group->ngroup; k++)
            if (mask[i] & groupbits[k])
            density_brick_types[k][mz][my][mx] += w;
            // density_brick_types[k][mz][my][mx] += w*B[nsplit*type + k];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition
   remap density from 3d brick decomposition to FFTdecomposition
   for dispersion with special case
------------------------------------------------------------------------- */

void TILD::brick2fft_none()
{
  int k,n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  for (k = 0; k<group->ngroup; k++) {
    n = 0;
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++)
          density_fft_types[k][n++] = density_brick_types[k][iz][iy][ix];
  }

  for (k=0; k<group->ngroup; k++)
    remap->perform(density_fft_types[k],density_fft_types[k],work1);
}

/* ----------------------------------------------------------------------
   Compute qopt for Coulomb interactions
------------------------------------------------------------------------- */
double TILD::compute_qopt()
{
  double qopt;
  if (differentiation_flag == 1) {
    qopt = compute_qopt_ad();
  } else {
    qopt = compute_qopt_ik();
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ik differentiation scheme and Coulomb interaction
------------------------------------------------------------------------- */
double TILD::compute_qopt_ik()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3,dot1,dot2;

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += sx*sy*sz*sx*sy*sz/dot2*4.0*4.0*MY_PI*MY_PI;
                sum2 += u2*sx*sy*sz*4.0*MY_PI/dot2*dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ad differentiation scheme and Coulomb interaction
------------------------------------------------------------------------- */
double TILD::compute_qopt_ad()
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;


  double unitkx = (2.0*MY_PI/xprd);
  double unitky = (2.0*MY_PI/yprd);
  double unitkz = (2.0*MY_PI/zprd_slab);

  int nx,ny,nz,kper,lper,mper;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double u2, sqk;
  double sum1,sum2,sum3,sum4,dot2;

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {

          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*pow(qx/g_ewald,2.0));
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*pow(qy/g_ewald,2.0));
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*pow(qz/g_ewald,2.0));
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order);

                dot2 = qx*qx+qy*qy+qz*qz;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += sx*sy*sz*sx*sy*sz/dot2*4.0*4.0*MY_PI*MY_PI;
                sum2 += sx*sy*sz * u2*4.0*MY_PI;
                sum3 += u2;
                sum4 += dot2*u2;
              }
            }
          }
          sum2 *= sum2;
          qopt += sum1 - sum2/(sum3*sum4);
        }
      }
    }
  }
  return qopt;
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */
void TILD::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                              const FFT_SCALAR &dz, int ord,
                             FFT_SCALAR **rho_c, FFT_SCALAR **r1d)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-ord)/2; k <= ord/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = ord-1; l >= 0; l--) {
      r1 = rho_c[l][k] + r1*dx;
      r2 = rho_c[l][k] + r2*dy;
      r3 = rho_c[l][k] + r3*dz;
    }
    r1d[0][k] = r1;
    r1d[1][k] = r2;
    r1d[2][k] = r3;
  }
}

void TILD::accumulate_gradient(){

  int Dim = domain->dimension;
  double rho0;
  output->thermo->evaluate_keyword("density",&rho0);
  bool do_fft = false;
  for (int i = 0; i < group->ngroup; i++){

    for (int j = 0; j < group->ngroup; j++){
      if (fabs(param[i][j]) > 1e-8){
        do_fft = true; 
        break;
      }
    }

    if (do_fft){
      fft1->compute(density_fft_types[i], ktmp, 1);

      for (int j = 0; j < Dim; j++) {
        for (int k = 0; k < nfft; k++) {
          gradWgroup[i][j][k] += grad_uG_hat[j][k] * ktmp[k];
        }

        fft1->compute(ktmp, tmp, -1);
      }

      for (int i2 = 0; i2 < group->ngroup; i++) {
        if (fabs(param[i][i2]) > 1e-8) {
          for (int j = 0; j < Dim; j++) {
            for (int k = 0; k < nfft; k++) {
              gradWgroup[i][j][k] += tmp[k] * param[i][i2] / rho0;
            }
          }
        }
      }
    }
  }
}