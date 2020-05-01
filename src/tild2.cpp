/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    Andrew Santos
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
#include <fstream>
#include "fft3d_wrap.h"
#include "pppm.h"
#include "group.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"
#include <fstream>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#define SMALL 0.00001
#define OFFSET 16384
#define PI 3.141592653589793238462643383279
#define MAXORDER   7
#define MAXPARAM   3 // Update as new potentials are introduced
#define MAX_GROUP 32

enum{REVERSE_RHO, REVERSE_RHO_NONE};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM, FORWARD_NONE};

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
  if (comm->me == 0) {
    if (screen) fprintf(screen, "TILD construction...\n");
    if (logfile) fprintf(logfile, "TILD construction...\n");
  }
  csumflag = 0;
  B = NULL;
  cii = NULL;
  csumi = NULL;
  peratom_allocate_flag = 0;

  group_allocate_flag = 0;

  nstyles = 2;  // total number of sytles

  pppmflag = 0;
  group_group_enable = 0;
  tildflag = 1;
  triclinic = domain->triclinic;
  triclinic_support = 0;
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
  density_brick_types = NULL;
  density_fft_types = NULL;
  kappa_density = NULL;
  density_of_types_fft_ed = NULL;
  density_of_potentials_fft_ed = NULL;
  gradWtype = NULL;
  potent = NULL;
  potent_hat = NULL;
  grad_potent = NULL;
  grad_potent_hat = NULL;
  setflag = NULL;
  potent_type_map = NULL;
  chi = NULL;
  a2 = NULL;
  rp = NULL;
  xi = NULL;
  temp = NULL;
  ktmp = NULL;
  ktmp2 = NULL;
  ktmpi = ktmpj = NULL;
  ktmp2i = ktmp2j = NULL;
  tmp = NULL;
  u_brick = NULL;
  v0_brick = v1_brick = v2_brick = v3_brick = v4_brick = v5_brick = NULL;
  greensfn = NULL;
  work1 = work2 = NULL;
  worki = workj = NULL;
  vg = NULL;
  vg_hat = NULL;
  fkx = fky = fkz = NULL;
  fkx2 = fky2 = fkz2 = NULL;

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
  specified_all_group = 0;

  rho0 = 0.0;
  set_rho0 = 0.0;
  old_volume = 0.0;
  nmax = 0;
  sub_flag  = 1;
  mix_flag  = 1;
  subtract_rho0 = 0;
  set_rho0_flag = 0;
  norm_flag = 1;
  part2grid = NULL;

  // define acons coefficients for estimation of kspace errors
  // see JCP 109, pg 7698 for derivation of coefficients
  // higher order coefficients may be computed if needed

/*
  for (int i =0; i < MAX_GROUP; i++){
    for (int j = 0; j < MAX_GROUP; j++) {
        group_with_potential[i][j] = -1;
    }
  }
*/
  int ntypes = atom->ntypes;
  memory->create(potent_type_map,nstyles+1,ntypes+1,ntypes+1,"tild:potent_type_map");  

  memory->create(chi,ntypes+1,ntypes+1,"tild:chi");
  memory->create(a2,ntypes+1,ntypes+1,"tild:a2"); // gaussian parameter
  memory->create(xi,ntypes+1,ntypes+1,"tild:xi"); // erfc parameter
  memory->create(rp,ntypes+1,ntypes+1,"tild:rp"); // erfc parameter

  for (int i = 0; i <= ntypes; i++) {
    for (int j = 0; j <= ntypes; j++) {
      potent_type_map[0][i][j] = 1; // style 0 is 1 if no tild potential is used
      for (int k = 1; k <= nstyles; k++) {
        potent_type_map[k][i][j] = 0; // style is set to 1 if it is used by type-type interaction
      }
      chi[i][j] = 0.0;
      a2[i][j] = 0;
      xi[i][j] = 0;
      rp[i][j] = 0;
    }
  }

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
    grid_res = fabs(force->numeric(FLERR, arg[0]));
}

TILD::~TILD(){
  delete [] factors;
  deallocate();
  deallocate_groups();
  deallocate_peratom();
  memory->destroy(part2grid);
  part2grid = NULL;

  // check whether cutoff and pair style are set

  triclinic = domain->triclinic;
  pair_check();
  memory->destroy(potent_type_map);
  memory->destroy(chi);
  memory->destroy(a2);
  memory->destroy(rp);
  memory->destroy(xi);

}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use TILD with 2d simulation");
  if (differentiation_flag != 0)
    error->all(FLERR, "Cannot use analytic differentiation with TILD.");
  if (comm->style != 0)
    error->universe_all(FLERR,"TILD can only currently be used with "
                        "comm_style brick");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with TILD ");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab TILD");
  }

  if (order > MAXORDER || order_6 > MAXORDER) {
    char str[128];
    sprintf(str,"TILD coulomb order cannot be greater than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated
      //fprintf(screen,"deall\n");
  deallocate();
  deallocate_groups();
  deallocate_peratom();

      //fprintf(screen,"set_grid\n");
  set_grid();

      //fprintf(screen,"setup_grid\n");
  setup_grid();

}

void TILD::setup(){

  // ADD A SECTION TO CALCULATE RHO0

  if (specified_all_group == 0){
    if (comm->me == 0)
      error->warning(
          FLERR,
          "No groups specified for the total density. Using the all group.");
    start_group_ind = 0;
  }
  else start_group_ind = 1;

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with TILD");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab TILD");
  }

  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab TILD
  // z dimension for 3d TILD is zprd since slab_volfactor = 1.0

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

  double per;
    int i, j, k, n;

    for (i = nxlo_fft; i <= nxhi_fft; i++) {
      per = i - nx_pppm*(2*i/nx_pppm);
      fkx[i] = unitkx*per;
      j = (nx_pppm - i) % nx_pppm;
      per = j - nx_pppm*(2*j/nx_pppm);
      fkx2[i] = unitkx*per;
    }

    for (i = nylo_fft; i <= nyhi_fft; i++) {
      per = i - ny_pppm*(2*i/ny_pppm);
      fky[i] = unitky*per;
      j = (ny_pppm - i) % ny_pppm;
      per = j - ny_pppm*(2*j/ny_pppm);
      fky2[i] = unitky*per;
    }

    for (i = nzlo_fft; i <= nzhi_fft; i++) {
      per = i - nz_pppm*(2*i/nz_pppm);
      fkz[i] = unitkz*per;
      j = (nz_pppm - i) % nz_pppm;
      per = j - nz_pppm*(2*j/nz_pppm);
      fkz2[i] = unitkz*per;
    }

  delvolinv = delxinv*delyinv*delzinv;

  int ind = 0;
  if (sub_flag == 1) subtract_rho0 = 1;
  if (norm_flag == 1) normalize_by_rho0 = 1;

  if (mix_flag == 1) {
    int ntypes = atom->ntypes;
    for (int itype = 1; itype <= ntypes; itype++) {
      for (int jtype = itype+1; jtype <= ntypes; jtype++) {
        if (potent_type_map[0][itype][itype] ==1 || potent_type_map[0][jtype][jtype] == 1 ) {
          potent_type_map[0][itype][jtype] = 1;
          for (int istyle = 1; istyle <= nstyles; istyle++) 
            potent_type_map[istyle][itype][jtype] = potent_type_map[istyle][itype][itype];
        } else {
          potent_type_map[0][itype][jtype] = 0;
          for (int istyle = 1; istyle <= nstyles; istyle++) 
            // assume it is of type istyle, but it does the convolution
            potent_type_map[istyle][itype][jtype] = -1; 
        }
      }
    }
  }
  init_cross_potentials();
  //init_gauss();
  vir_func_init();

  return;
}

void TILD::vir_func_init() {
  int z, y, x, i, j, n;
  int Dim = domain->dimension;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double zprd_slab = zprd * slab_volfactor;
  double k[Dim];
  double scale_inv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);
  double delx = xprd/nx_pppm;
  double dely = yprd/ny_pppm;
  double delz = zprd/nz_pppm;
  int ntypes = atom->ntypes;

  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      // Skip if type cross-interaction does not use density/tild
      if ( potent_type_map[0][itype][jtype] == 1) continue;
      int loc = itype*(jtype+1)/2;

      n = 0;
      for (z = nzlo_fft; z <= nzhi_fft; z++) {
        // PBC
        if (double(z) < double(nz_pppm) / 2.0)
          k[2] = double(z) * delz;
        else
          k[2] = -double(nz_pppm - z) * delz;

        for (y = nylo_fft; y <= nyhi_fft; y++) {
          if (double(y) < double(ny_pppm) / 2.0)
            k[1] = double(y) * dely;
          else
            k[1] = -double(ny_pppm - y) * dely;

          for (x = nxlo_fft; x <= nxhi_fft; x++) {
            if (double(x) < double(nx_pppm) / 2.0)
              k[0] = double(x) * delx;
            else
              k[0] = -double(nx_pppm - x) * delx;

            vg[loc][0][n] = k[0] * -grad_potent[loc][0][n];
            vg[loc][1][n] = k[1] * -grad_potent[loc][1][n];
            vg[loc][2][n] = k[2] * -grad_potent[loc][2][n];
            vg[loc][3][n] = k[1] * -grad_potent[loc][0][n];
            vg[loc][4][n] = k[2] * -grad_potent[loc][0][n];
            vg[loc][5][n] = k[2] * -grad_potent[loc][1][n];

            }
          }
        n++;
      }

      for (i = 0; i < 6; i++){
        int loc = itype*(jtype+1)/2;
        n=0;
        for (j = 0; j < nfft; j++){
          work1[n++] = vg[loc][i][j];
          work1[n++] = ZEROF;
        }
        fft1->compute(work1, vg_hat[loc][i], 1);
       //fprintf(screen,"vghatting %f %f\n", vg_hat[loc][i][0], grad_potent[loc][1][0]);
        for (j = 0; j < 2 * nfft; j++) {
          vg_hat[loc][i][j] *= scale_inv;
        }
      }
    }
  }
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


  // set_grid_global();
  // set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

  if ( set_rho0_flag ) {
    rho0 = set_rho0;
    old_volume = domain->xprd * domain->yprd * domain->zprd;
  } else {
    rho0 = calculate_rho0();
  }

  // pre-compute volume-dependent coeffs
  compute_rho_coeff(rho_coeff, drho_coeff, order);
  cg->ghost_notify();
  if (overlap_allowed == 0 && cg->ghost_overlap())
    error->all(FLERR,"PPPM grid stencil extends "
               "beyond nearest neighbor processor");
  cg->setup();

  setup();
}

/*
void TILD::precompute_density_fft() {
  int n = 0;
  double scale_inv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);
  int ntypes = atom->ntypes;

  
      fprintf(screen,"precomp\n");
  for ( int ktype = 0; ktype < ntypes; ktype++) {

    memset(&(density_of_potentials_fft_ed[ktype][0]), 0,
           ngrid * sizeof(FFT_SCALAR));

    if ( potent_type_map[0][ktype][ktype] == 1) continue;
    for (int k = 0; k < nfft; k++) {
      work1[n++] = density_fft_types[ktype][k];
      work1[n++] = ZEROF;
    }

    // FFT the density to k-space
    fft1->compute(work1, work1, 1);

    n = 0;
    for (int k = 0; k < nfft; k++) {
      work1[n] *= scale_inv;
      density_of_types_fft_ed[ktype][n] += work1[n];
      n++;
      work1[n] *= scale_inv;
      density_of_types_fft_ed[ktype][n] += work1[n];
      n++;
    }
  }
}

*/
void TILD::compute(int eflag, int vflag){

  double density;
  if (domain->box_change) {
    if ( set_rho0_flag ) {
      rho0 *=  old_volume / (domain->xprd * domain->yprd * domain->zprd);
    } else {
      rho0 = calculate_rho0();
    }
    old_volume = domain->xprd * domain->yprd * domain->zprd;
  }
  
  int i; 

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // convert atoms from box to lamda coords

  
      //fprintf(screen,"e v flag %d %d\n", eflag, vflag);
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;
      //fprintf(screen,"e v flag %d %d\n", eflag_global, vflag_global);

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

      //fprintf(screen,"make rho none\n");
  make_rho_none();

      //fprintf(screen,"revcom\n");
  cg->reverse_comm(this, REVERSE_RHO_NONE);
      //fprintf(screen,"brick2fft_none\n");

  brick2fft_none();

      //fprintf(screen,"accumulate_gradient\n");
  accumulate_gradient();

  cg->forward_comm(this, FORWARD_NONE);

      //fprintf(screen,"force\n");
  fieldforce_param();

  //fieldtorque_param();

  if (eflag_global){
    double energy_all;
    double volume = domain->xprd * domain->yprd * domain->zprd;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all * volume;

  }

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    double volume = domain->xprd * domain->yprd * domain->zprd;
    for (i = 0; i < 6; i++) virial[i] = volume * virial_all[i]; // DOUBLE CHECK THIS CALCULATION
  }
  if (atom->natoms != natoms_original) {
    natoms_original = atom->natoms;
  }
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

  if (triclinic) domain->lamda2x(atom->nlocal);
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

  int Dim = domain->dimension;
  int (*procneigh)[2] = comm->procneigh;
  
  
  int ntypes = atom->ntypes;
  int ntypecross = ntypes*(ntypes+1)/2;


  // style coeffs
  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;

  //memory->create(potential_type_list,ntypes,ntypes+1,"pppm:potential_type_list"); // lets you know if a type has a potential

  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick");
 
  memory->create4d_offset(density_brick_types,ntypes+1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick_types");
  memory->create5d_offset(gradWtype,ntypecross+1, 0, Dim-1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:gradWtype");

  memory->create(density_fft,nfft_both,"pppm:density_fft");
  memory->create(greensfn,nfft_both,"pppm:greensfn");
  memory->create(worki,2*nfft_both,"pppm:worki");
  memory->create(workj,2*nfft_both,"pppm:workj");
  memory->create(work1,2*nfft_both,"pppm:work1");
  memory->create(work2,2*nfft_both,"pppm:work2");
  memory->create(ktmp,2*nfft_both,"tild:ktmp");
  memory->create(ktmpi,2*nfft_both,"tild:ktmpi");
  memory->create(ktmpj,2*nfft_both,"tild:ktmpj");
  memory->create(ktmp2,2*nfft_both, "tild:ktmp2");
  memory->create(ktmp2i,2*nfft_both, "tild:ktmp2i");
  memory->create(ktmp2j,2*nfft_both, "tild:ktmp2j");
  memory->create(tmp, nfft, "tild:tmp");
  memory->create(vg,ntypecross+1,6,nfft_both,"pppm:vg");
  memory->create(vg_hat,ntypecross+1,6,2*nfft_both,"pppm:vg_hat");
  memory->create(density_fft_types,ntypes+1, nfft_both, "pppm:density_fft_types");
  memory->create(kappa_density, nfft_both, "pppm:kappa_density");
  memory->create(density_of_types_fft_ed,ntypes+1, nfft_both, "pppm:density_of_types_fft_ed");
  memory->create(potent,ntypecross+1,nfft_both,"pppm:potent"); // Voignot 
  memory->create(potent_hat,ntypecross+1,2*nfft_both,"pppm:potent_hat");
  memory->create(grad_potent,ntypecross+1,domain->dimension,nfft_both,"pppm:grad_potent");
  memory->create(grad_potent_hat,ntypecross+1, domain->dimension,2*nfft_both,"pppm:grad_potent_hat");
  // for (int ind = 0; ind < MAX_GROUP; ind++) assigned_pot[ind] = -1;

  // determine if a type has a density function description

  // defines whether atom type has an associated density descripion 


  if (triclinic == 0) {
    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm:fkz");
    memory->create1d_offset(fkx2,nxlo_fft,nxhi_fft,"pppm:fkx2");
    memory->create1d_offset(fky2,nylo_fft,nyhi_fft,"pppm:fky2");
    memory->create1d_offset(fkz2,nzlo_fft,nzhi_fft,"pppm:fkz2");
  } else {
    memory->create(fkx,nfft_both,"pppm:fkx");
    memory->create(fky,nfft_both,"pppm:fky");
    memory->create(fkz,nfft_both,"pppm:fkz");
    memory->create(fkx2,nfft_both,"pppm:fkx2");
    memory->create(fky2,nfft_both,"pppm:fky2");
    memory->create(fkz2,nfft_both,"pppm:fkz2");
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


  if (differentiation_flag == 1)
    cg = new GridComm(lmp,world,ntypes+1,ntypes+1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  else
    cg = new GridComm(lmp,world,3*(ntypes+1),ntypes+1,
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
  int ntypes = atom->ntypes;
  memory->destroy(setflag);

  //memory->destroy(potential_type_list)
  //memory->destroy(potent_type_map);
  //memory->destroy(pot_type_map);

  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset(density_brick_types,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(kappa_density);
  memory->destroy(density_fft_types);
  memory->destroy(density_of_potentials_fft_ed);
  memory->destroy(density_of_types_fft_ed);
  memory->destroy(ktmp);
  memory->destroy(ktmpi);
  memory->destroy(ktmpj);
  memory->destroy(ktmp2);
  memory->destroy(ktmp2i);
  memory->destroy(ktmp2j);
  memory->destroy(tmp);

  memory->destroy(vg);
  memory->destroy(vg_hat);
  memory->destroy(potent);
  memory->destroy(potent_hat);
  memory->destroy(grad_potent);
  memory->destroy(grad_potent_hat);

  memory->destroy5d_offset(gradWtype, 0,
                          nzlo_out,nylo_out,
                          nxlo_out);

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
  memory->destroy(worki);
  memory->destroy(workj);
  memory->destroy(work1);
  memory->destroy(work2);

  if (triclinic == 0) {
    memory->destroy1d_offset(fkx,nxlo_fft);
    memory->destroy1d_offset(fky,nylo_fft);
    memory->destroy1d_offset(fkz,nzlo_fft);
    memory->destroy1d_offset(fkx2,nxlo_fft);
    memory->destroy1d_offset(fky2,nylo_fft);
    memory->destroy1d_offset(fkz2,nzlo_fft);
  } else {
    memory->destroy(fkx);
    memory->destroy(fky);
    memory->destroy(fkz);
    memory->destroy(fkx2);
    memory->destroy(fky2);
    memory->destroy(fkz2);
  }

  memory->destroy(gf_b);
  if (stagger_flag) gf_b = NULL;
  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(drho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);

  memory->destroy(temp);

  delete fft1;
  delete fft2;
  delete remap;
  delete cg;
  fft1 = fft2 = NULL;
  remap = NULL;
  cg = NULL;
  //memory->destroy(assigned_pot);
}

/* ----------------------------------------------------------------------
   set size of FFT grid  and g_ewald_6
   for Dispersion interactions
------------------------------------------------------------------------- */

void TILD::set_grid_6()
{
  // Calculate csum
  if (!gridflag_6) set_n_pppm_6();
  while (!factorable(nx_pppm_6)) nx_pppm_6++;
  while (!factorable(ny_pppm_6)) ny_pppm_6++;
  while (!factorable(nz_pppm_6)) nz_pppm_6++;

}

/* ----------------------------------------------------------------------
   calculate nx_pppm, ny_pppm, nz_pppm for dispersion interaction
   ---------------------------------------------------------------------- */

void TILD::set_n_pppm_6()
{
  bigint natoms = atom->natoms;

  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;
  double h, h_x,h_y,h_z;

  double acc_kspace = accuracy;
  if (accuracy_kspace_6 > 0.0) acc_kspace = accuracy_kspace_6;

  // initial value for the grid spacing
  h = h_x = h_y = h_z = 4.0/g_ewald_6;
  // decrease grid spacing untill required precision is obtained
  int count = 0;
  while(1) {

    // set grid dimension
    nx_pppm_6 = static_cast<int> (xprd/h_x);
    ny_pppm_6 = static_cast<int> (yprd/h_y);
    nz_pppm_6 = static_cast<int> (zprd_slab/h_z);

    if (nx_pppm_6 <= 1) nx_pppm_6 = 2;
    if (ny_pppm_6 <= 1) ny_pppm_6 = 2;
    if (nz_pppm_6 <= 1) nz_pppm_6 = 2;

    //set local grid dimension
    int npey_fft,npez_fft;
    if (nz_pppm_6 >= nprocs) {
      npey_fft = 1;
      npez_fft = nprocs;
    } else procs2grid2d(nprocs,ny_pppm_6,nz_pppm_6,&npey_fft,&npez_fft);

    int me_y = me % npey_fft;
    int me_z = me / npey_fft;

    nxlo_fft_6 = 0;
    nxhi_fft_6 = nx_pppm_6 - 1;
    nylo_fft_6 = me_y*ny_pppm_6/npey_fft;
    nyhi_fft_6 = (me_y+1)*ny_pppm_6/npey_fft - 1;
    nzlo_fft_6 = me_z*nz_pppm_6/npez_fft;
    nzhi_fft_6 = (me_z+1)*nz_pppm_6/npez_fft - 1;

    double qopt = compute_qopt_6();

    double df_kspace = sqrt(qopt/natoms)*csum/(xprd*yprd*zprd_slab);

    count++;

    // break loop if the accuracy has been reached or too many loops have been performed
    if (df_kspace <= acc_kspace) break;
    if (count > 500) error->all(FLERR, "Could not compute grid size for Dispersion");
    h *= 0.95;
    h_x = h_y = h_z = h;
  }
}

/* ----------------------------------------------------------------------
   Compute qopt for Dispersion interactions
------------------------------------------------------------------------- */

double TILD::compute_qopt_6()
{
  double qopt;
  if (differentiation_flag == 1) {
    qopt = compute_qopt_6_ad();
  } else {
    qopt = compute_qopt_6_ik();
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   Compute qopt for the ik differentiation scheme and Dispersion interaction
------------------------------------------------------------------------- */

double TILD::compute_qopt_6_ik()
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
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(MY_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz +
                       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew);
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*MY_PI*MY_PI*MY_PI/9.0 * dot2;
                sum2 += -u2*term*MY_PI*rtpi/3.0*dot1;
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
   Compute qopt for the ad differentiation scheme and Dispersion interaction
------------------------------------------------------------------------- */

double TILD::compute_qopt_6_ad()
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
  double sum1,sum2,sum3,sum4;
  double dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1/inv2ew;
  double rtpi = sqrt(MY_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);

        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) +
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {

          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz +
                       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew);
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*MY_PI*MY_PI*MY_PI/9.0 * dot2;
                sum2 += -term*MY_PI*rtpi/3.0 * u2 * dot2;
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
   set size of FFT grid (nx,ny,nz_pppm) and g_ewald
   for Coulomb interactions
------------------------------------------------------------------------- */

void TILD::set_grid()
{
  double q2 = qsqsum * force->qqrd2e;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
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

  double h, h_x,h_y,h_z;
  bigint natoms = atom->natoms;

  // if (!gewaldflag) {
  //   g_ewald = accuracy*sqrt(natoms*cutoff*xprd*yprd*zprd) / (2.0*q2);
  //   if (g_ewald >= 1.0)
  //     error->all(FLERR,"KSpace accuracy too large to estimate G vector");
  //   g_ewald = sqrt(-log(g_ewald)) / cutoff;
  // }

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {
    h = h_x = h_y = h_z = grid_res;
    int count = 0;

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

    }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;

  char str[128];
  if (comm->me == 0) {
    if (screen){
      fprintf(screen,"  grid points x dim: %d\n",nx_pppm);
      fprintf(screen,"  grid points y dim: %d\n",ny_pppm);
      fprintf(screen,"  grid points z dim: %d\n",nz_pppm);
}
    if (logfile){
      fprintf(logfile,"  grid points x dim: %d\n",nx_pppm);
      fprintf(logfile,"  grid points y dim: %d\n",ny_pppm);
      fprintf(logfile,"  grid points z dim: %d\n",nz_pppm);
    }
  }
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

      h = h_x = h_y = h_z = grid_res;
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

        // double qopt = compute_qopt();

        // double df_kspace = compute_df_kspace();

        count++;

        // break loop if the accuracy has been reached or
        // too many loops have been performed

        // if (df_kspace <= accuracy) break;
        if (count > 500) error->all(FLERR, "Could not compute grid size");
        h *= 0.95;
        h_x = h_y = h_z = h;
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
    error->all(FLERR," grid is too large");
    
  char str[128];
  if (comm->me == 0) {
    if (screen){
      fprintf(screen,"  grid points x dim: %d\n",nx_pppm);
      fprintf(screen,"  grid points y dim: %d\n",ny_pppm);
      fprintf(screen,"  grid points z dim: %d\n",nz_pppm);
}
    if (logfile){
      fprintf(logfile,"  grid points x dim: %d\n",nx_pppm);
      fprintf(logfile,"  grid points y dim: %d\n",ny_pppm);
      fprintf(logfile,"  grid points z dim: %d\n",nz_pppm);
    }
  }
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void TILD::allocate_peratom()
{

  int (*procneigh)[2] = comm->procneigh;

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
  // memory->create5d_offset(gradWtype,group->ngroup, 0, Dim,
  //                         nzlo_out,nzhi_out,nylo_out,nyhi_out,
  //                         nxlo_out,nxhi_out,"tild:gradWtype");
  // memory->create(gradWtype, group->ngroup, Dim, , "tild:gradWtype");

  // create ghost grid object for rho and electric field communication

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
  cg_peratom = NULL;
}

void TILD::init_cross_potentials(){
  
  // Represents the 0 to N-1 points that are captured on this grid.
  int Dim = domain->dimension;
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  double mdr2;

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
  int ntypes = atom->ntypes;
  double scale_inv = 1.0/ nx_pppm/ ny_pppm/ nz_pppm;
  n = 0;


  double vole = domain->xprd * domain->yprd * domain->zprd; // Note: the factor of V comes from the FFT
  // Loop over potental styles
  for (int itype = 1; itype <= ntypes; itype++) {
    // Skip if type cross-interaction does not use density/tild

    for (int jtype = itype; jtype <= ntypes; jtype++) {
      // Skip if type cross-interaction does not use density/tild
      if ( potent_type_map[0][itype][jtype] == 1) continue;
      fprintf(screen,"i j %d %d\n", itype, jtype);
      int loc = itype*(jtype+1)/2;


      // If both parameters are Gaussian, just do analytical convolution
      if (potent_type_map[1][itype][jtype] == 1 or (potent_type_map[1][itype][itype] == 1 && potent_type_map[1][jtype][jtype] == 1) ){
        // mixing
        double a2_mix;
        if (mix_flag == 1) {
          a2_mix = a2[itype][itype] + a2[jtype][jtype];
        } else {
          a2_mix = a2[itype][jtype] + a2[itype][jtype];
        }
        const double* p = &a2_mix;
        init_potential(potent[loc], 1, p);

        int j = 0;
        for (int i = 0; i < nfft; i++) {
          ktmp[j++] = potent[loc][i];
          ktmp[j++] = ZEROF;
        }
  
        fft1->compute(ktmp, potent_hat[loc], 1);

        for (int i = 0; i < 2 * nfft; i++) {
          potent_hat[loc][i] *= scale_inv;
        }
/*
        init_potential_ft(potent_hat[loc], 1, p);
*/
      } 
      // Computational Convolution
      else {

        // calculate 1st and 2nd potentials to convolv
        int style;
      fprintf(screen,"i j %d %d %d\n", itype, jtype, mix_flag);
        if (mix_flag == 1) {
          calc_work(work1, itype, itype);
          if ( itype == jtype) {
            for (int i = 0; i < 2*nfft; i++) work2[i] = work1[i];
          } else {
            calc_work(work2, jtype, jtype);
          }
        } else {
          calc_work(work1, itype, jtype);
          for (int i = 0; i < 2*nfft; i++) work2[i] = work1[i];
        }

        // Convolution of potentials
        n = 0;
        for (int i = 0; i < nfft; i++) {
          complex_multiply(work1, work2, potent_hat[loc], n);
          n += 2;
        }

        fft2->compute(potent_hat[loc], ktmp, -1);

        n = 0;
        for (int j = 0; j < nfft; j++){
          potent[loc][j] = ktmp[n];
          n += 2;
        }
      }

      get_k_alias(potent_hat[loc], grad_potent_hat[loc]);
      for (int i=0; i < Dim; i ++){
        fft2->compute(grad_potent_hat[loc][i], work2, -1);
        n = 0;
        for (int j = 0; j < nfft; j++){
          grad_potent[loc][i][j] = -work2[n];
          n+=2;
        }
      } 

      // output
      std::string fnameU = "U_lammps_"+std::to_string(itype)+"-"+std::to_string(jtype)+".txt";
      std::string fnamegradU = "gradU_lammps_"+std::to_string(itype)+"-"+std::to_string(jtype)+".txt";
      std::string fnamegradUhat = "gradUhat_lammps_"+std::to_string(itype)+"-"+std::to_string(jtype)+".txt";
      std::string fnamegradUhatI = "gradUhatimag_lammps_"+std::to_string(itype)+"-"+std::to_string(jtype)+".txt";
      ofstream fileU(fnameU);
      ofstream filegradU(fnamegradU);
      ofstream filegradUhat(fnamegradUhat);
      ofstream filegradUhatI(fnamegradUhatI);
      n=0;
      for (int j = 0; j < nfft; j++) {
          fileU << j << '\t' << potent[loc][j] << std::endl;
          filegradU << j << '\t' << grad_potent[loc][0][j] << '\t' << grad_potent[loc][1][j] << '\t' << grad_potent[loc][2][j] << std::endl;
          filegradUhat << j << '\t' << grad_potent_hat[loc][0][n] << '\t' << grad_potent_hat[loc][1][n] << '\t' << grad_potent_hat[loc][2][n] << std::endl;
          filegradUhatI << j << '\t' << grad_potent_hat[loc][0][n+1] << '\t' << grad_potent_hat[loc][1][n+1] << '\t' << grad_potent_hat[loc][2][n+1] << std::endl;
          n += 2;
      }

      fileU.close();
      filegradU.close();
      filegradUhat.close();
      filegradUhatI.close();
    }
  }

}

int TILD::get_style( const int i, const int j) {
  for (int istyle = 1; istyle <= nstyles; istyle++) { 
    if ( potent_type_map[istyle][i][j] == 1 ) return istyle;
  }
  return 0;
}

void TILD::calc_work(double* work, const int itype, const int jtype){
  double scale_inv = 1.0/ nx_pppm/ ny_pppm/ nz_pppm;

  // needs work is this right for cross terms of the same potential?
  
  double params[4];
  for (int i = 0; i < 4; i++) params[i] = 0.0;

  int style = get_style(itype, jtype);
  if (style == 1) {
    params[0] = a2[itype][jtype];
    init_potential_ft(work, style, params);
  } else if (style == 2) {
    params[0] = rp[itype][jtype];
    params[1] = xi[itype][jtype];
// need another check to return the analytical value for gaussian FFT
    init_potential(tmp, style, params);

    int j = 0;
    for (int i = 0; i < nfft; i++) {
      work[j++] = tmp[i];
      work[j++] = ZEROF;
    }

    fft1->compute(work, work, 1);

    for (int i = 0; i < 2 * nfft; i++) {
      work[i] *= scale_inv;
    }
  }
}

// analytical fourier transform of potential
void TILD::init_potential_ft(FFT_SCALAR *wk1, const int type, const double* parameters){

  int n = 0;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  double zper,yper,xper;
  double k2;
  double factor = 4.0 * MY_PI * MY_PI;

  double vole = domain->xprd * domain->yprd * domain->zprd; // Note: the factor of V comes from the FFT
  int Dim = domain->dimension;

  if (type == 1) {
                                                                // should be 3/2 right?
    for (int z = nzlo_fft; z <= nzhi_fft; z++) {
      zper = static_cast<double>(z) / zprd;
      if (z >= nz_pppm / 2.0) {
        zper -= static_cast<double>(nz_pppm) / zprd;
      }
      double zper2 = factor * zper * zper; 
      //cout << 'z' << '\t' << z << '\t' << zprd << '\t' << nz_pppm << '\t' << zper2 << endl;

      for (int y = nylo_fft; y <= nyhi_fft; y++) {
        yper = static_cast<double>(y) / yprd;
        if (y >= ny_pppm / 2.0) {
          yper -= static_cast<double>(ny_pppm) / yprd;
        }
        double yper2 = factor * yper * yper; 
      //cout << 'y' << '\t' << y << '\t' << yprd << '\t' << ny_pppm << '\t' << yper2 << endl;

        for (int x = nxlo_fft; x <= nxhi_fft; x++) {
          xper = static_cast<double>(x) / xprd;
          if (x >= nx_pppm / 2.0) {
            xper -= static_cast<double>(nx_pppm) / xprd;
          }
      //cout << 'x' << '\t' << x << '\t' << xprd << '\t' << nx_pppm << '\t' << (factor * xper * xper) << endl;

          k2 = (factor * xper * xper) + yper2 + zper2;
          cout << n << '\t' << k2 << '\t' << exp( -k2 * 0.5 * parameters[0]) << endl;
          wk1[n++] = exp(-k2 * 0.5 * parameters[0]);
          wk1[n++] = ZEROF;
        
        }
      }
    }
  }
}

void TILD::init_potential(FFT_SCALAR *wk1, const int type, const double* parameters){

  int m,l,k;
  int n = 0;
  double xprd=domain->xprd;
  double yprd=domain->yprd;
  double zprd=domain->zprd;
  double zper,yper,xper;
  double mdr2;

  double vole = domain->xprd * domain->yprd * domain->zprd; // Note: the factor of V comes from the FFT
  int Dim = domain->dimension;

  if (type == 1) {
                                                                // should be 3/2 right?
    double pref = vole / (pow( sqrt(2.0 * PI * (parameters[0]) ), Dim));
    for (m = nzlo_fft; m <= nzhi_fft; m++) {
      zper = zprd * (static_cast<double>(m) / nz_pppm);
      if (zper >= zprd / 2.0) {
        zper = zprd - zper;
      }

      for (l = nylo_fft; l <= nyhi_fft; l++) {
        yper = yprd * (static_cast<double>(l) / ny_pppm);
        if (yper >= yprd / 2.0) {
          yper = yprd - yper;
        }

        for (k = nxlo_fft; k <= nxhi_fft; k++) {
          xper = xprd * (static_cast<double>(k) / nx_pppm);
          if (xper >= xprd / 2.0) {
            xper = xprd - xper;
          }

          mdr2 = xper * xper + yper * yper + zper * zper;
          wk1[n++] = exp(-mdr2 * 0.5 / parameters[0]) * pref;
        
        }
      }
    }
  }
  else if (type == 2) {
    for (m = nzlo_fft; m <= nzhi_fft; m++) {
      zper = zprd * (static_cast<double>(m) / nz_pppm);
      if (zper >= zprd / 2.0) {
        zper = zprd - zper;
      }

      for (l = nylo_fft; l <= nyhi_fft; l++) {
        yper = yprd * (static_cast<double>(l) / ny_pppm);
        if (yper >= yprd / 2.0) {
          yper = yprd - yper;
        }

        for (k = nxlo_fft; k <= nxhi_fft; k++) {
          xper = xprd * (static_cast<double>(k) / nx_pppm);
          if (xper >= xprd / 2.0) {
            xper = xprd - xper;
          }

          mdr2 = xper * xper + yper * yper + zper * zper;
                                    // missing pararentheses, xi should be inside erf()
          wk1[n++] = rho0 * 0.5 * (1.0 - erf((sqrt(mdr2) - parameters[0])/parameters[1])) * vole;
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

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

void TILD::get_k_alias(FFT_SCALAR* wk1, FFT_SCALAR **out){
  int Dim = domain->dimension; 
  double k[Dim];
  int x, y, z;
  int n=0;
  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  // id2 = unstack_stack( id ) ;
  // unstack(id2, n);
  for (z = nzlo_fft; z <= nzhi_fft; z++) {
    if (nz_pppm % 2 == 0 && z == nz_pppm / 2)
      k[2] = 0.0;
    else if (double(z) < double(nz_pppm) / 2.0)
      k[2] = 2 * PI * double(z) / zprd;
    else
      k[2] = 2 * PI * double(z - nz_pppm) / zprd;

    for (y = nylo_fft; y <= nyhi_fft; y++) {
      if (ny_pppm % 2 == 0 && y == ny_pppm / 2)
        k[1] = 0.0;
      else if (double(y) < double(ny_pppm) / 2.0)
        k[1] = 2 * PI * double(y) / yprd;
      else
        k[1] = 2 * PI * double(y - ny_pppm) / yprd;

      for (x = nxlo_fft; x <= nxhi_fft; x++) {
        if (nx_pppm % 2 == 0 && x == nx_pppm / 2)
          k[0] = 0.0;
        else if (double(x) < double(nx_pppm) / 2.0)
          k[0] = 2 * PI * double(x) / xprd;
        else
          k[0] = 2 * PI * double(x - nx_pppm) / xprd;

/*
        out[0][n] = 0.5*(fkx[x]-fkx2[x])*wk1[n+1];
        out[0][n+1] = -0.5*(fkx[x]-fkx2[x])*wk1[n];
        out[1][n] = 0.5*(fky[y]-fky2[y])*wk1[n+1];
        out[1][n+1] = -0.5*(fky[y]-fky2[y])*wk1[n];
        out[2][n] = 0.5*(fkz[z]-fkz2[z])*wk1[n+1];
        out[2][n+1] = -0.5*(fkz[z]-fkz2[z])*wk1[n];
*/
        out[0][n] = -wk1[n + 1] * k[0];
        out[0][n + 1] = wk1[n] * k[0];
        out[1][n] = -wk1[n + 1] * k[1];
        out[1][n + 1] = wk1[n] * k[1];
        out[2][n] = -wk1[n + 1] * k[2];
        out[2][n + 1] = wk1[n] * k[2];
/*
        out[0][n] = 0.0;
        out[0][n + 1] = wk1[n+1] * k[0];
        out[0][n] = 0.0;
        out[1][n + 1] = wk1[n+1] * k[1];
        out[0][n] = 0.0;
        out[2][n + 1] = wk1[n+1] * k[2];
*/
        n += 2;
          
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
        nz+nlow < nzlo || nz+nup > nzhi){
      flag = 1;
      std::cout << i << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n"; 
      
      }
  }

 if (flag) error->one(FLERR,"Out of range atoms - cannot compute TILD");
}


int TILD::modify_param(int narg, char** arg)
{
  int i;
  int ntypes = atom->ntypes;
  if (strcmp(arg[0], "tild/chi") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, "TILD command before simulation box is defined");

      if (narg != 4) error->all(FLERR, "Illegal kspace_modify tild command");

      int ilo,ihi,jlo,jhi;
      force->bounds(FLERR,arg[1],ntypes,ilo,ihi);
      force->bounds(FLERR,arg[2],ntypes,jlo,jhi);
      double chi_one = force->numeric(FLERR,arg[3]);
      int count = 0;
      for (int i = ilo; i <= ihi; i++) {
        for (int j = MAX(jlo,i); j <= jhi; j++) {
          chi[i][j] = chi_one;
          count++;
        }
      }
 
  } else if (strcmp(arg[0], "tild/kappa") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal kspace_modify tild command");
      kappa = atof(arg[1]);
 
  } else if (strcmp(arg[0], "tild/coeff") == 0) {
      if (narg < 3) error->all(FLERR, "Illegal kspace_modify tild command");

      int ilo,ihi,jlo,jhi;
      force->bounds(FLERR,arg[1],ntypes,ilo,ihi);
      force->bounds(FLERR,arg[2],ntypes,jlo,jhi);
      if (strcmp(arg[3], "gaussian") == 0) {
        if (narg < 4) error->all(FLERR, "Illegal kspace_modify tild command");

        double a2_one = force->numeric(FLERR,arg[4])*force->numeric(FLERR,arg[4]);
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo,i); j <= jhi; j++) {
            potent_type_map[1][i][j] = 1;
            potent_type_map[0][i][j] = 0;
            fprintf(screen,"types %d %d %d %d %d\n", i,j, potent_type_map[0][i][j], potent_type_map[1][i][j], ntypes);
            a2[i][j] = a2_one;
          }
        //nstyles++;
        }

      } else if (strcmp(arg[3], "erfc") == 0) {
        if (narg < 5) error->all(FLERR, "Illegal kspace_modify tild command");
        double rp_one = force->numeric(FLERR,arg[4]);
        double xi_one = force->numeric(FLERR,arg[5]);
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo,i); j <= jhi; j++) {
            potent_type_map[2][i][j] = 1;
            potent_type_map[0][i][j] = 0;
            fprintf(screen,"types %d %d %d %d %d\n", i,j, potent_type_map[0][i][j], potent_type_map[1][i][j], potent_type_map[2][i][j]);
            rp[i][j] = rp_one;
            xi[i][j] = xi_one;
          }
        //nstyles++; // eventually will be part of kspace_style hybrid
        }
      } else if (strcmp(arg[3], "none") == 0) {
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo,i); j <= jhi; j++) {
            potent_type_map[0][i][j] = 1;
            for (int istyle = 1; istyle <= nstyles; istyle++) 
              potent_type_map[istyle][i][j] = 0;
          }
        }
      } else error->all(FLERR, "Illegal kspace_modify tild/coeff density function argument");

  } else if (strcmp(arg[0], "mix") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      mix_flag = 1;
      if (strcmp(arg[1], "convolution") == 0)  mix_flag = 1;
      else if (strcmp(arg[1], "define") == 0) mix_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify tild subtract_rho0 argument");
  } else if (strcmp(arg[0], "set_rho0") == 0) {
      if (narg < 2 ) error->all(FLERR, "Illegal kspace_modify tild command");
      if (strcmp(arg[1], "yes") == 0) {
        set_rho0_flag = 1;
        set_rho0 = force->numeric(FLERR,arg[2]);
      } 
      else error->all(FLERR, "Illegal kspace_modify tild subtract_rho0 argument");
  } else if (strcmp(arg[0], "subtract_rho0") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      if (strcmp(arg[1], "yes") == 0) sub_flag = 1;
      else if (strcmp(arg[1], "no") == 0) sub_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify tild subtract_rho0 argument");
  } else if (strcmp(arg[0], "normalize_by_rho0") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      if (strcmp(arg[1], "yes") == 0) norm_flag = 1;
      else if (strcmp(arg[1], "no") == 0) norm_flag = 0;
      else 
        error->all(FLERR, "Illegal kspace_modify tild normalize_by_rho0 argument");

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

void TILD::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  int Dim = domain->dimension;

  if (flag == FORWARD_NONE){
    for (int k = 0; k <= atom->ntypes; k++) {
    for (int j = 0; j < Dim; j++) {
      FFT_SCALAR *src = &gradWtype[k][j][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        buf[n++] = src[list[i]];
    }
    }
  }
}

void TILD::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  int Dim = domain->dimension;

  if (flag == FORWARD_NONE){
    for (int k = 0; k <= atom->ntypes; k++) {
    for (int j = 0; j < Dim; j++) {
      FFT_SCALAR *dest = &gradWtype[k][j][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        dest[list[i]] = buf[n++];
      }
    }
  }
}

void TILD::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list) {
  int n = 0;
      //fprintf(screen,"pack_reverse\n");
  if (flag == REVERSE_RHO_NONE) {
    for (int k = 0; k <= atom->ntypes; k++) {
      FFT_SCALAR *src = &density_brick_types[k][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        buf[n++] = src[list[i]];
    }
  }
}

void TILD::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
      //fprintf(screen,"unpack_reverse\n");
  if (flag == REVERSE_RHO_NONE) {
    for (int k = 0; k <= atom->ntypes; k++) {
      FFT_SCALAR *dest = &density_brick_types[k][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        dest[list[i]] += buf[n++];
    }
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
  double cuthalf = 0.5*neighbor->skin; // removed qdist no offset needed
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


void TILD::make_rho_none()
{
  int k,l,m,n,nx,ny,nz,mx,my,mz;
  int ntypes = atom->ntypes;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  // clear 3d density array
  // for (k = 0; k < nsplit_alloc; k++)
  //   memset(&(density_brick_none[k][nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
  //          ngrid_6*sizeof(FFT_SCALAR));

  for (int k = 0; k <= ntypes; k++) {
    memset(&(density_brick_types[k][nzlo_out][nylo_out][nxlo_out]),0,
           ngrid*sizeof(FFT_SCALAR));
  }


  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double *mass = atom->mass;
  double temp_mass = 0;

  for (int i = 0; i < nlocal; i++) {

    temp_mass = mass[type[i]];
    if (temp_mass == 0) continue;
    // Skip if the particle type has no TILD interaction
    if (potent_type_map[0][type[i]][type[i]] == 1) continue;

    // do the following for all 4 grids
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;
    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);
    z0 = delvolinv;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          w = x0*rho1d[0][l];
          if (mask[i]) {
            density_brick_types[type[i]][mz][my][mx] += w * temp_mass;
            density_brick_types[0][mz][my][mx] += density_brick_types[type[i]][mz][my][mx];
          }
        }
      }
    }
  }
}

void TILD::brick2fft_none()
{
  int k,n,ix,iy,iz;
  int ntypes = atom->ntypes;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  //  std::ofstream rhoA("rhot_lammps.txt");
  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++){
        density_fft_types[0][n++] = density_brick_types[0][iz][iy][ix];
        for (int k = 1; k <= ntypes; k++) {
          if ( potent_type_map[0][k][k] == 1) continue;
          density_fft_types[k][n++] = density_brick_types[k][iz][iy][ix];
        }
          // rhoA<<ix<<'\t'<<iy<<'\t'<<iz<<'\t'<<density_brick_types[k][iz][iy][ix] <<std::endl;;
      }

  remap->perform(density_fft_types[0],density_fft_types[0],work1);
  for (int k = 1; k <= ntypes; k++) {
    if ( potent_type_map[0][k][k] == 1) continue;
    remap->perform(density_fft_types[k],density_fft_types[k],work1);
  }
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

void TILD::accumulate_gradient() {
  int Dim = domain->dimension;
  double temp_param;
  double scale_inv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);
  int type1 = 0, type2 = 0;
  int n = 0;
  int itype, jtype;
  int ntypes = atom->ntypes;

  for ( int ktype = 1; ktype <= ntypes; ktype++)
    for (int i = 0; i < Dim; i++)
      memset(&(gradWtype[ktype][i][nzlo_out][nylo_out][nxlo_out]), 0,
             ngrid * sizeof(FFT_SCALAR));

  //typedef std::tuple<int, int, double, int, int> tup_iidii;

  // This part is for the specified chi interactions.
  // Kappa (incompressibility) is later

  double tmp_kappa = kappa;

  FFT_SCALAR tmp_sub = 0.0;
  if (subtract_rho0 == 1)
    tmp_sub = rho0;

  if ( tmp_kappa != 0 ) {
    n = 0;
    for (int k = 0; k < nfft; k++) {
      kappa_density[k] = density_fft_types[0][k] - tmp_sub;
      work1[n++] = kappa_density[k];
      work1[n++] = ZEROF;
    }
    fft1->compute(work1, work1, 1);
    for (int k = 0; k < 2*nfft; k++) {
      work1[k] *= scale_inv;
    }
  }
  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      if ( potent_type_map[0][itype][jtype] == 1) continue;
      //fprintf(screen,"i j %d %d\n", itype,jtype);

      int loc = itype*(jtype+1)/2;
      double tmp_chi = chi[itype][jtype];

      if ( tmp_chi != 0 && tmp_kappa != 0 ) continue;

      if (normalize_by_rho0 == 1) {
        tmp_chi /= rho0;
        tmp_kappa = kappa/rho0;
        if (itype == jtype) {
          tmp_kappa /= 2.0;
        }
      }

      if ( tmp_chi != 0) {
        n = 0;
        for (int k = 0; k < nfft; k++) {
          worki[n] = density_fft_types[itype][k];
          workj[n] = density_fft_types[jtype][k];
          n++;
          worki[n] = ZEROF;
          workj[n] = ZEROF;
          n++;
        }
        fft1->compute(worki, worki, 1);
        fft1->compute(workj, workj, 1);
        for (int k = 0; k < 2*nfft; k++) {
          worki[k] *= scale_inv;
          workj[k] *= scale_inv;
        }
      }

      //fprintf(screen,"pre ev_calculation\n");
      if (eflag_global || vflag_global) {
        ev_calculation(kappa_density, work1, worki, itype, jtype);
      }
      //fprintf(screen,"post ev_calculation\n");

      for (int i = 0; i < Dim; i++) {
/*
        n = 0;
        for (int k = 0; k < nfft; k++) {
          work2[n] = grad_potent[loc][i][n];
          //work2[n] = grad_potent_hat[loc][i][n];
          n++;
          work2[n] = grad_potent[loc][i][n];
          //work2[n] = grad_potent_hat[loc][i][n];
          n++;
        }
      //fprintf(screen,"complex mult\n");
*/

        n = 0;
        for (int k = 0; k < nfft; k++) {
/*
          if ( tmp_kappa != 0 ) {
            complex_multiply(grad_potent_hat[loc][i], work1, ktmp2, n);
            //complex_multiply(work1, work2, ktmp2, n);
          }
          if ( tmp_chi != 0) {
            //complex_multiply(worki, work2, ktmp2i, n);
            //complex_multiply(workj, work2, ktmp2j, n);
          }
*/
          complex_multiply(grad_potent_hat[loc][i], worki, ktmp2i, n);
          complex_multiply(grad_potent_hat[loc][i], workj, ktmp2j, n);
          n += 2;
        }
      //fprintf(screen,"compute\n");
        if ( tmp_kappa != 0 ) {
          fft2->compute(ktmp2, ktmp, -1);
        }
        if ( tmp_chi != 0) {
          fft2->compute(ktmp2i, ktmpi, -1);
          fft2->compute(ktmp2j, ktmpj, -1);
        }

        n = 0;
        for (int k = nzlo_in; k <= nzhi_in; k++)
          for (int m = nylo_in; m <= nyhi_in; m++)
            for (int o = nxlo_in; o <= nxhi_in; o++) {
/*
              if ( tmp_kappa != 0 ) {
                gradWtype[loc][i][k][m][o] += ktmp[n] * tmp_kappa;
              }
              if ( tmp_chi != 0) {
                gradWtype[loc][i][k][m][o] += ktmpi[n] * tmp_chi;
                gradWtype[loc][i][k][m][o] += ktmpj[n] * tmp_chi;
              }
*/
              gradWtype[loc][i][k][m][o] += ktmpi[n] * (tmp_chi + tmp_kappa);
              gradWtype[loc][i][k][m][o] += ktmpj[n] * (tmp_chi + tmp_kappa);
              n += 2;
            }
      }
    }
  }
}

void TILD::fieldforce_param(){
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;
  int ngroups = group -> ngroup;
  int dim=domain->dimension;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;

  int nlocal = atom->nlocal;

  // Convert field to force per particle
  ofstream delvolinv_forces("delv_forces_lammps.txt");
  // ofstream grid_vol_forces("gridvol_forces_lammps.txt");
  // ofstream both_forces("delv_gridvol_forces_lammps.txt");
  int *type = atom->type;
  double *mass = atom->mass;
  double temp_mass = 0;

  for (int i = 0; i < nlocal; i++) {

    temp_mass = mass[type[i]];
    if (temp_mass == 0) continue;
    int temp_type = type[i];
    if (temp_type == -1) continue;
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

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
          // gradWpotential 0 to 1 ?????
          ekx -= x0 * gradWtype[temp_type][0][mz][my][mx];
          eky -= x0 * gradWtype[temp_type][1][mz][my][mx];
          ekz -= x0 * gradWtype[temp_type][2][mz][my][mx];
        }
      }
    }

    // convert field to force
    
    f[i][0] += ekx;
    f[i][1] += eky;
    f[i][2] += ekz;
    // f[i][0] += grid_vol * ekx;
    // f[i][1] += grid_vol * eky;
    // f[i][2] += grid_vol * ekz;
    delvolinv_forces << i <<'\t'<<  f[i][0] <<'\t'<< f[i][1] <<'\t'<<  f[i][2] << endl;
    //delvolinv_forces << i <<'\t'<<  ekx <<'\t'<< eky <<'\t'<<  ekz << endl;
    // grid_vol_forces<<i<<'\t'<< "0 " <<'\t'<<  grid_vol*ekx << endl;
    // grid_vol_forces<<i<<'\t'<< "1 " <<'\t'<<  grid_vol*eky << endl;
    // grid_vol_forces<<i<<'\t'<< "2 " <<'\t'<<  grid_vol*ekz << endl;
    // both_forces<<i<<'\t'<< "0 " <<'\t'<<  delvolinv*grid_vol*ekx << endl;
    // both_forces<<i<<'\t'<< "1 " <<'\t'<<  delvolinv*grid_vol*eky << endl;
    // both_forces<<i<<'\t'<< "2 " <<'\t'<<  delvolinv*grid_vol*ekz << endl;
    // forces<< "Y " <<'\t'<< delvolinv*eky << "\t";
    // forces<< "Z " <<'\t'<< delvolinv*ekz << std::endl;
  }
delvolinv_forces.close();
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
void TILD::compute_rho_coeff(FFT_SCALAR **coeff , FFT_SCALAR **dcoeff,
                                 int ord)
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;
  memory->create2d_offset(a,ord,-ord,ord,"pppm/disp:a");

  for (k = -ord; k <= ord; k++)
    for (l = 0; l < ord; l++)
      a[l][k] = 0.0;

  a[0][0] = 1.0;
  for (j = 1; j < ord; j++) {
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

  m = (1-ord)/2;
  for (k = -(ord-1); k < ord; k += 2) {
    for (l = 0; l < ord; l++)
      coeff[l][m] = a[l][k];
    for (l = 1; l < ord; l++)
      dcoeff[l-1][m] = l*a[l][k];
    m++;
  }

  memory->destroy2d_offset(a,-ord);
}


inline void TILD::complex_multiply(FFT_SCALAR *in1,FFT_SCALAR  *in2,FFT_SCALAR  *out, int n){
  out[n] = ((in1[n] * in2[n]) - (in1[n + 1] * in2[n + 1]));
  out[n + 1] = ((in1[n + 1] * in2[n]) + (in1[n] * in2[n + 1]));
}

inline void TILD::complex_multiply(double *in1,double  *in2, int n){
  double temp = (in1[n] * in2[n] - in1[n + 1] * in2[n + 1]);
  in2[n + 1] = (in1[n + 1] * in2[n] + in1[n] * in2[n + 1]);
  in2[n] = temp;
}

//void TILD::ev_calculation() {
void TILD::ev_calculation(FFT_SCALAR *rho1, double *wk1, double *wki, const int itype, const int jtype) {
  int n = 0;
  double eng;
  double vtmp;
  double scale_inv = 1.0/(nx_pppm *ny_pppm * nz_pppm);
  double s2 = scale_inv * scale_inv;
  int ntypes = atom->ntypes;
  FFT_SCALAR *dummy, *wk2;

  double tmp_rho_div = 1.0;
  if (normalize_by_rho0 == 1) tmp_rho_div = rho0*rho0;
  // Convolve kspace-potential with fft(den_group)
  double tmp_chi = chi[itype][jtype];
/*
  if ( tmp_chi != 0) {
    fft1->compute(wki,wki,1);
    
    for (int k = 0; k < 2*nfft; k++) {
      wki[n++] *= scale_inv;
    }
  }
*/
  // take jtype density into k-space
  
  // convolve itype-jtype interaction potential and itype density
  int loc = itype*(jtype+1)/2;
  n = 0;
  //fft1->compute(potent[loc], ktmp, 1);
  for (int k = 0; k < nfft; k++) {
      complex_multiply(potent_hat[loc], wki, ktmp2i, n);
/*
    if ( tmp_chi != 0) {
    }
    if ( kappa != 0 ) {
      complex_multiply(potent_hat[loc], wk1, ktmp2, n);
    }
*/
    n += 2;
  }
  // IFFT the convolution
/*
  if ( kappa != 0 ) {
    fft2->compute(ktmp2, ktmp2, -1);
  }
  
  if ( tmp_chi != 0) {
    // not sure about this
  }
*/
    fft2->compute(ktmp2i, ktmp2i, -1);

  double tmp_kappa = kappa;
  if (itype != jtype) {
    tmp_kappa *= 2.0;
  }
  double factor = 1.0 / delvolinv / tmp_rho_div / 2.0;
  //double factor = 1.0 / delvolinv * same_group_factor * 0.5 / tmp_rho_div;
  //double factor = (chi[itype][jtype] + kappa) / delvolinv * same_group_factor * 0.5 / tmp_rho_div;

  // chi and kappa energy contribution
  //fprintf(screen,"energy %f %f %f %f %f\n", energy, delvolinv, same_group_factor, 0.5, tmp_rho_div);
  if (eflag_global) {
    n = 0;
    eng = 0.0;
    for (int k = 0; k < nfft; k++) {
      // rho_i * u_ij * rho_j * chi * prefactor
/*
      if ( tmp_chi != 0) energy += tmp_chi * ktmp2i[n] * density_fft_types[jtype][k] * factor;
      // rho_i * u_ij * rho_j * prefactor * kappa
      if ( kappa != 0 ) energy += kappa * ktmp2[n] * rho1[k] * factor;
*/
      eng += ktmp2i[n] * density_fft_types[jtype][k];
      //engk += ktmp2i[n] * rho0
      n += 2;
    }
    energy += eng * factor * (tmp_chi + tmp_kappa); // - ( tmp_kappa *engk
  }
  //fprintf(screen,"energy %f\n", energy);

  // pressure tensor calculation
  if (vflag_global) {
    // loop over stress tensor
   //fprintf(screen,"virial %f %f %f %f %f %f\n", virial[0], virial[1], virial[2], virial[3], virial[4], virial[5]);
    for (int i = 0; i < 6; i++) {
      double factor2 = factor * ( i/3 ? 2.0 : 1.0);
      n = 0;
      for (int k = 0; k < nfft; k++) {
/*
        if ( tmp_chi != 0) complex_multiply(vg_hat[loc][i], wki, ktmpi, n);
        if ( kappa != 0 ) complex_multiply(vg_hat[loc][i], wk1, ktmp, n);
*/
        complex_multiply(vg_hat[loc][i], wki, ktmpi, n);
        n += 2;
      }
      fft2->compute(ktmp, ktmp, -1);
      n=0;
      vtmp = 0.0;
      for (int k = 0; k < nfft; k++) {
        // rho_i * u_ij * rho_j * chi * prefactor
/*
        if ( tmp_chi != 0) virial[i] += tmp_chi * ktmpi[n] * density_fft_types[jtype][k] * factor2;
       //fprintf(screen,"factor2 %f %f %f %f %f\n", kappa, vg_hat[loc][i], ktmp[n], wk1[k], factor2);
        // rho_i * u_ij * rho_j * prefactor * kappa
        if ( kappa != 0 ) virial[i] += kappa * ktmp[n] * rho1[k] * factor2;
*/
        vtmp += ktmpi[n] * density_fft_types[jtype][k];
        n+=2;
      }
      virial[i] += vtmp *(tmp_chi + tmp_kappa) * factor2;
    }
   //fprintf(screen,"virial %f %f %f %f %f %f\n", virial[0], virial[1], virial[2], virial[3], virial[4], virial[5]);
  }
  
}
/*
  double tmp_rho_sub = 0.0;
  if (subtract_rho0 == 1) tmp_rho_sub = rho0;
  // Convolve kspace-potential with fft(den_group)
  for ( int itype = 1; itype <= ntypes; itype++) {
    for ( int jtype = itype; jtype <= ntypes; jtype++) {
      if ( potent_type_map[0][itype][jtype] == 1) continue;
      if (fabs(chi[itype][jtype]) == 0  && kappa == 0) continue;

      // take itype density into k-space
      for (int k = 0; k < nfft; k++) {
         work1[n++] = density_fft_types[itype][k] - tmp_rho_sub;
         work1[n++] = ZEROF;
       }
      fft1->compute(work1,work1,1);
  
      n=0;
      for (int k = 0; k < nfft; k++) {
        work1[n++] *= scale_inv;
        work1[n++] *= scale_inv;
      }
      // take jtype density into k-space
  
      // convolve itype-jtype interaction potential and itype density
      int loc = itype*(jtype+1)/2;
      n = 0;
      for (int k = 0; k < nfft; k++) {
        complex_multiply(potent[loc], work1, ktmp2, n);
        n += 2;
      }
      // IFFT the convolution
      fft1->compute(ktmp2, ktmp2, -1);
  
      double same_group_factor = 1.0;
      if (itype == jtype) {
        same_group_factor = 2.0;
        work2 = work1;
      } else {
        // take jtype density into k-space
        for (int k = 0; k < nfft; k++) {
           work2[n++] = density_fft_types[jtype][k] - tmp_rho_sub;
           work2[n++] = ZEROF;
         }
        fft1->compute(work1,work1,1);
      
        n=0;
        for (int k = 0; k < nfft; k++) {
          work2[n++] *= scale_inv;
          work2[n++] *= scale_inv;
        }
      }

      double factor = (chi[itype][jtype] + kappa) / delvolinv * same_group_factor * 0.5 / tmp_rho_div;

      // chi and kappa energy contribution
      if (eflag_global) {
      fprintf(screen,"accu grad E\n");
        n = 0;
        for (int k = 0; k < nfft; k++) {
          // rho_i * u_ij * rho_j * prefactor
          energy += ktmp2[n] * work2[k] * factor;
          n += 2;
        }
      }

      // pressure tensor calculation
      if (vflag_global) {
      fprintf(screen,"accu grad P\n");
        for (int i = 0; i < 6; i++) {
          double factor2 = factor * ( i/3 ? 2.0 : 1.0);
          n = 0;
          for (int k = 0; k < nfft; k++) {
            complex_multiply(vg_hat[loc][i], work1, ktmp, n);
            n += 2;
          }
          fft1->compute(ktmp, ktmp, -1);
          n=0;
          for (int k = 0; k < nfft; k++) {
            virial[i] += ktmp[n] * work2[k] * factor2;
            n+=2;
          }
        }
      }
*/

/* ----------------------------------------------------------------------
  Calculates the rho0 needed for this system instead of the rho0 that 
  comes from the particle densities. 
------------------------------------------------------------------------- */
double TILD::calculate_rho0(){
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;
  int particles_not_tild = 0;
  vector<int> count_per_type (ntypes,0);
  int k;
  double local_rho0 = 0;
  int temp1, temp2;
  double lmass = 0, lmass_all = 0;
  double *mass = atom->mass;

  for (int i = 0; i < nlocal; i++) {
    if ( potent_type_map[0][type[i]][type[i]] == 1) {
      particles_not_tild++;
    } else {
      count_per_type[type[i]]++;
    }
  }

  if (screen) {
    fprintf(screen, "Found %d particles without a TILD potential.\n", particles_not_tild);
  }
  if (logfile) {
    fprintf(logfile, "Found %d particles without a TILD potential.\n", particles_not_tild);
  }

  for ( int itype = 1; itype <= ntypes; itype++) {
    if (mass[itype] == 0) continue;

    // gaussian potential, has no radius
    if (potent_type_map[1][itype][itype] == 1) {
      lmass += count_per_type[itype] * mass[itype];
    } else if (potent_type_map[2][itype][itype] == 1) {
      double volume = 4.0 * PI / 3 * rp[itype][itype] * rp[itype][itype] * rp[itype][itype];
      lmass += count_per_type[itype] * mass[itype] * volume;
    }
  }
  MPI_Allreduce(&lmass, &lmass_all, 1, MPI_DOUBLE, MPI_SUM, world);

  double vole = domain->xprd * domain->yprd * domain->zprd;
  rho0 = lmass_all / vole;
  return rho0;
}

/*
int TILD::identify_potential_for_type(int particle_type){
  return identify_potential_for_type(particle_type, true);
}

int TILD::identify_potential_for_type(int particle_type, bool print_error){

  auto itr = std::find_if( types_and_potentials.cbegin(), 
                           types_and_potentials.cend(), 
                           [&particle_type](std::pair<int, int> const &types_and_potentials) 
                           { return types_and_potentials.first == particle_type; }
                           );

  if (itr != types_and_potentials.cend()) {
	  return distance(types_and_potentials.cbegin(), itr);
	}
	else if (print_error) {
    char str[128];
    sprintf(str,"Unable to find a potential for type %d",particle_type);
    error->all(FLERR,str);
	}
	return -1;

}
*/

/*
void TILD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&chi[i][j],sizeof(double),1,fp);
        fwrite(&a2[i][j],sizeof(double),1,fp);
        fwrite(&rp[i][j],sizeof(double),1,fp);
        fwrite(&xi[i][j],sizeof(double),1,fp);
      }
    }
}

*/
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

/*
void TILD::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&chi[i][j],sizeof(double),1,fp);
          fread(&a2[i][j],sizeof(double),1,fp);
          fread(&rp[i][j],sizeof(double),1,fp);
          fread(&xi[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&chi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rp[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&xi[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}
*/

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

/*
void TILD::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}
*/

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

/*
void TILD::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}
*/


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

/*
void TILD::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,chi[i][i],a2[i][i],rp[i][j],xi[i][j]);
}
*/

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

/*
void TILD::write_data_all(FILE *fp)
{
  int ntypes = atom->ntypes;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,chi[i][j],a2[i][j],rp[i][j],xi[i][j]);
}

*/

