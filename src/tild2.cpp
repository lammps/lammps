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
#include "update.h"
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
  vg(NULL), factors(NULL), fkx(NULL), fky(NULL), fkz(NULL), 
  work1(NULL), work2(NULL), rho1d(NULL),  rho_coeff(NULL), drho_coeff(NULL), 
  density_A_brick(NULL), density_B_brick(NULL), density_A_fft(NULL), density_B_fft(NULL), 
  fft1(NULL), fft2(NULL), remap(NULL), cg(NULL), cg_peratom(NULL),  part2grid(NULL), boxlo(NULL)
 {
  if (comm->me == 0) {
    if (screen) fprintf(screen, "TILD construction...\n");
    if (logfile) fprintf(logfile, "TILD construction...\n");
  }
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
  density_brick_types = NULL;
  avg_density_brick_types = NULL;
  density_fft_types = NULL;
  density_hat_fft_types = NULL;
  kappa_density = NULL;
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
  ktmp = NULL;
  ktmp2 = NULL;
  ktmpi = ktmpj = NULL;
  ktmp2i = ktmp2j = NULL;
  tmp = NULL;
  work1 = work2 = NULL;
  vg = NULL;
  vg_hat = NULL;
  fkx = fky = fkz = NULL;
  fkx2 = fky2 = fkz2 = NULL;

  density_A_brick = density_B_brick = NULL;
  density_A_fft = density_B_fft = NULL;

  rho1d = rho_coeff = drho_coeff = NULL;

  fft1 = fft2 = NULL;
  remap = NULL;
  cg = NULL;
  cg_peratom = NULL;
  specified_all_group = 0;

  rho0 = 0.0;
  nmax = 0;
  sub_flag  = 1;
  mix_flag  = 1;
  ave_grid_flag  = 0;
  nvalid_last = -1;
  nvalid = 0;
  nevery = 0;
  irepeat = 0;
  nrepeat = 0;
  peratom_freq = 0;
  write_grid_flag  = 0;
  set_rho0 = 1.0;
  subtract_rho0 = 0;
  norm_flag = 1;
  part2grid = NULL;

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

  if (order > MAXORDER) {
    char str[128];
    sprintf(str,"TILD coulomb order cannot be greater than %d",MAXORDER);
    error->all(FLERR,str);
  }

  // free all arrays previously allocated
  deallocate();
  deallocate_groups();
  deallocate_peratom();

  set_grid();

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
    int i, j;

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

  if (sub_flag == 1) {
    subtract_rho0 = 1;
  } else {
    subtract_rho0 = 0;
  }  
  if (norm_flag == 1) {
    normalize_by_rho0 = 1;
  } else {
    normalize_by_rho0 = 0;
  }

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
  vir_func_init();

  return;
}

void TILD::vir_func_init() {
  int z, y, x, n;
  int Dim = domain->dimension;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  //double zprd_slab = zprd * slab_volfactor;
  double k[Dim];
  double scale_inv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);
  double delx = xprd/nx_pppm;
  double dely = yprd/ny_pppm;
  double delz = zprd/nz_pppm;
  int ntypes = atom->ntypes;

  int loc = 0;
  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      // Skip if type cross-interaction does not use density/tild
      if ( potent_type_map[0][itype][jtype] == 1) continue;

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
            //fprintf(screen,"virial %d %f %f %f %f %f %f\n", n, k[0], k[1], k[2], grad_potent[loc][0][n], grad_potent[loc][1][n], grad_potent[loc][2][n]);
            n++;

            }
          }
      }

      for (int i = 0; i < 6; i++){
        n = 0;
        for (int j = 0; j < nfft; j++) {
          ktmp[n++] = vg[loc][i][j];
          ktmp[n++] = ZEROF;
        }
  
        fft1->compute(ktmp, ktmp2, 1);

        for (int j = 0; j < 2 * nfft; j++) {
          ktmp2[j] *= scale_inv;
          vg_hat[loc][i][j] = ktmp2[j];
        }

      }
/*
      std::string fname = "vg_"+std::to_string(itype)+"-"+std::to_string(jtype)+"_lammps.txt";
      std::ofstream file(fname);
      n = 0;
      for (int j = 0; j <  nfft; j++) {
        file << j <<'\t'<< vg[loc][0][j]  <<'\t'<< vg[loc][1][j]  <<'\t'<< vg[loc][2][j]  <<'\t'<< vg[loc][3][j]  <<'\t'<< vg[loc][4][j] <<'\t'<< vg[loc][5][j] <<'\t'<< vg_hat[loc][0][n]  <<'\t'<< vg_hat[loc][1][n]  <<'\t'<< vg_hat[loc][2][n]  <<'\t'<< vg_hat[loc][3][n]  <<'\t'<< vg_hat[loc][4][n] <<'\t'<< vg_hat[loc][5][n] <<'\t'<< vg_hat[loc][0][n+1]  <<'\t'<< vg_hat[loc][1][n+1]  <<'\t'<< vg_hat[loc][2][n+1]  <<'\t'<< vg_hat[loc][3][n+1]  <<'\t'<< vg_hat[loc][4][n+1] <<'\t'<< vg_hat[loc][5][n+1] << endl;
        n += 2;
      }
*/
      loc++;
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


  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate_peratom(), compute() will allocate when needed

  allocate();

  // pre-compute volume-dependent coeffs
  compute_rho_coeff(rho_coeff, drho_coeff, order);

  rho0 = calculate_rho0();

  cg->ghost_notify();
  if (overlap_allowed == 0 && cg->ghost_overlap())
    error->all(FLERR,"PPPM grid stencil extends "
               "beyond nearest neighbor processor");
  cg->setup();

  setup();
}

void TILD::precompute_density_hat_fft() {
  int n = 0;
  double scale_inv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);
  int ntypes = atom->ntypes;

  for ( int ktype = 0; ktype <= ntypes; ktype++) {
    n = 0;
    for (int k = 0; k < nfft; k++) {
      //fprintf(screen,"rho %d %d %f\n", ktype, k,density_fft_types[ktype][k]);
      work1[n++] = density_fft_types[ktype][k];
      work1[n++] = ZEROF;
    }

    // FFT the density to k-space
    fft1->compute(work1, work1, 1);

    for (int k = 0; k < 2*nfft; k++) {
      work1[k] *= scale_inv;
      density_hat_fft_types[ktype][k] = work1[k];
    }

/*
    std::string fname = "rho_"+std::to_string(ktype)+"_lammps.txt";
    std::ofstream rhof(fname);
    n = 0;
    for (int k = 0; k < nfft; k++){
      rhof << k << '\t' << density_fft_types[ktype][k] << '\t' << density_hat_fft_types[ktype][n] << '\t' << density_hat_fft_types[ktype][n+1] << endl;
      n += 2;
    }
*/
  }
}

void TILD::compute(int eflag, int vflag){

  if (domain->box_change) {
    rho0 = calculate_rho0();
  }
  
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

  make_rho();

  cg->reverse_comm(this, REVERSE_RHO_NONE);

  brick2fft();

  accumulate_gradient();

  cg->forward_comm(this, FORWARD_NONE);

  fieldforce_param();

  if ( write_grid_flag == 1 ) {
    if ( (update->ntimestep % grid_data_output_freq) ==  0) {
      write_grid_data(grid_data_filename, 0);
    }
  }
  if ( ave_grid_flag == 1 ) ave_grid();

/*
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
*/
  if (atom->natoms != natoms_original) {
    natoms_original = atom->natoms;
  }
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
  int ntypecross = ((ntypes-1)*ntypes) - int(0.5*(ntypes-2)*(ntypes-1));

  // style coeffs
  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;

  //memory->create(potential_type_list,ntypes,ntypes+1,"pppm:potential_type_list"); // lets you know if a type has a potential

  memory->create4d_offset(density_brick_types,ntypes+1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:density_brick_types");
  memory->create4d_offset(avg_density_brick_types,ntypes+1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:avg_density_brick_types");
  memory->create5d_offset(gradWtype,ntypes+1, 0, Dim-1,
                          nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"tild:gradWtype");

  memory->create(work1,2*nfft_both,"tild:work1");
  memory->create(work2,2*nfft_both,"tild:work2");
  memory->create(ktmp,2*nfft_both,"tild:ktmp");
  memory->create(ktmpi,2*nfft_both,"tild:ktmpi");
  memory->create(ktmpj,2*nfft_both,"tild:ktmpj");
  memory->create(ktmp2,2*nfft_both, "tild:ktmp2");
  memory->create(ktmp2i,2*nfft_both, "tild:ktmp2i");
  memory->create(ktmp2j,2*nfft_both, "tild:ktmp2j");
  memory->create(tmp, nfft, "tild:tmp");
  memory->create(vg,ntypecross+1,6,nfft_both,"tild:vg");
  memory->create(vg_hat,ntypecross+1,6,2*nfft_both,"tild:vg_hat");
  memory->create(density_fft_types,ntypes+1, nfft_both, "tild:density_fft_types");
  memory->create(density_hat_fft_types,ntypes+1, 2*nfft_both, "tild:density_hat_fft_types");
  memory->create(kappa_density, nfft_both, "tild:kappa_density");
  memory->create(potent,ntypecross+1,nfft_both,"tild:potent"); // Voignot 
  memory->create(potent_hat,ntypecross+1,2*nfft_both,"tild:potent_hat");
  memory->create(grad_potent,ntypecross+1,domain->dimension,nfft_both,"tild:grad_potent");
  memory->create(grad_potent_hat,ntypecross+1, domain->dimension,2*nfft_both,"tild:grad_potent_hat");
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

  // summation coeffs

  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,"pppm:drho_coeff");

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
  memory->destroy(setflag);

  //memory->destroy(potent_type_map);

  memory->destroy4d_offset(density_brick_types,nzlo_out,nylo_out,nxlo_out);
  memory->destroy4d_offset(avg_density_brick_types,nzlo_out,nylo_out,nxlo_out);
  memory->destroy(kappa_density);
  memory->destroy(density_fft_types);
  memory->destroy(density_hat_fft_types);
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

  int order_allocated = order;
  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);

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
   set size of FFT grid (nx,ny,nz_pppm) and g_ewald
   for Coulomb interactions
------------------------------------------------------------------------- */

void TILD::set_grid()
{
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
  //bigint natoms = atom->natoms;

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

  // int Dim = domain->dimension;
  // memory->create5d_offset(gradWtype,group->ngroup, 0, Dim,
  //                         nzlo_out,nzhi_out,nylo_out,nyhi_out,
  //                         nxlo_out,nxhi_out,"tild:gradWtype");
  // memory->create(gradWtype, group->ngroup, Dim,, "tild:gradWtype");

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

  delete cg_peratom;
  cg_peratom = NULL;
}

void TILD::init_cross_potentials(){
  

  // decomposition of FFT mesh
  // global indices range from 0 to N-1
  // proc owns entire x-dimension, clumps of columns in y,z dimensions
  // npey_fft,npez_fft = # of procs in y,z dims
  // if nprocs is small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own

  int Dim = domain->dimension;
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

  int ntypes = atom->ntypes;
  double scale_inv = 1.0/ nx_pppm/ ny_pppm/ nz_pppm;
  int n = 0;

  // Loop over potental styles
  int loc = 0;
  for (int itype = 1; itype <= ntypes; itype++) {
    // Skip if type cross-interaction does not use density/tild

    for (int jtype = itype; jtype <= ntypes; jtype++) {
      // Skip if type cross-interaction does not use density/tild
      if ( potent_type_map[0][itype][jtype] == 1) continue;


      // If both parameters are Gaussian, just do analytical convolution
      if (potent_type_map[1][itype][jtype] == 1 or (mix_flag == 1 && potent_type_map[1][itype][itype] == 1 && potent_type_map[1][jtype][jtype] == 1) ){
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
  
        fft1->compute(ktmp, ktmp2, 1);

        for (int i = 0; i < 2 * nfft; i++) {
          ktmp2[i] *= scale_inv;
          potent_hat[loc][i] = ktmp2[i];
        }
        //init_potential_ft(potent_hat[loc], 1, p);
      } 
      // Computational Convolution
      else {

        // calculate 1st and 2nd potentials to convolv
        //fprintf(screen,"i j %d %d %d\n", itype, jtype, mix_flag);
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
          complex_multiply(work1, work2, ktmp2, n);
          potent_hat[loc][n] = ktmp2[n];
          potent_hat[loc][n+1] = ktmp2[n+1];
          n += 2;
        }

        fft1->compute(ktmp2, ktmp, -1);

        n = 0;
        for (int j = 0; j < nfft; j++){
          potent[loc][j] = ktmp[n];
          n += 2;
        }
      }

      get_k_alias(potent_hat[loc], grad_potent_hat[loc]);
      for (int i=0; i < Dim; i ++){
        for (int j = 0; j < 2*nfft; j++) {
          work1[j] = grad_potent_hat[loc][i][j];
        }
        fft2->compute(work1, work2, -1);
        n = 0;
        for (int j = 0; j < nfft; j++) {
          grad_potent[loc][i][j] = -work2[n]; // still not sure about minus sign ..
          n+=2;
        }
      } 

      // output
/*
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
          fileU << j << '\t' << potent[loc][j] << '\t' << potent_hat[loc][n] << '\t' << potent_hat[loc][n+1] << std::endl;
          filegradU << j << '\t' << grad_potent[loc][0][j] << '\t' << grad_potent[loc][1][j] << '\t' << grad_potent[loc][2][j] << std::endl;
          filegradUhat << j << '\t' << grad_potent_hat[loc][0][n] << '\t' << grad_potent_hat[loc][1][n] << '\t' << grad_potent_hat[loc][2][n] << std::endl;
          filegradUhatI << j << '\t' << grad_potent_hat[loc][0][n+1] << '\t' << grad_potent_hat[loc][1][n+1] << '\t' << grad_potent_hat[loc][2][n+1] << std::endl;
          n += 2;
      }

      fileU.close();
      filegradU.close();
      filegradUhat.close();
      filegradUhatI.close();
*/
      loc++;
    }
  }

}

int TILD::get_style( const int i, const int j) {
  for (int istyle = 1; istyle <= nstyles; istyle++) { 
    if ( potent_type_map[istyle][i][j] == 1 ) return istyle;
  }
  return 0;
}

void TILD::calc_work(FFT_SCALAR *wk, const int itype, const int jtype){
  double scale_inv = 1.0/ nx_pppm/ ny_pppm/ nz_pppm;

  // needs work is this right for cross terms of the same potential?
  
  double params[4];
  for (int i = 0; i < 4; i++) params[i] = 0.0;

  int style = get_style(itype, jtype);
  if (style == 1) {
    params[0] = a2[itype][jtype];
    init_potential_ft(wk, style, params);
  } else if (style == 2) {
    params[0] = rp[itype][jtype];
    params[1] = xi[itype][jtype];
    init_potential(tmp, style, params);

    int j = 0;
    for (int i = 0; i < nfft; i++) {
      ktmp[j++] = tmp[i];
      ktmp[j++] = ZEROF;
    }

    fft1->compute(ktmp, wk, 1);

    for (int i = 0; i < 2 * nfft; i++) {
      wk[i] *= scale_inv;
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

  if (type == 1) {
                                                                // should be 3/2 right?
    for (int z = nzlo_fft; z <= nzhi_fft; z++) {
      zper = static_cast<double>(z) / zprd;
      if (z >= nz_pppm / 2.0) {
        zper -= static_cast<double>(nz_pppm) / zprd;
      }
      double zper2 = factor * zper * zper; 

      for (int y = nylo_fft; y <= nyhi_fft; y++) {
        yper = static_cast<double>(y) / yprd;
        if (y >= ny_pppm / 2.0) {
          yper -= static_cast<double>(ny_pppm) / yprd;
        }
        double yper2 = factor * yper * yper; 

        for (int x = nxlo_fft; x <= nxhi_fft; x++) {
          xper = static_cast<double>(x) / xprd;
          if (x >= nx_pppm / 2.0) {
            xper -= static_cast<double>(nx_pppm) / xprd;
          }

          k2 = (factor * xper * xper) + yper2 + zper2;
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
    double pref = vole / (pow( sqrt(2.0 * MY_PI * (parameters[0]) ), Dim));
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

  int has_nyquist = 0;
  if ( nz_pppm % 2 == 0 || ny_pppm % 2 == 0 || nx_pppm % 2 == 0 ) {
     has_nyquist = 1;
  }

  // id2 = unstack_stack( id ) ;
  // unstack(id2, n);
  for (z = nzlo_fft; z <= nzhi_fft; z++) {
    if (nz_pppm % 2 == 0 && z == nz_pppm / 2)
      k[2] = 0.0;
    else if (double(z) < double(nz_pppm) / 2.0)
      k[2] = 2 * MY_PI * double(z) / zprd;
    else
      k[2] = 2 * MY_PI * double(z - nz_pppm) / zprd;

    for (y = nylo_fft; y <= nyhi_fft; y++) {
      // skip to kill off Nyquist modes
      
      //if (k[2] == 0.0 || (ny_pppm % 2 == 0 && y == ny_pppm / 2)) {
      //  k[2] = 0.0;
      if (ny_pppm % 2 == 0 && y == ny_pppm / 2) {
        k[1] = 0.0;
      } else if (double(y) < double(ny_pppm) / 2.0)
        k[1] = 2 * MY_PI * double(y) / yprd;
      else
        k[1] = 2 * MY_PI * double(y - ny_pppm) / yprd;

      for (x = nxlo_fft; x <= nxhi_fft; x++) {
        //if (k[2] == 0.0 || k[1] == 0.0 || (nx_pppm % 2 == 0 && x == nx_pppm / 2)) {
        //  k[2] = 0.0;
        //  k[1] = 0.0;
        if (nx_pppm % 2 == 0 && x == nx_pppm / 2) {
          k[0] = 0.0;
        } else if (double(x) < double(nx_pppm) / 2.0)
          k[0] = 2 * MY_PI * double(x) / xprd;
        else
          k[0] = 2 * MY_PI * double(x - nx_pppm) / xprd;
/*
        if (has_nyquist && n>0) {
          if (k[2] == 0.0 || k[1] == 0.0 || k[0] == 0.0) {
            k[2] = 0.0;
            k[1] = 0.0;
            k[0] = 0.0;
          }
        }
*/
        out[0][n] = -wk1[n + 1] * k[0];
        out[0][n + 1] = wk1[n] * k[0];
        out[1][n] = -wk1[n + 1] * k[1];
        out[1][n + 1] = wk1[n] * k[1];
        out[2][n] = -wk1[n + 1] * k[2];
        out[2][n + 1] = wk1[n] * k[2];
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
      
      }
  }

 if (flag) error->one(FLERR,"Out of range atoms - cannot compute TILD");
}


int TILD::modify_param(int narg, char** arg)
{
  int ntypes = atom->ntypes;
  if (strcmp(arg[0], "tild/chi") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, "TILD command before simulation box is defined");

      if (narg != 4) error->all(FLERR, "Illegal kspace_modify tild command");

      int ilo,ihi,jlo,jhi;
      force->bounds(FLERR,arg[1],ntypes,ilo,ihi);
      force->bounds(FLERR,arg[2],ntypes,jlo,jhi);
      double chi_one = force->numeric(FLERR,arg[3]);
      for (int i = ilo; i <= ihi; i++) {
        for (int j = MAX(jlo,i); j <= jhi; j++) {
          chi[i][j] = chi_one;
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

  } else if (strcmp(arg[0], "tild/mix") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      mix_flag = 1;
      if (strcmp(arg[1], "convolution") == 0)  mix_flag = 1;
      else if (strcmp(arg[1], "define") == 0) mix_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify tild mix argument");
 } else if (strcmp(arg[0], "tild/set_rho0") == 0) {
     if (narg < 2 ) error->all(FLERR, "Illegal kspace_modify tild command");
     set_rho0 = force->numeric(FLERR,arg[1]);
  } else if (strcmp(arg[0], "tild/subtract_rho0") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      if (strcmp(arg[1], "yes") == 0) sub_flag = 1;
      else if (strcmp(arg[1], "no") == 0) sub_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify tild subtract_rho0 argument");
  } else if (strcmp(arg[0], "tild/normalize_by_rho0") == 0) {
      if (narg != 2) error->all(FLERR, "Illegal kspace_modify tild command");
      if (strcmp(arg[1], "yes") == 0) norm_flag = 1;
      else if (strcmp(arg[1], "no") == 0) norm_flag = 0;
      else 
        error->all(FLERR, "Illegal kspace_modify tild normalize_by_rho0 argument");

  } else if (strcmp(arg[0], "tild/write_grid_data") == 0) {
      if (narg != 3) error->all(FLERR, "Illegal kspace_modify tild command");
      write_grid_flag = 1;
      grid_data_output_freq = force->inumeric(FLERR,arg[1]);
      strcpy(grid_data_filename,arg[2]);

  } else if (strcmp(arg[0], "tild/ave/grid") == 0) {
      if (narg != 5) error->all(FLERR, "Illegal kspace_modify tild command");
      ave_grid_flag = 1;
      nevery = force->inumeric(FLERR,arg[1]);
      nrepeat = force->inumeric(FLERR,arg[2]);
      peratom_freq = force->inumeric(FLERR,arg[3]);
      strcpy(ave_grid_filename,arg[4]);
      nvalid = nextvalid();
      if (nevery <= 0 || nrepeat <= 0 || peratom_freq <= 0)
        error->all(FLERR,"Illegal fix tild/ave/grid command");
      if (peratom_freq % nevery || nrepeat*nevery > peratom_freq)
        error->all(FLERR,"Illegal kspace_modify tild/ave/grid command");
  } else
    error->all(FLERR, "Illegal kspace_modify tild command");

  return narg;
}

void TILD::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  int Dim = domain->dimension;

  if (flag == FORWARD_NONE){
    for (int ktype = 0; ktype <= atom->ntypes; ktype++) {
      for (int j = 0; j < Dim; j++) {
        FFT_SCALAR *src = &gradWtype[ktype][j][nzlo_out][nylo_out][nxlo_out];
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
    for (int ktype = 0; ktype <= atom->ntypes; ktype++) {
      for (int j = 0; j < Dim; j++) {
        FFT_SCALAR *dest = &gradWtype[ktype][j][nzlo_out][nylo_out][nxlo_out];
        for (int i = 0; i < nlist; i++)
          dest[list[i]] = buf[n++];
      }
    }
  }
}

void TILD::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list) {
  int n = 0;
  if (flag == REVERSE_RHO_NONE) {
    for (int ktype = 0; ktype <= atom->ntypes; ktype++) {
      FFT_SCALAR *src = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
      for (int i = 0; i < nlist; i++)
        buf[n++] = src[list[i]];
    }
  }
}

void TILD::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  if (flag == REVERSE_RHO_NONE) {
    for (int ktype = 0; ktype <= atom->ntypes; ktype++) {
      FFT_SCALAR *dest = &density_brick_types[ktype][nzlo_out][nylo_out][nxlo_out];
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

void TILD::make_rho()
{
  int nx,ny,nz,mx,my,mz;
  int ntypes = atom->ntypes;
  FFT_SCALAR dx,dy,dz,x0,y0,z0,w;

  for (int k = 0; k <= ntypes; k++) {
    memset(&(density_brick_types[k][nzlo_out][nylo_out][nxlo_out]),0,
           ngrid*sizeof(FFT_SCALAR));
  }


  // loop over my particles, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {
    
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
    // loop over the stencil
    for (int n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (int m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (int l = nlower; l <= nupper; l++) {
          mx = l+nx;
          w = x0*rho1d[0][l];

          density_brick_types[type[i]][mz][my][mx] += w;
          density_brick_types[0][mz][my][mx] += w;
        }
      }
    }
  }
}

void TILD::brick2fft()
{
  int ix,iy,iz;
  int ntypes = atom->ntypes;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  for (int k = 0; k <= ntypes; k++) {
    int n = 0;
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++){
          density_fft_types[k][n++] = density_brick_types[k][iz][iy][ix];
          //rhof<<n-1<<'\t'<<density_fft_types[k][n-1] <<std::endl;
        }
  }

  for (int k = 0; k <= ntypes; k++) {
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
  int n = 0;
  int ntypes = atom->ntypes;
  
  precompute_density_hat_fft(); // calculate rho_hat for each type

  for ( int ktype = 1; ktype <= ntypes; ktype++)
    for (int i = 0; i < Dim; i++)
      memset(&(gradWtype[ktype][i][nzlo_out][nylo_out][nxlo_out]), 0,
             ngrid * sizeof(FFT_SCALAR));

  //typedef std::tuple<int, int, double, int, int> tup_iidii;

  // This part is for the specified chi interactions.
  // Kappa (incompressibility) is later

  double tmp_kappa = kappa;
  double tmp2_kappa;

  FFT_SCALAR tmp_sub = 0.0;
  if (subtract_rho0 == 1)
    tmp_sub = rho0;

  if (normalize_by_rho0 == 1) {
    tmp_kappa = kappa/rho0;
  }

  int loc = 0;
  for (int itype = 1; itype <= ntypes; itype++) {
    for (int jtype = itype; jtype <= ntypes; jtype++) {
      if ( potent_type_map[0][itype][jtype] == 1) continue;
      //fprintf(screen,"i j %d %d\n", itype,jtype);

      double tmp_chi = chi[itype][jtype];
      if (normalize_by_rho0 == 1) tmp_chi /= rho0;

      if ( tmp_chi == 0 && tmp_kappa == 0 ) continue;

      int calc_kappa_rho0 = 0;
      int diff_type = 1;
      if (itype == jtype) {
        diff_type = 0;
        tmp2_kappa = tmp_kappa/2.0;
      } else {
        tmp2_kappa = tmp_kappa;
        if ( subtract_rho0 == 1 ) {
          calc_kappa_rho0 = 1;
        }
      }

      if (eflag_global || vflag_global) {
        ev_calculation(loc, itype, jtype);
      }

      for (int i = 0; i < Dim; i++) {

        n = 0;
        for (int k = 0; k < nfft; k++) {
          complex_multiply(grad_potent_hat[loc][i], density_hat_fft_types[itype], ktmp2i, n);
          //fprintf(screen,"%d %d %d %d %f %f %f %f %f %f\n", itype, jtype, k, i, grad_potent_hat[loc][i][n], grad_potent_hat[loc][i][n+1], density_hat_fft_types[itype][n], density_hat_fft_types[itype][n+1], ktmp2i[n], ktmp2i[n+1]);
          if (diff_type) complex_multiply(grad_potent_hat[loc][i], density_hat_fft_types[jtype], ktmp2j, n);
          n += 2;
        }
        fft1->compute(ktmp2i, ktmpi, -1);
        if (diff_type) fft1->compute(ktmp2j, ktmpj, -1);

        n = 0;
        int j = 0;
        for (int k = nzlo_in; k <= nzhi_in; k++)  {
          for (int m = nylo_in; m <= nyhi_in; m++) {
            for (int o = nxlo_in; o <= nxhi_in; o++) {
              gradWtype[jtype][i][k][m][o] += ktmpi[n] * (tmp_chi + tmp_kappa);
              if (diff_type) gradWtype[itype][i][k][m][o] += ktmpj[n] * (tmp_chi + tmp_kappa);
              if (calc_kappa_rho0 == 1) {
                gradWtype[itype][i][k][m][o] -= grad_potent[loc][i][j] * rho0 * tmp_kappa;
                if (diff_type) gradWtype[jtype][i][k][m][o] -= grad_potent[loc][i][j] * rho0 * tmp_kappa;
              }
              n += 2;
              j++;
            }
          }
        }
      }
      loc++;
    }
  }
/*
  for (int itype = 1; itype <= ntypes; itype++) {
    std::string fname = "gradw_"+std::to_string(itype)+"_lammps.txt";
    std::ofstream gradtype(fname);
    n = 0;
    for (int k = nzlo_in; k <= nzhi_in; k++)
      for (int m = nylo_in; m <= nyhi_in; m++)
        for (int o = nxlo_in; o <= nxhi_in; o++)
          gradtype << n++ <<'\t'<< gradWtype[itype][0][k][m][o] <<'\t'<< gradWtype[itype][1][k][m][o] <<'\t'<<  gradWtype[itype][2][k][m][o]<< endl;
  }
*/
}

void TILD::fieldforce_param(){
  int nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  // Convert field to force per particle
  //ofstream delvolinv_forces("delv_forces_lammps.txt");
  // ofstream grid_vol_forces("gridvol_forces_lammps.txt");
  // ofstream both_forces("delv_gridvol_forces_lammps.txt");
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {

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
    for (int n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (int m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (int l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          // gradWpotential 0 to 1 ?????
          ekx += x0 * gradWtype[temp_type][0][mz][my][mx];
          eky += x0 * gradWtype[temp_type][1][mz][my][mx];
          ekz += x0 * gradWtype[temp_type][2][mz][my][mx];
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
    //delvolinv_forces << i <<'\t'<<  f[i][0] <<'\t'<< f[i][1] <<'\t'<<  f[i][2] << endl;
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
//delvolinv_forces.close();
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
                                 const int ord)
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a; // = NULL;
  memory->create2d_offset(a,ord,-ord,ord,"pppm:a");

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

void TILD::complex_multiply(FFT_SCALAR *in1, FFT_SCALAR  *in2, FFT_SCALAR  *out, int n){
  out[n] = ((in1[n] * in2[n]) - (in1[n + 1] * in2[n + 1]));
  out[n + 1] = ((in1[n + 1] * in2[n]) + (in1[n] * in2[n + 1]));
}

void TILD::ev_calculation(const int loc, const int itype, const int jtype) {
  int n = 0;
  double eng, engk;
  double vtmp, vtmpk;
  double scale_inv = 1.0/(nx_pppm *ny_pppm * nz_pppm);
  double V = domain->xprd * domain->yprd * domain->zprd;
  //int ntypes = atom->ntypes;
  //int loc = (jtype - itype) + ((itype-1)*ntypes) - (0.5*(itype-2)*(itype-1));

  double tmp_rho_div = 1.0;
  if (normalize_by_rho0 == 1) tmp_rho_div = rho0;

  double tmp_kappa = kappa;
  int calc_kappa_rho0 = 0;
  if (itype == jtype) {
    tmp_kappa /= 2.0;
  } else {
    if ( subtract_rho0 == 1 ) {
      calc_kappa_rho0 = 1;
    }
  }
 
  double factor = scale_inv / tmp_rho_div;
  if (eflag_global) {
    // convolve itype-jtype interaction potential and itype density
    n = 0;
    for (int k = 0; k < nfft; k++) {
      complex_multiply(density_hat_fft_types[itype], potent_hat[loc], ktmpi, n);
      n += 2;
    }
    // IFFT the convolution
    fft1->compute(ktmpi, ktmp2i, -1);

    // chi and kappa energy contribution
    n = 0;
    eng = 0.0;
    engk = 0.0;
    for (int k = 0; k < nfft; k++) {
      // rho_i * u_ij * rho_j * chi * prefactor
      eng += ktmp2i[n] * density_fft_types[jtype][k];
      if ( calc_kappa_rho0 == 1 ) engk -= ktmp2i[n];
      n += 2;
    }
    energy += eng * factor * (chi[itype][jtype] + tmp_kappa) * V; 
    //fprintf(screen,"ene %d %d %f %f\n", itype, jtype, eng * factor * (tmp_kappa) * V, eng * factor * (chi[itype][jtype]) * V);
    if ( calc_kappa_rho0 == 1 ) { 
      energy += engk * rho0 * factor * kappa;
      energy += kappa * rho0 * rho0 * factor * nfft / 2.0;
    }
  }

  // pressure tensor calculation
  if (vflag_global) {
    // loop over stress tensor
    for (int i = 0; i < 6; i++) {
      double factor2 = factor;//  * ( i/3 ? 2.0 : 1.0);
      n = 0;
      for (int k = 0; k < nfft; k++) {
        complex_multiply(density_hat_fft_types[itype], vg_hat[loc][i], ktmp2i, n);
        n += 2;
      }
      fft1->compute(ktmp2i, ktmpi, -1);
      vtmp = 0.0;
      vtmpk = 0.0;
      n=0;
      for (int k = 0; k < nfft; k++) {
        vtmp += ktmpi[n] * density_fft_types[jtype][k];
        if ( calc_kappa_rho0 == 1 ) vtmpk -= ktmpi[n];
        n+=2;
      }
      //fprintf(screen,"virial %d %d %d %f\n", i, itype, jtype, vtmp * (chi[itype][jtype] + tmp_kappa) * factor2);
      virial[i] += vtmp * (chi[itype][jtype] + tmp_kappa) * factor2;
      if ( i < 3 ) { // if diagonal member of the matrix, not sure why
        if ( calc_kappa_rho0 == 1 ) {
          virial[i] += vtmpk * rho0 * tmp_kappa * factor2;
          virial[i] += rho0 * kappa;
        }
      }
    }
    //fprintf(screen,"virial %d %d %f %f %f %f %f %f\n", itype, jtype, virial[0], virial[1], virial[2], virial[3], virial[4], virial[5]);
  }
  
}

/* ----------------------------------------------------------------------
  Calculates the rho0 needed for this system instead of the rho0 that 
  comes from the particle densities. 
------------------------------------------------------------------------- */
double TILD::calculate_rho0(){
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;
  int particles_not_tild = 0;
  vector<int> count_per_type (ntypes+1,0);
  double lmass = 0, lmass_all = 0;

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

    // gaussian potential, has no radius
    if (potent_type_map[1][itype][itype] == 1) {
      lmass += count_per_type[itype];
    } else if (potent_type_map[2][itype][itype] == 1) {
      double volume = (4.0 * MY_PI / 3.0) * rp[itype][itype] * rp[itype][itype] * rp[itype][itype] * set_rho0;
      lmass += count_per_type[itype] * volume;
    }
  }
  MPI_Allreduce(&lmass, &lmass_all, 1, MPI_DOUBLE, MPI_SUM, world);

  double vole = domain->xprd * domain->yprd * domain->zprd;

  rho0 = lmass_all / vole;
  if (logfile) {
    fprintf(logfile, "User set rho0 = %f; actual rho0 = %f for TILD potential.\n", set_rho0, rho0);
  }
  return rho0;
}

void TILD::write_grid_data( char *filename, const int avg ) {
  int ntypes = atom->ntypes;
  if (comm->me == 0) {
    otp = fopen( filename, "w" ) ;

    // header
    fprintf( otp, "# x y z") ;
    for (int itype = 1; itype <= ntypes; itype++)
      fprintf( otp, " rho_%d", itype) ;
    fprintf( otp, "\n") ;
  }
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = ntypes + 3;

  int sendrow = nfft_both;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"TILD:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"TILD:buf");

  // pack my atom data into buf

  if (avg == 1) pack_avg_grid_data(buf);
  else pack_grid_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

     for (int n = 0; n < recvrow; n++) {
         fprintf( otp, "%lf %lf %lf", buf[n][0], buf[n][1], buf[n][2]);
         for (int itype = 1; itype <= ntypes; itype++) {
           fprintf( otp, " %1.16e", buf[n][2+itype]);
         }
         fprintf( otp, "\n" ) ;
       }
    } 
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
  if (me == 0) fclose( otp ) ;
}

void TILD::pack_grid_data(double **buf)
{
  int ntypes = atom->ntypes;
  double fx = domain->xprd/nx_pppm;
  double fy = domain->yprd/ny_pppm;
  double fz = domain->zprd/nz_pppm;
  int n = nzlo_in + nylo_in + nxlo_in;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n][0] = ix * fx;
        buf[n][1] = iy * fy;
        buf[n][2] = iz * fz;
        for (int itype = 1; itype <= ntypes; itype++) {
          buf[n][2+itype] = density_brick_types[itype][iz][iy][ix];
        }
        n++;
      }
    }
  }  
}

void TILD::pack_avg_grid_data(double **buf)
{
  int ntypes = atom->ntypes;
  double fx = domain->xprd/nx_pppm;
  double fy = domain->yprd/ny_pppm;
  double fz = domain->zprd/nz_pppm;
  int n = nzlo_in + nylo_in + nxlo_in;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        buf[n][0] = ix * fx;
        buf[n][1] = iy * fy;
        buf[n][2] = iz * fz;
        for (int itype = 1; itype <= ntypes; itype++) {
          buf[n][2+itype] = avg_density_brick_types[itype][iz][iy][ix];
        }
        n++;
      }
    }
  }  
}

void TILD::sum_grid_data()
{
  int ntypes = atom->ntypes;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        for (int itype = 1; itype <= ntypes; itype++) {
          avg_density_brick_types[itype][iz][iy][ix] += density_brick_types[itype][iz][iy][ix];
        }
      }
    }
  }  
}

void TILD::multiply_ave_grid_data(const double factor)
{
  int ntypes = atom->ntypes;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++) {
    for (int iy = nylo_in; iy <= nyhi_in; iy++) {
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        for (int itype = 1; itype <= ntypes; itype++) {
          avg_density_brick_types[itype][iz][iy][ix] *= factor;
        }
      }
    }
  }  
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint TILD::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}

void TILD::ave_grid()
{
  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/atom");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  if (irepeat == 0)
    multiply_ave_grid_data( 0.0 ); 

  // accumulate results of attributes,computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  sum_grid_data(); 

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;

  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+peratom_freq - (nrepeat-1)*nevery;

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  
  multiply_ave_grid_data( 1.0 / repeat );
  write_grid_data(ave_grid_filename, 1);
}

