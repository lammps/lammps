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
#include "error.h"
#include <iostream>
#include "fft3d_wrap.h"
#include "pppm.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ------------------------------------------------------------ */

TILD::TILD(LAMMPS *lmp) : PPPM(lmp)  {
    if (screen) fprintf(screen,"TILD construction...\n");
    if (logfile) fprintf(logfile,"TILD construction...\n");
  };

TILD::~TILD(){
  return;
}

void TILD::settings(int narg, char **arg)
{
    if (narg < 1) error->all(FLERR,"Illegal kspace_style tild command");
    accuracy_relative = fabs(force->numeric(FLERR,arg[0]));
//   std::cout<<"help me "<< std::endl;
}

void TILD::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"TILD initialization ...\n");
    if (logfile) fprintf(logfile,"TILD initialization ...\n");
  }
}


void TILD::setup(){
  return;
}

void TILD::compute(int i1, int i2){
  return;
}

void TILD::compute_group_group(int, int, int){
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

  // fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
  //                  nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
  //                  nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
  //                  0,0,&tmp,collective_flag);

  // fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
  //                  nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
  //                  nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
  //                  0,0,&tmp,collective_flag);

  // remap = new Remap(lmp,world,
  //                   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
  //                   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
  //                   1,0,0,FFT_PRECISION,collective_flag);

  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;

  // if (differentiation_flag == 1)
  //   cg = new GridComm(lmp,world,1,1,
  //                     nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
  //                     nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
  //                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
  //                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  // else
  //   cg = new GridComm(lmp,world,3,1,
  //                     nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
  //                     nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
  //                     procneigh[0][0],procneigh[0][1],procneigh[1][0],
  //                     procneigh[1][1],procneigh[2][0],procneigh[2][1]);
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


// Need to create functionality that would loop over all places owned by 
// this processor and use that to generate the position for the gaussian. and assign them as well.
void TILD::init_gauss(){
  
  // Represents the 0 to N-1 points that are captured on this grid.
  int nlo_in, nhi_in;
  int nrange = nhi_in - nlo_in +1;
  int n_loc=0;
  int Dim = domain->dimension;
  double mdr2;
  double V = domain->xprd * domain->yprd * domain->zprd ;


  double pref = V / ( pow(2.0* sqrt(MY_PI * gauss_a2), Dim ));

  // for (double zloc = nzlo_in; zloc < nzhi_in; zloc++) {
  //   for (double yloc = nylo_in; yloc < nyhi_in; yloc++) {
  //     for (double xloc = nxlo_in; xloc < nxhi_in; xloc++) {
  //       mdr2 = xloc ^ 2 + yloc ^ 2 + zloc ^ 2;
  //       uG[n_loc] = exp(-mdr2 / 4.0 / a_squared);
  //     }
  //   }
  // }
}
