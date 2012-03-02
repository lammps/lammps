/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_gpu.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER 7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

// External functions from cuda library for atom decomposition

#ifdef FFT_SINGLE
#define PPPM_GPU_API(api)  pppm_gpu_ ## api ## _f
#else
#define PPPM_GPU_API(api)  pppm_gpu_ ## api ## _d
#endif

FFT_SCALAR* PPPM_GPU_API(init)(const int nlocal, const int nall, FILE *screen,
			       const int order, const int nxlo_out, 
			       const int nylo_out, const int nzlo_out,
			       const int nxhi_out, const int nyhi_out,
			       const int nzhi_out, FFT_SCALAR **rho_coeff,
			       FFT_SCALAR **_vd_brick, 
			       const double slab_volfactor,
			       const int nx_pppm, const int ny_pppm,
			       const int nz_pppm, const bool split, 
			       int &success);
void PPPM_GPU_API(clear)(const double poisson_time);
int PPPM_GPU_API(spread)(const int ago, const int nlocal, const int nall,
                      double **host_x, int *host_type, bool &success,
                      double *host_q, double *boxlo, const double delxinv,
                      const double delyinv, const double delzinv);
void PPPM_GPU_API(interp)(const FFT_SCALAR qqrd2e_scale);
double PPPM_GPU_API(bytes)();
void PPPM_GPU_API(forces)(double **f);

/* ---------------------------------------------------------------------- */

PPPMGPU::PPPMGPU(LAMMPS *lmp, int narg, char **arg) : PPPM(lmp, narg, arg)
{
  if (narg != 1) error->all(FLERR,"Illegal kspace_style pppm/gpu command");

  density_brick_gpu = vd_brick = NULL;
  kspace_split = false;
  im_real_space = false;

  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error); 
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPMGPU::~PPPMGPU()
{
  PPPM_GPU_API(clear)(poisson_time);
}

/* ----------------------------------------------------------------------
   called once before run 
------------------------------------------------------------------------- */

void PPPMGPU::init()
{
  PPPM::init();
  
  if (strcmp(update->integrate_style,"verlet/split") == 0)
    kspace_split=true;

  if (kspace_split && universe->iworld == 0) {
    im_real_space = true;
    return;
  }

  // GPU precision specific init

  if (order>8)
    error->all(FLERR,"Cannot use order greater than 8 with pppm/gpu.");
  PPPM_GPU_API(clear)(poisson_time);

  int success;
  FFT_SCALAR *data, *h_brick;
  h_brick = PPPM_GPU_API(init)(atom->nlocal, atom->nlocal+atom->nghost, screen,
			       order, nxlo_out, nylo_out, nzlo_out, nxhi_out,
			       nyhi_out, nzhi_out, rho_coeff, &data, 
			       slab_volfactor,nx_pppm,ny_pppm,nz_pppm,
			       kspace_split,success);

  GPU_EXTRA::check_flag(success,error,world);

  density_brick_gpu =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
		     nxlo_out,nxhi_out,"pppm:density_brick_gpu",h_brick,1);
  vd_brick =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
		     nxlo_out,nxhi_out,"pppm:vd_brick",data,4);

  poisson_time=0;
}

/* ----------------------------------------------------------------------
   compute the PPPMGPU long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPMGPU::compute(int eflag, int vflag)
{
  if (im_real_space) return;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = eflag_global = vflag_global = eflag_atom = vflag_atom = 0;

  if (!peratom_allocate_flag && (eflag_atom || vflag_atom)) {
    allocate_peratom();
    peratom_allocate_flag = 1;
  }

  bool success = true;
  int flag=PPPM_GPU_API(spread)(neighbor->ago, atom->nlocal, atom->nlocal + 
			     atom->nghost, atom->x, atom->type, success,
			     atom->q, domain->boxlo, delxinv, delyinv,
			     delzinv);
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");
  if (flag != 0)
    error->one(FLERR,"Out of range atoms - cannot compute PPPM");

  // If need per-atom energies/virials, also do particle map on host
  // concurrently with GPU calculations

  if (evflag_atom) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
    particle_map();
  }

  int i,j;

  // convert atoms from box to lamda coords
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  double t3=MPI_Wtime();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  
  poisson();

  // all procs communicate E-field values to fill ghost cells
  //   surrounding their 3d bricks

  fillbrick();

  poisson_time+=MPI_Wtime()-t3;

  // calculate the force on my particles

  FFT_SCALAR qscale = force->qqrd2e * scale;
  PPPM_GPU_API(interp)(qscale);

  // Compute per-atom energy/virial on host if requested

  if (evflag_atom) {
    fillbrick_peratom();
    fieldforce_peratom();
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
	eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum / 
	  (g_ewald*g_ewald*volume);
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*q[i]*qscale;
    }
  }

  // sum energy across procs and add in volume-dependent term

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
   
    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/1.772453851 +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // 2d slab correction

  if (slabflag) slabcorr();

  // convert atoms back from lamda to box coords
  
  if (triclinic) domain->lamda2x(atom->nlocal);

  if (kspace_split) PPPM_GPU_API(forces)(atom->f);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMGPU::allocate()
{
  memory->create(density_fft,nfft_both,"pppm:density_fft");
  memory->create(greensfn,nfft_both,"pppm:greensfn");
  memory->create(work1,2*nfft_both,"pppm:work1");
  memory->create(work2,2*nfft_both,"pppm:work2");
  memory->create(vg,nfft_both,6,"pppm:vg");

  memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm:fkx");
  memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm:fky");
  memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm:fkz");

  memory->create(buf1,nbuf,"pppm:buf1");
  memory->create(buf2,nbuf,"pppm:buf2");

  // summation coeffs

  memory->create(gf_b,order,"pppm:gf_b");
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decompostion
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   0,0,&tmp);

  fft2 = new FFT3d(lmp,world,nx_pppm,ny_pppm,nz_pppm,
		   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
		   0,0,&tmp);

  remap = new Remap(lmp,world,
		    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
		    nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
		    1,0,0,FFT_PRECISION);
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMGPU::deallocate()
{
  destroy_3d_offset(density_brick_gpu,nzlo_out,nylo_out);
  destroy_3d_offset(vd_brick,nzlo_out,nylo_out);

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  memory->destroy1d_offset(fkx,nxlo_fft);
  memory->destroy1d_offset(fky,nylo_fft);
  memory->destroy1d_offset(fkz,nzlo_fft);

  memory->destroy(buf1);
  memory->destroy(buf2);

  memory->destroy(gf_b);
  memory->destroy2d_offset(rho1d,-order/2);
  memory->destroy2d_offset(rho_coeff,(1-order)/2);

  delete fft1;
  delete fft2;
  delete remap;
}


/* ----------------------------------------------------------------------
   ghost-swap to accumulate full density in brick decomposition 
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPMGPU::brick2fft()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my ghosts for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in+1; ix <= nxhi_out; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix < nxlo_in+nxlo_ghost; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix < nxlo_in; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxhi_in-nxhi_ghost+1; ix <= nxhi_in; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // pack my ghosts for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // pack my ghosts for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my real cells

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	buf1[n++] = density_brick_gpu[iz][iy][ix];

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_brick_gpu[iz][iy][ix] += buf2[n++];

  // remap from 3d brick decomposition to FFT decomposition
  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
	density_fft[n++] = density_brick_gpu[iz][iy][ix];

  remap->perform(density_fft,density_fft,work1);
}

/* ----------------------------------------------------------------------
   ghost-swap to fill ghost cells of my brick with field values
------------------------------------------------------------------------- */

void PPPMGPU::fillbrick()
{
  int i,n,ix,iy,iz;
  MPI_Request request;
  MPI_Status status;

  // pack my real cells for +z processor
  // pass data to self or +z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  int x_lo = nxlo_in * 4;
  int x_hi = nxhi_in * 4 + 1;
  for (iz = nzhi_in-nzhi_ghost+1; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[2][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz < nzlo_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }

  // pack my real cells for -z processor
  // pass data to self or -z processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_in; iz < nzlo_in+nzlo_ghost; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[2][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[2][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[2][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzhi_in+1; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }

  // pack my real cells for +y processor
  // pass data to self or +y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in-nyhi_ghost+1; iy <= nyhi_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[1][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy < nylo_in; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }

  // pack my real cells for -y processor
  // pass data to self or -y processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_in; iy < nylo_in+nylo_ghost; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[1][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[1][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[1][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nyhi_in+1; iy <= nyhi_out; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }

  // pack my real cells for +x processor
  // pass data to self or +x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  x_lo = (nxhi_in-nxhi_ghost+1) * 4;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[0][1] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  x_lo = nxlo_out * 4;
  x_hi = nxlo_in * 4;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }

  // pack my real cells for -x processor
  // pass data to self or -x processor
  // unpack and sum recv data into my ghost cells

  n = 0;
  x_lo = x_hi;
  x_hi = (nxlo_in+nxlo_ghost) * 4;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	buf1[n++] = vd_brick[iz][iy][ix];
	buf1[n++] = vd_brick[iz][iy][ix+1];
	buf1[n++] = vd_brick[iz][iy][ix+2];
      }

  if (comm->procneigh[0][0] == me)
    for (i = 0; i < n; i++) buf2[i] = buf1[i];
  else {
    MPI_Irecv(buf2,nbuf,MPI_FFT_SCALAR,comm->procneigh[0][1],0,world,&request);
    MPI_Send(buf1,n,MPI_FFT_SCALAR,comm->procneigh[0][0],0,world);
    MPI_Wait(&request,&status);
  }

  n = 0;
  x_lo = (nxhi_in + 1) * 4;
  x_hi = nxhi_out * 4 + 1;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = x_lo; ix < x_hi; ix+=4) {
	vd_brick[iz][iy][ix] = buf2[n++];
	vd_brick[iz][iy][ix+1] = buf2[n++];
	vd_brick[iz][iy][ix+2] = buf2[n++];
      }
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver 
------------------------------------------------------------------------- */

void PPPMGPU::poisson()
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k) 

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }
 
  fft1->compute(work1,work1,1);

  // if requested, compute energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
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

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fkx[i]*work1[n+1];
	work2[n+1] = -fkx[i]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  int x_hi = nxhi_in * 4 + 3;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4; i < x_hi; i+=4) {
	vd_brick[k][j][i] = work2[n];
	n += 2;
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fky[j]*work1[n+1];
	work2[n+1] = -fky[j]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4 + 1; i < x_hi; i+=4) {
	vd_brick[k][j][i] = work2[n];
	n += 2;
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	work2[n] = fkz[k]*work1[n+1];
	work2[n+1] = -fkz[k]*work1[n];
	n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4 + 2; i < x_hi; i+=4) {
	vd_brick[k][j][i] = work2[n];
	n += 2;
      }
}

/* ----------------------------------------------------------------------
   Create array using offsets from pinned memory allocation
------------------------------------------------------------------------- */

FFT_SCALAR ***PPPMGPU::create_3d_offset(int n1lo, int n1hi, int n2lo, int n2hi,
				     int n3lo, int n3hi, const char *name,
				     FFT_SCALAR *data, int vec_length)
{
  int i,j;
  int n1 = n1hi - n1lo + 1;
  int n2 = n2hi - n2lo + 1;
  int n3 = n3hi - n3lo + 1;

  FFT_SCALAR **plane = (FFT_SCALAR **)
    memory->smalloc(n1*n2*sizeof(FFT_SCALAR *),name);
  FFT_SCALAR ***array = (FFT_SCALAR ***)
    memory->smalloc(n1*sizeof(FFT_SCALAR **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3*vec_length;
    }
  }

  for (i = 0; i < n1*n2; i++) array[0][i] -= n3lo*vec_length;
  for (i = 0; i < n1; i++) array[i] -= n2lo;
  return array-n1lo;
}

/* ----------------------------------------------------------------------
   3d memory offsets
------------------------------------------------------------------------- */

void PPPMGPU::destroy_3d_offset(FFT_SCALAR ***array, int n1_offset,
				 int n2_offset)
{
  if (array == NULL) return;
  memory->sfree(&array[n1_offset][n2_offset]);
  memory->sfree(array + n1_offset);
}


/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double PPPMGPU::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * 
    (nzhi_out-nzlo_out+1);
  bytes += 4 * nbrick * sizeof(FFT_SCALAR);
  bytes += 6 * nfft_both * sizeof(double);
  bytes += nfft_both * sizeof(double);
  bytes += nfft_both*5 * sizeof(FFT_SCALAR);
  bytes += 2 * nbuf * sizeof(double);

  if (peratom_allocate_flag) {
    bytes += 7 * nbrick * sizeof(FFT_SCALAR);
    bytes += 2 * nbuf_peratom * sizeof(FFT_SCALAR);
  }

  return bytes + PPPM_GPU_API(bytes)();
}

/* ----------------------------------------------------------------------
   perform and time the 4 FFTs required for N timesteps
------------------------------------------------------------------------- */

void PPPMGPU::timing(int n, double &time3d, double &time1d) {
  if (im_real_space) {
    time3d = 0.0;
    time1d = 0.0;
    return;
  }
  PPPM::timing(n,time3d,time1d);
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void PPPMGPU::setup()
{
  if (im_real_space) return;
  PPPM::setup();
} 
