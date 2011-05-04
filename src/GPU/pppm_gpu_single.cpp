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
   Contributing authors: Mike Brown (ORNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_gpu_single.h"
#include "lmptype.h"
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
#include "memory.h"
#include "error.h"
#include "gpu_extra.h"

#define grdtyp float

// External functions from cuda library for atom decomposition

grdtyp* pppm_gpu_init_f(const int nlocal, const int nall, FILE *screen,
		        const int order, const int nxlo_out, 
			const int nylo_out, const int nzlo_out,
			const int nxhi_out, const int nyhi_out,
			const int nzhi_out, double **rho_coeff,
			grdtyp **_vd_brick, const double slab_volfactor,
			const int nx_pppm, const int ny_pppm,
			const int nz_pppm, int &success);
void pppm_gpu_clear_f(const double poisson_time);
int pppm_gpu_spread_f(const int ago, const int nlocal, const int nall,
                      double **host_x, int *host_type, bool &success,
                      double *host_q, double *boxlo, const double delxinv,
                      const double delyinv, const double delzinv);
void pppm_gpu_interp_f(const grdtyp qqrd2e_scale);
double pppm_gpu_bytes_f();

using namespace LAMMPS_NS;

#define MAXORDER 7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PPPMGPUSingle::PPPMGPUSingle(LAMMPS *lmp, int narg, char **arg) :
  PPPMGPU<grdtyp>(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

PPPMGPUSingle::~PPPMGPUSingle()
{
  pppm_gpu_clear_f(poisson_time);
}

/* ----------------------------------------------------------------------
   called once before run 
------------------------------------------------------------------------- */

void PPPMGPUSingle::init()
{
  base_init();

  if (order>8)
    error->all("Cannot use order greater than 8 with pppm/gpu.");
  pppm_gpu_clear_f(poisson_time);

  int success;
  grdtyp *data, *h_brick;
  h_brick = pppm_gpu_init_f(atom->nlocal, atom->nlocal+atom->nghost, screen,
			    order, nxlo_out, nylo_out, nzlo_out, nxhi_out,
			    nyhi_out, nzhi_out, rho_coeff, &data, 
			    slab_volfactor,nx_pppm,ny_pppm,nz_pppm,success);

  GPU_EXTRA::check_flag(success,error,world);

  density_brick =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
		     nxlo_out,nxhi_out,"pppm:density_brick",h_brick,1);
  vd_brick =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
		     nxlo_out,nxhi_out,"pppm:vd_brick",data,4);

  poisson_time=0;
}

/* ----------------------------------------------------------------------
   compute the PPPMGPU long-range force, energy, virial 
------------------------------------------------------------------------- */

void PPPMGPUSingle::compute(int eflag, int vflag)
{
  bool success = true;
  int flag=pppm_gpu_spread_f(neighbor->ago, atom->nlocal, atom->nlocal + 
			     atom->nghost, atom->x, atom->type, success,
			     atom->q, domain->boxlo, delxinv, delyinv,
			     delzinv);
  if (!success)
    error->one("Out of memory on GPGPU");
  if (flag != 0)
    error->one("Out of range atoms - cannot compute PPPM");

  int i;

  // convert atoms from box to lamda coords
  
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  energy = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  double t3=MPI_Wtime();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  
  poisson(eflag,vflag);

  // all procs communicate E-field values to fill ghost cells
  //   surrounding their 3d bricks

  fillbrick();

  poisson_time+=MPI_Wtime()-t3;

  // calculate the force on my particles

  grdtyp qqrd2e_scale=qqrd2e*scale;
  pppm_gpu_interp_f(qqrd2e_scale);

  // sum energy across procs and add in volume-dependent term

  if (eflag) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
   
    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/1.772453851 +
      0.5*PI*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qqrd2e*scale;
  }

  // sum virial across procs

  if (vflag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qqrd2e*scale*volume*virial_all[i];
  }

  // 2d slab correction

  if (slabflag) slabcorr(eflag);

  // convert atoms back from lamda to box coords
  
  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double PPPMGPUSingle::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * 
    (nzhi_out-nzlo_out+1);
  bytes += 4 * nbrick * sizeof(grdtyp);
  bytes += 6 * nfft_both * sizeof(double);
  bytes += nfft_both*6 * sizeof(double);
  bytes += 2 * nbuf * sizeof(double);
  return bytes + pppm_gpu_bytes_f();
}
