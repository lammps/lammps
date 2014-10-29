/***************************************************************************
                                 pppm_ext.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Functions for LAMMPS access to PPPM acceleration routines

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_pppm.h"

using namespace std;
using namespace LAMMPS_AL;

static PPPM<PRECISION,ACC_PRECISION,float,_lgpu_float4> PPPMF;
static PPPM<PRECISION,ACC_PRECISION,double,_lgpu_double4> PPPMD;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
template <class grdtyp, class memtyp>
grdtyp * pppm_gpu_init(memtyp &pppm, const int nlocal, const int nall,
                       FILE *screen, const int order, const int nxlo_out, 
                       const int nylo_out, const int nzlo_out,
                       const int nxhi_out, const int nyhi_out,
                       const int nzhi_out, grdtyp **rho_coeff,
                       grdtyp **vd_brick, const double slab_volfactor,
                       const int nx_pppm, const int ny_pppm, const int nz_pppm,
                       const bool split, int &success) {
  pppm.clear(0.0);
  int first_gpu=pppm.device->first_device();
  int last_gpu=pppm.device->last_device();
  int world_me=pppm.device->world_me();
  int gpu_rank=pppm.device->gpu_rank();
  int procs_per_gpu=pppm.device->procs_per_gpu();

  pppm.device->init_message(screen,"pppm",first_gpu,last_gpu);

  bool message=false;
  if (pppm.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  success=0;
  grdtyp * host_brick=NULL;
  if (world_me==0)
    host_brick=pppm.init(nlocal,nall,screen,order,nxlo_out,nylo_out,nzlo_out,
                         nxhi_out,nyhi_out,nzhi_out,rho_coeff,vd_brick,
                         slab_volfactor,nx_pppm,ny_pppm,nz_pppm,split,success);

  pppm.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      host_brick=pppm.init(nlocal,nall,screen,order,nxlo_out,nylo_out,
                           nzlo_out,nxhi_out,nyhi_out,nzhi_out,rho_coeff,
                           vd_brick,slab_volfactor,nx_pppm,ny_pppm,nz_pppm,
                           split,success);

    pppm.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
  return host_brick;
}

float * pppm_gpu_init_f(const int nlocal, const int nall, FILE *screen,
                        const int order, const int nxlo_out, 
                        const int nylo_out, const int nzlo_out,
                        const int nxhi_out, const int nyhi_out,
                        const int nzhi_out, float **rho_coeff,
                        float **vd_brick, const double slab_volfactor,
                        const int nx_pppm, const int ny_pppm, const int nz_pppm,
                        const bool split, const bool respa, int &success) {
  float *b=pppm_gpu_init(PPPMF,nlocal,nall,screen,order,nxlo_out,nylo_out,
                         nzlo_out,nxhi_out,nyhi_out,nzhi_out,rho_coeff,vd_brick,
                         slab_volfactor,nx_pppm,ny_pppm,nz_pppm,split,success);
  if (split==false && respa==false)
    PPPMF.device->set_single_precompute(&PPPMF);                         
  return b;
}

void pppm_gpu_clear_f(const double cpu_time) {
  PPPMF.clear(cpu_time);
}

int pppm_gpu_spread_f(const int ago, const int nlocal, const int nall,
                     double **host_x, int *host_type, bool &success,
                     double *host_q, double *boxlo, const double delxinv,
                     const double delyinv, const double delzinv) {
  return PPPMF.spread(ago,nlocal,nall,host_x,host_type,success,host_q,boxlo,
                      delxinv,delyinv,delzinv);
}

void pppm_gpu_interp_f(const float qqrd2e_scale) {
  PPPMF.interp(qqrd2e_scale);
}

double pppm_gpu_bytes_f() {
  return PPPMF.host_memory_usage();
}

void pppm_gpu_forces_f(double **f) {
  double etmp;
  PPPMF.atom->data_unavail();
  PPPMF.ans->get_answers(f,NULL,NULL,NULL,NULL,etmp);
}

double * pppm_gpu_init_d(const int nlocal, const int nall, FILE *screen,
                         const int order, const int nxlo_out, 
                         const int nylo_out, const int nzlo_out,
                         const int nxhi_out, const int nyhi_out,
                         const int nzhi_out, double **rho_coeff,
                         double **vd_brick, const double slab_volfactor,
                         const int nx_pppm, const int ny_pppm,
                         const int nz_pppm, const bool split, 
                         const bool respa, int &success) {
  double *b=pppm_gpu_init(PPPMD,nlocal,nall,screen,order,nxlo_out,nylo_out,
                          nzlo_out,nxhi_out,nyhi_out,nzhi_out,rho_coeff,
                          vd_brick,slab_volfactor,nx_pppm,ny_pppm,nz_pppm,
                          split,success);                        
  if (split==false && respa==false)
    PPPMD.device->set_double_precompute(&PPPMD);                         
  return b;
}

void pppm_gpu_clear_d(const double cpu_time) {
  PPPMD.clear(cpu_time);
}

int pppm_gpu_spread_d(const int ago, const int nlocal, const int nall,
                      double **host_x, int *host_type, bool &success,
                      double *host_q, double *boxlo, const double delxinv,
                      const double delyinv, const double delzinv) {
  return PPPMD.spread(ago,nlocal,nall,host_x,host_type,success,host_q,boxlo,
                      delxinv,delyinv,delzinv);
}

void pppm_gpu_interp_d(const double qqrd2e_scale) {
  PPPMD.interp(qqrd2e_scale);
}

double pppm_gpu_bytes_d() {
  return PPPMD.host_memory_usage();
}

void pppm_gpu_forces_d(double **f) {
  double etmp;
  PPPMD.atom->data_unavail();
  PPPMD.ans->get_answers(f,NULL,NULL,NULL,NULL,etmp);
}

