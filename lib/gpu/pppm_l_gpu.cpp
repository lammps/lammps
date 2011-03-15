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
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#include <iostream>
#include <cassert>
#include <math.h>

#include "pppm_gpu_memory.h"

using namespace std;

static PPPMGPUMemory<PRECISION,ACC_PRECISION> PPPMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
PRECISION * pppm_gpu_init(const int nlocal, const int nall, FILE *screen,
                          const int order, const int nxlo_out, 
                          const int nylo_out, const int nzlo_out,
                          const int nxhi_out, const int nyhi_out,
                          const int nzhi_out, double **rho_coeff,
                          PRECISION **vdx_brick, PRECISION **vdy_brick,
                          PRECISION **vdz_brick, bool &success) {
  PPPMF.clear();
  int first_gpu=PPPMF.device->first_device();
  int last_gpu=PPPMF.device->last_device();
  int world_me=PPPMF.device->world_me();
  int gpu_rank=PPPMF.device->gpu_rank();
  int procs_per_gpu=PPPMF.device->procs_per_gpu();

  PPPMF.device->init_message(screen,"pppm",first_gpu,last_gpu);

  bool message=false;
  if (PPPMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  success=true;
  PRECISION * host_brick=NULL;
  if (world_me==0) {
    host_brick=PPPMF.init(nlocal,nall,screen,order,nxlo_out,nylo_out,
                                  nzlo_out,nxhi_out,nyhi_out,nzhi_out,rho_coeff,
                                  vdx_brick, vdy_brick, vdz_brick, success);
    if (!success)
      return host_brick;
  }

  PPPMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing GPU %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing GPUs %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0) {
      host_brick=PPPMF.init(nlocal,nall,screen,order,nxlo_out,nylo_out,
                            nzlo_out,nxhi_out,nyhi_out,nzhi_out,rho_coeff,
                            vdx_brick, vdy_brick, vdz_brick, success);
      if (!success)
        return host_brick;
    }
    PPPMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
  return host_brick;
}

void pppm_gpu_clear() {
  PPPMF.clear();
}

int pppm_gpu_spread(const int ago, const int nlocal, const int nall,
                    double **host_x, int *host_type, bool &success,
                    double *host_q, double *boxlo, const double delxinv,
                    const double delyinv, const double delzinv) {
  return PPPMF.spread(ago,nlocal,nall,host_x,host_type,success,host_q,boxlo,
                      delxinv,delyinv,delzinv);
}

void pppm_gpu_interp(const PRECISION qqrd2e_scale) {
  return PPPMF.interp(qqrd2e_scale);
}

double pppm_gpu_bytes() {
  return PPPMF.host_memory_usage();
}

PRECISION * pppm_gpu_force() { return PPPMF.force_temp.begin(); }

