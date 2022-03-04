/***************************************************************************
                              vashishta_ext.cpp
                             -------------------
                            Anders Hafreager (UiO)

  Class for acceleration of the vashishta pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Mon June 12, 2017
    email                : andershaf@gmail.com
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_vashishta.h"
using namespace LAMMPS_AL;

static Vashishta<PRECISION,ACC_PRECISION> VashishtaMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int vashishta_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                const double cell_size, int &gpu_mode, FILE *screen,
                int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                const double* cutsq, const double* r0,
                const double* gamma, const double* eta,
                const double* lam1inv, const double* lam4inv,
                const double* zizj, const double* mbigd,
                const double* dvrc, const double* big6w,
                const double* heta, const double* bigh,
                const double* bigw, const double* c0,
                const double* costheta, const double* bigb,
                const double* big2b, const double* bigc) {
  VashishtaMF.clear();
  gpu_mode=VashishtaMF.device->gpu_mode();
  double gpu_split=VashishtaMF.device->particle_split();
  int first_gpu=VashishtaMF.device->first_device();
  int last_gpu=VashishtaMF.device->last_device();
  int world_me=VashishtaMF.device->world_me();
  int gpu_rank=VashishtaMF.device->gpu_rank();
  int procs_per_gpu=VashishtaMF.device->procs_per_gpu();

  // disable host/device split for now
  if (gpu_split != 1.0)
    return -8;

  VashishtaMF.device->init_message(screen,"vashishta/gpu",first_gpu,last_gpu);

  bool message=false;
  if (VashishtaMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=VashishtaMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split, screen,
                      host_map, nelements, host_elem2param, nparams,
                      cutsq, r0, gamma, eta, lam1inv,
                      lam4inv, zizj, mbigd, dvrc, big6w, heta, bigh, bigw,
                      c0, costheta, bigb, big2b, bigc);

  VashishtaMF.device->world_barrier();
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
      init_ok=VashishtaMF.init(ntypes, inum, nall, max_nbors, cell_size, gpu_split, screen,
                        host_map, nelements, host_elem2param, nparams,
                        cutsq, r0, gamma, eta, lam1inv,
                        lam4inv, zizj, mbigd, dvrc, big6w, heta, bigh, bigw,
                        c0, costheta, bigb, big2b, bigc);

    VashishtaMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    VashishtaMF.estimate_gpu_overhead();

  return init_ok;
}

void vashishta_gpu_clear() {
  VashishtaMF.clear();
}

int ** vashishta_gpu_compute_n(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success) {
  return VashishtaMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success);
}

void vashishta_gpu_compute(const int ago, const int nlocal, const int nall,
                    const int nlist, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom,
                    int &host_start, const double cpu_time, bool &success) {
  VashishtaMF.compute(ago,nlocal,nall,nlist,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double vashishta_gpu_bytes() {
  return VashishtaMF.host_memory_usage();
}


