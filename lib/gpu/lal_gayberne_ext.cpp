/***************************************************************************
                              gayberne_ext.cpp
                             -------------------
                               W. Michael Brown

  LAMMPS Wrappers for Gay-Berne Acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_gayberne.h"

using namespace std;
using namespace LAMMPS_AL;

static GayBerne<PRECISION,ACC_PRECISION> GBMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int gb_gpu_init(const int ntypes, const double gamma,
                const double upsilon, const double mu, double **shape,
                double **well, double **cutsq, double **sigma,
                double **epsilon, double *host_lshape, int **form,
                double **host_lj1, double **host_lj2, double **host_lj3,
                double **host_lj4, double **offset, double *special_lj,
                const int inum, const int nall, const int max_nbors, 
                const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen) {
  GBMF.clear();
  gpu_mode=GBMF.device->gpu_mode();
  double gpu_split=GBMF.device->particle_split();
  int first_gpu=GBMF.device->first_device();
  int last_gpu=GBMF.device->last_device();
  int world_me=GBMF.device->world_me();
  int gpu_rank=GBMF.device->gpu_rank();
  int procs_per_gpu=GBMF.device->procs_per_gpu();

  GBMF.device->init_message(screen,"gayberne",first_gpu,last_gpu);

  bool message=false;
  if (GBMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=GBMF.init(ntypes, gamma, upsilon, mu, shape, well, cutsq, 
                      sigma, epsilon, host_lshape, form, host_lj1, 
                      host_lj2, host_lj3, host_lj4, offset, special_lj, 
                      inum, nall, max_nbors, maxspecial, cell_size, gpu_split,
                      screen);

  GBMF.device->world_barrier();
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
    if (gpu_rank==i && world_me!=0)
      init_ok=GBMF.init(ntypes, gamma, upsilon, mu, shape, well, cutsq,  sigma,
                        epsilon, host_lshape, form, host_lj1, host_lj2,
                        host_lj3, host_lj4, offset, special_lj,  inum, nall,
                        max_nbors, maxspecial, cell_size, gpu_split,  screen);

    GBMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    GBMF.estimate_gpu_overhead();
  return init_ok;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
void gb_gpu_clear() {
  GBMF.clear();
}

  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, int *tag, int **nspecial,
                int **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                int **ilist, int **numj, const double cpu_time, bool &success,
                double **host_quat);

int** gb_gpu_compute_n(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, double *sublo,
                       double *subhi, int *tag, int **nspecial, int **special,
                       const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int &host_start, int **ilist,
                       int **jnum, const double cpu_time, bool &success,
                       double **host_quat) {
  return GBMF.compute(ago, inum_full, nall, host_x, host_type, sublo, subhi, 
                      tag, nspecial, special, eflag, vflag, eatom, vatom, 
                      host_start, ilist, jnum, cpu_time, success, host_quat);
}

int * gb_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double **host_quat) {
  return GBMF.compute(ago, inum_full, nall, host_x, host_type, ilist,
                      numj, firstneigh, eflag, vflag, eatom, vatom, host_start,
                      cpu_time, success, host_quat);
}

// ---------------------------------------------------------------------------
// Return memory usage
// ---------------------------------------------------------------------------
double gb_gpu_bytes() {
  return GBMF.host_memory_usage();
}

