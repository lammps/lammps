/***************************************************************************
                                 amoeba_ext.cpp
                             -------------------
                           Trung Dac Nguyen (Northwestern)

  Functions for LAMMPS access to amoeba acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_amoeba.h"

using namespace std;
using namespace LAMMPS_AL;

static Amoeba<PRECISION,ACC_PRECISION> AMOEBAMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int amoeba_gpu_init(const int ntypes, const int max_amtype,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double aewald, const double felec,
                    const double off2, const double polar_dscale,
                    const double polar_uscale, int& tep_size) {
  AMOEBAMF.clear();
  gpu_mode=AMOEBAMF.device->gpu_mode();
  double gpu_split=AMOEBAMF.device->particle_split();
  int first_gpu=AMOEBAMF.device->first_device();
  int last_gpu=AMOEBAMF.device->last_device();
  int world_me=AMOEBAMF.device->world_me();
  int gpu_rank=AMOEBAMF.device->gpu_rank();
  int procs_per_gpu=AMOEBAMF.device->procs_per_gpu();

  tep_size=sizeof(PRECISION);

  AMOEBAMF.device->init_message(screen,"amoeba",first_gpu,last_gpu);

  bool message=false;
  if (AMOEBAMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=AMOEBAMF.init(ntypes, max_amtype, host_pdamp, host_thole, host_dirdamp,
                          host_special_polar_wscale, host_special_polar_piscale,
                          host_special_polar_pscale, nlocal, nall, max_nbors,
                          maxspecial, maxspecial15, cell_size, gpu_split, screen,
                          aewald, felec, off2, polar_dscale, polar_uscale);

  AMOEBAMF.device->world_barrier();
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
      init_ok=AMOEBAMF.init(ntypes, max_amtype, host_pdamp, host_thole, host_dirdamp,
                            host_special_polar_wscale, host_special_polar_piscale,
                            host_special_polar_pscale, nlocal, nall, max_nbors,
                            maxspecial, maxspecial15, cell_size, gpu_split, screen,
                            aewald, felec, off2, polar_dscale, polar_uscale);

    AMOEBAMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    AMOEBAMF.estimate_gpu_overhead();
  return init_ok;
}

void amoeba_gpu_clear() {
  AMOEBAMF.clear();
}

int** amoeba_gpu_compute_polar_real(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup,
                           double **host_rpole, double **host_uind, double **host_uinp,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd, void **tep_ptr) {
  return AMOEBAMF.compute_polar_real(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, host_uind, host_uinp,
                          sublo, subhi, tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom,
                          vatom, host_start, ilist, jnum, cpu_time, success,
                          host_q, boxlo, prd, tep_ptr);
}

int** amoeba_gpu_compute_udirect2b(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd, void **fieldp_ptr) {
  return AMOEBAMF.compute_udirect2b(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, sublo, subhi,
                          tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, host_q, boxlo, prd, fieldp_ptr);
}

double amoeba_gpu_bytes() {
  return AMOEBAMF.host_memory_usage();
}
