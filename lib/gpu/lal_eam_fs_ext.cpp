/***************************************************************************
                                 eam_fs_ext.cpp
                             -------------------
                   Trung Dac Nguyen, W. Michael Brown (ORNL)

  Functions for LAMMPS access to buck acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_eam.h"

using namespace std;
using namespace LAMMPS_AL;

static EAM<PRECISION,ACC_PRECISION> EAMFSMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int eam_fs_gpu_init(const int ntypes, double host_cutforcesq,
                 int **host_type2rhor, int **host_type2z2r, int *host_type2frho,
                 double ***host_rhor_spline, double ***host_z2r_spline,
                 double ***host_frho_spline, double** host_cutsq,
                 double rdr, double rdrho, double rhomax, int nrhor,
                 int nrho, int nz2r, int nfrho, int nr,
                 const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen, int &fp_size) {
  EAMFSMF.clear();
  gpu_mode=EAMFSMF.device->gpu_mode();
  double gpu_split=EAMFSMF.device->particle_split();
  int first_gpu=EAMFSMF.device->first_device();
  int last_gpu=EAMFSMF.device->last_device();
  int world_me=EAMFSMF.device->world_me();
  int gpu_rank=EAMFSMF.device->gpu_rank();
  int procs_per_gpu=EAMFSMF.device->procs_per_gpu();

  // disable host/device split for now
  if (gpu_split != 1.0)
    return -8;

  fp_size=sizeof(PRECISION);

  EAMFSMF.device->init_message(screen,"eam/fs",first_gpu,last_gpu);

  bool message=false;
  if (EAMFSMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=EAMFSMF.init(ntypes, host_cutforcesq, host_type2rhor, host_type2z2r,
                       host_type2frho, host_rhor_spline, host_z2r_spline,
                       host_frho_spline, host_cutsq, rdr, rdrho, rhomax, nrhor, nrho, nz2r,
                       nfrho, nr, nlocal, nall, max_nbors, maxspecial, cell_size,
                       gpu_split, screen);

  EAMFSMF.device->world_barrier();
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
      init_ok=EAMFSMF.init(ntypes, host_cutforcesq, host_type2rhor, host_type2z2r,
                         host_type2frho, host_rhor_spline, host_z2r_spline,
                         host_frho_spline, host_cutsq, rdr, rdrho, rhomax, nrhor, nrho,
                         nz2r, nfrho, nr, nlocal, nall, max_nbors, maxspecial,
                         cell_size, gpu_split, screen);

    EAMFSMF.device->serialize_init();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    EAMFSMF.estimate_gpu_overhead();
  return init_ok;
}

void eam_fs_gpu_clear() {
  EAMFSMF.clear();
}

int ** eam_fs_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,  const double cpu_time,
                         bool &success, int &inum, void **fp_ptr) {
  return EAMFSMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                       subhi, tag, nspecial, special, eflag, vflag, eatom,
                       vatom, host_start, ilist, jnum, cpu_time, success,
                       inum, fp_ptr);
}

void eam_fs_gpu_compute(const int ago, const int inum_full, const int nlocal,
                     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag,
                     const bool vflag, const bool eatom, const bool vatom,
                     int &host_start, const double cpu_time, bool &success,
                     void **fp_ptr) {
  EAMFSMF.compute(ago,inum_full,nlocal,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                fp_ptr);
}

void eam_fs_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom) {
  EAMFSMF.compute2(ilist, eflag, vflag, eatom, vatom);
}


double eam_fs_gpu_bytes() {
  return EAMFSMF.host_memory_usage();
}


