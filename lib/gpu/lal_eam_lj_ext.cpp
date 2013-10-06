// **************************************************************************
//                               lal_eam_ext.cpp
//                             -------------------
//                     W. Michael Brown, Trung Dac Nguyen (ORNL)
//
//  Class for acceleration of the eam pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov nguyentd@ornl.gov
// ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_eam_lj.h"

using namespace std;
using namespace LAMMPS_AL;

static EAMLJ<PRECISION,ACC_PRECISION> EAMLJMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int eam_lj_gpu_init(const int ntypes, double host_cutforcesq, 
                 int **host_type2rhor, int **host_type2z2r, int *host_type2frho,
                 double ***host_rhor_spline, double ***host_z2r_spline,
                 double ***host_frho_spline, double **host_cutljsq, 
                 double **host_lj1, double **host_lj2, double **host_lj3, 
                 double **host_lj4, double **offset, double *special_lj,
                 double rdr, double rdrho, int nrhor, 
                 int nrho, int nz2r, int nfrho, int nr, 
                 const int nlocal, const int nall, const int max_nbors, 
                 const int maxspecial, const double cell_size, 
                 int &gpu_mode, FILE *screen, int &fp_size) {
  EAMLJMF.clear();
  gpu_mode=EAMLJMF.device->gpu_mode();
  double gpu_split=EAMLJMF.device->particle_split();
  int first_gpu=EAMLJMF.device->first_device();
  int last_gpu=EAMLJMF.device->last_device();
  int world_me=EAMLJMF.device->world_me();
  int gpu_rank=EAMLJMF.device->gpu_rank();
  int procs_per_gpu=EAMLJMF.device->procs_per_gpu();

  // disable host/device split for now
  if (gpu_split != 1.0) 
    return -8;
    
  fp_size=sizeof(PRECISION);
    
  EAMLJMF.device->init_message(screen,"eam_lj",first_gpu,last_gpu);

  bool message=false;
  if (EAMLJMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=EAMLJMF.init(ntypes, host_cutforcesq,
                    host_type2rhor, host_type2z2r, host_type2frho,
                    host_rhor_spline, host_z2r_spline,
                    host_frho_spline, host_cutljsq, host_lj1,
                    host_lj2, host_lj3, host_lj4,
                    offset, special_lj, rdr, rdrho, nrhor, 
                    nrho, nz2r, nfrho, nr, 
                    nlocal, nall, 300, maxspecial,
                    cell_size, gpu_split, screen);

  EAMLJMF.device->world_barrier();
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
      init_ok=EAMLJMF.init(ntypes, host_cutforcesq,
                    host_type2rhor, host_type2z2r, host_type2frho,
                    host_rhor_spline, host_z2r_spline,
                    host_frho_spline, host_cutljsq, host_lj1,
                    host_lj2, host_lj3, host_lj4,
                    offset, special_lj, rdr, rdrho, nrhor, 
                    nrho, nz2r, nfrho, nr, 
                    nlocal, nall, 300, maxspecial,
                    cell_size, gpu_split, screen);

    EAMLJMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    EAMLJMF.estimate_gpu_overhead();
  return init_ok;
}

void eam_lj_gpu_clear() {
  EAMLJMF.clear();
}

int ** eam_lj_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,  const double cpu_time,
                         bool &success, int &inum, void **fp_ptr) {
  return EAMLJMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        inum, fp_ptr);
}  

void eam_lj_gpu_compute(const int ago, const int inum_full, const int nlocal, 
                     const int nall, double **host_x, int *host_type, 
                     int *ilist, int *numj, int **firstneigh, const bool eflag,
                     const bool vflag, const bool eatom, const bool vatom,
                     int &host_start, const double cpu_time, bool &success,
                     void **fp_ptr) {
  EAMLJMF.compute(ago,inum_full,nlocal,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                fp_ptr);
}

void eam_lj_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom) {
  EAMLJMF.compute2(ilist, eflag, vflag, eatom, vatom);
}


double eam_lj_gpu_bytes() {
  return EAMLJMF.host_memory_usage();
}


