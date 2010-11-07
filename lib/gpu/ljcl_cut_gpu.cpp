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

#include "ljcl_cut_gpu_memory.h"

using namespace std;

static LJCL_GPU_Memory<PRECISION,ACC_PRECISION> LJCLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
bool ljcl_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                   double **host_lj2, double **host_lj3, double **host_lj4,
                   double **offset, double *special_lj, const int inum,
                   const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen,
                   double **host_cut_ljsq, double host_cut_coulsq,
                   double *host_special_coul, const double qqrd2e,
                   const double g_ewald) {
  LJCLMF.clear();
  gpu_mode=LJCLMF.device->gpu_mode();
  double gpu_split=LJCLMF.device->particle_split();
  int first_gpu=LJCLMF.device->first_device();
  int last_gpu=LJCLMF.device->last_device();
  int world_me=LJCLMF.device->world_me();
  int gpu_rank=LJCLMF.device->gpu_rank();
  int procs_per_gpu=LJCLMF.device->procs_per_gpu();

  LJCLMF.device->init_message(screen,"lj/cut/coul/long",first_gpu,last_gpu);

  bool message=false;
  if (world_me==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  if (world_me==0) {
    bool init_ok=LJCLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                            host_lj4, offset, special_lj, inum, nall, 300,
                            maxspecial, cell_size, gpu_split, screen,
                            host_cut_ljsq, host_cut_coulsq, host_special_coul,
                            qqrd2e,g_ewald);
    if (!init_ok)
      return false;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing GPU %d on core %d...",gpu_rank,i);
      else
        fprintf(screen,"Initializing GPUs %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0) {
      bool init_ok=LJCLMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                              host_lj4, offset, special_lj, inum, nall, 300,
                              maxspecial, cell_size, gpu_split,
			      screen, host_cut_ljsq, host_cut_coulsq,
                              host_special_coul, qqrd2e, g_ewald);
      if (!init_ok)
        return false;
    }
    MPI_Barrier(LJCLMF.device->gpu_comm);
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");
  return true;
}

void ljcl_gpu_clear() {
  LJCLMF.clear();
}

int * ljcl_gpu_compute_n(const int timestep, const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *boxlo, double *boxhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         const double cpu_time, bool &success, double *host_q) {
  return LJCLMF.compute(timestep, ago, inum_full, nall, host_x, host_type, boxlo,
                        boxhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, cpu_time, success, host_q);
}  
			
void ljcl_gpu_compute(const int timestep, const int ago, const int inum_full,
	 	     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
		     const bool eflag, const bool vflag, const bool eatom,
                     const bool vatom, int &host_start, const double cpu_time,
                     bool &success, double *host_q) {
  LJCLMF.compute(timestep,ago,inum_full,nall,host_x,host_type,ilist,numj,
                firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success,
                host_q);
}

double ljcl_gpu_bytes() {
  return LJCLMF.host_memory_usage();
}


