/***************************************************************************
                                 hippo_ext.cpp
                             -------------------
                           Trung Dac Nguyen (Northwestern)

  Functions for LAMMPS access to hippo acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_hippo.h"

using namespace std;
using namespace LAMMPS_AL;

static Hippo<PRECISION,ACC_PRECISION> HIPPOMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int hippo_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp, const int *host_amtype2class,
                    const double *host_special_hal,
                    const double *host_special_repel,
                    const double *host_special_disp,
                    const double *host_special_mpole,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const double *host_csix, const double *host_adisp,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double polar_dscale, const double polar_uscale,
                    int& tep_size) {
  HIPPOMF.clear();
  gpu_mode=HIPPOMF.device->gpu_mode();
  double gpu_split=HIPPOMF.device->particle_split();
  int first_gpu=HIPPOMF.device->first_device();
  int last_gpu=HIPPOMF.device->last_device();
  int world_me=HIPPOMF.device->world_me();
  int gpu_rank=HIPPOMF.device->gpu_rank();
  int procs_per_gpu=HIPPOMF.device->procs_per_gpu();

  tep_size=sizeof(ACC_PRECISION); // tep_size=sizeof(PRECISION);

  HIPPOMF.device->init_message(screen,"HIPPO",first_gpu,last_gpu);

  bool message=false;
  if (HIPPOMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing GPU and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=HIPPOMF.init(ntypes, max_amtype, max_amclass,
                          host_pdamp, host_thole, host_dirdamp,
                          host_amtype2class, host_special_hal,
                          host_special_repel, host_special_disp,
                          host_special_mpole, host_special_polar_wscale,
                          host_special_polar_piscale, host_special_polar_pscale,
                          host_csix, host_adisp, nlocal, nall, max_nbors,
                          maxspecial, maxspecial15, cell_size, gpu_split,
                          screen, polar_dscale, polar_uscale);

  HIPPOMF.device->world_barrier();
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
      init_ok=HIPPOMF.init(ntypes, max_amtype, max_amclass,
                            host_pdamp, host_thole, host_dirdamp,
                            host_amtype2class, host_special_hal,
                            host_special_repel, host_special_disp,
                            host_special_mpole, host_special_polar_wscale,
                            host_special_polar_piscale, host_special_polar_pscale,
                            host_csix, host_adisp, nlocal, nall, max_nbors,
                            maxspecial, maxspecial15, cell_size, gpu_split,
                            screen, polar_dscale, polar_uscale);

    HIPPOMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    HIPPOMF.estimate_gpu_overhead();
  return init_ok;
}

void hippo_gpu_clear() {
  HIPPOMF.clear();
}

int** hippo_gpu_compute_dispersion_real(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, const double aewald, const double off2,
                           double *host_q, double *boxlo, double *prd) {                             
  return HIPPOMF.compute_dispersion_real(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, sublo, subhi,
                          tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, off2, host_q, boxlo, prd);
}

int** hippo_gpu_compute_multipole_real(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, const double aewald, const double felec, const double off2,
                           double *host_q, double *boxlo, double *prd, void **tep_ptr) {
  return HIPPOMF.compute_multipole_real(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, sublo, subhi,
                          tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, felec, off2, host_q, boxlo, prd, tep_ptr);
}

int** hippo_gpu_compute_udirect2b(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double **host_uind, double **host_uinp,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success,  const double aewald, const double off2, double *host_q,
                           double *boxlo, double *prd, void **fieldp_ptr) {
  return HIPPOMF.compute_udirect2b(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, host_uind, host_uinp,
                          sublo, subhi, tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, off2, host_q, boxlo, prd, fieldp_ptr);
}

int** hippo_gpu_compute_umutual2b(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double **host_uind, double **host_uinp,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, const double aewald, const double off2, double *host_q,
                           double *boxlo, double *prd, void **fieldp_ptr) {
  return HIPPOMF.compute_umutual2b(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, host_uind, host_uinp,
                          sublo, subhi, tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, off2, host_q, boxlo, prd, fieldp_ptr);
}

int** hippo_gpu_compute_polar_real(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup,
                           double **host_rpole, double **host_uind, double **host_uinp,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, const double aewald, const double felec, const double off2,
                           double *host_q, double *boxlo, double *prd, void **tep_ptr) {
  return HIPPOMF.compute_polar_real(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, host_uind, host_uinp,
                          sublo, subhi, tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, felec, off2, host_q, boxlo, prd, tep_ptr);
}

double hippo_gpu_bytes() {
  return HIPPOMF.host_memory_usage();
}
