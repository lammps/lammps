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
int amoeba_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
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
                    const double polar_dscale, const double polar_uscale) {
  AMOEBAMF.clear();
  gpu_mode=AMOEBAMF.device->gpu_mode();
  double gpu_split=AMOEBAMF.device->particle_split();
  int first_gpu=AMOEBAMF.device->first_device();
  int last_gpu=AMOEBAMF.device->last_device();
  int world_me=AMOEBAMF.device->world_me();
  int gpu_rank=AMOEBAMF.device->gpu_rank();
  int procs_per_gpu=AMOEBAMF.device->procs_per_gpu();

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
    init_ok=AMOEBAMF.init(ntypes, max_amtype, max_amclass,
                          host_pdamp, host_thole, host_dirdamp,
                          host_amtype2class, host_special_hal,
                          host_special_repel, host_special_disp,
                          host_special_mpole, host_special_polar_wscale,
                          host_special_polar_piscale, host_special_polar_pscale,
                          host_csix, host_adisp, nlocal, nall, max_nbors,
                          maxspecial, maxspecial15, cell_size, gpu_split,
                          screen, polar_dscale, polar_uscale);

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
      init_ok=AMOEBAMF.init(ntypes, max_amtype, max_amclass,
                            host_pdamp, host_thole, host_dirdamp,
                            host_amtype2class, host_special_hal,
                            host_special_repel, host_special_disp,
                            host_special_mpole, host_special_polar_wscale,
                            host_special_polar_piscale, host_special_polar_pscale,
                            host_csix, host_adisp, nlocal, nall, max_nbors,
                            maxspecial, maxspecial15, cell_size, gpu_split,
                            screen, polar_dscale, polar_uscale);

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

int** amoeba_gpu_precompute(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, int *host_amtype,
                            int *host_amgroup, double **host_rpole, double ** /*host_uind*/,
                            double ** /*host_uinp*/, double * /*host_pval*/,
                            double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special,
                            int *nspecial15, tagint **special15,
                            const bool eflag_in, const bool vflag_in,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo, double *prd) {
  return AMOEBAMF.precompute(ago, inum_full, nall, host_x, host_type,
                             host_amtype, host_amgroup, host_rpole,
                             nullptr, nullptr, nullptr, sublo, subhi, tag,
                             nspecial, special, nspecial15, special15,
                             eflag_in, vflag_in, eatom, vatom,
                             host_start, ilist, jnum, cpu_time,
                             success, host_q, boxlo, prd);
}


void amoeba_gpu_compute_multipole_real(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, const double aewald, const double felec, const double off2,
                           double *host_q, double *boxlo, double *prd, void **tep_ptr) {
  AMOEBAMF.compute_multipole_real(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole, nullptr, sublo, subhi,
                          tag, nspecial, special, nspecial15, special15,
                          eflag, vflag, eatom, vatom, host_start, ilist, jnum,
                          cpu_time, success, aewald, felec, off2, host_q, boxlo, prd, tep_ptr);
}

void amoeba_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                           double **host_uind, double **host_uinp,
                           const double aewald, const double off2, void **fieldp_ptr) {
  AMOEBAMF.compute_udirect2b(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, nullptr,
                             aewald, off2, fieldp_ptr);
}

void amoeba_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                   double **host_uind, double **host_uinp,
                                   const double aewald, const double off2, void **fieldp_ptr) {
  AMOEBAMF.compute_umutual2b(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, nullptr,
                             aewald, off2, fieldp_ptr);
}

void amoeba_gpu_update_fieldp(void **fieldp_ptr) {
  AMOEBAMF.update_fieldp(fieldp_ptr);
}

void amoeba_gpu_compute_polar_real(int *host_amtype, int *host_amgroup, double **host_rpole,
                                    double **host_uind, double **host_uinp,
                                    const bool eflag_in, const bool vflag_in,
                                    const bool eatom, const bool vatom,
                                    const double aewald, const double felec, const double off2,
                                    void **tep_ptr) {
  AMOEBAMF.compute_polar_real(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, nullptr,
                              eflag_in, vflag_in, eatom, vatom, aewald, felec, off2, tep_ptr);
}

void amoeba_gpu_precompute_kspace(const int inum_full, const int bsorder,
                          double ***host_thetai1, double ***host_thetai2,
                          double ***host_thetai3, int** igrid,
                          const int nzlo_out, const int nzhi_out,
                          const int nylo_out, const int nyhi_out,
                          const int nxlo_out, const int nxhi_out) {
   AMOEBAMF.precompute_kspace(inum_full, bsorder, host_thetai1, host_thetai2, host_thetai3, igrid,
                              nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out);
}

void amoeba_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                          void **host_fdip_phi2, void **host_fdip_sum_phi) {
   AMOEBAMF.compute_fphi_uind(host_grid_brick, host_fdip_phi1,
                              host_fdip_phi2, host_fdip_sum_phi);
}

void amoeba_gpu_fphi_mpole(double ***host_grid_brick, void **host_fphi, const double felec) {
   AMOEBAMF.compute_fphi_mpole(host_grid_brick, host_fphi, felec);
}

void amoeba_setup_fft(const int numel, const int element_type) {
  AMOEBAMF.setup_fft(numel, element_type);
}

void amoeba_compute_fft1d(void* in, void* out, const int numel, const int mode) {
  AMOEBAMF.compute_fft1d(in, out, numel, mode);
}

double amoeba_gpu_bytes() {
  return AMOEBAMF.host_memory_usage();
}
