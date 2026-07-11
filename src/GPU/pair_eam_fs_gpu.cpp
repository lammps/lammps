/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Trung Dac Nguyen (ORNL), W. Michael Brown (ORNL)
------------------------------------------------------------------------- */

#include "pair_eam_fs_gpu.h"

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int eam_fs_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                    int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                    double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                    double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r,
                    int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                    int &fp_size);
void eam_fs_gpu_clear();
int **eam_fs_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start, int **ilist,
                           int **jnum, const double cpu_time, bool &success, int &inum,
                           void **fp_ptr, double *prd, int *periodicity);
void eam_fs_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                        const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                        int &host_start, const double cpu_time, bool &success, void **fp_ptr);
void eam_fs_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom);
double eam_fs_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairEAMFSGPU::PairEAMFSGPU(LAMMPS *lmp) : PairEAMGPU(lmp)
{
  fileformat = FS;
  one_coeff = 1;

  gpu_init_fn = eam_fs_gpu_init;
  gpu_clear_fn = eam_fs_gpu_clear;
  gpu_compute_n_fn = eam_fs_gpu_compute_n;
  gpu_compute_fn = eam_fs_gpu_compute;
  gpu_compute_force_fn = eam_fs_gpu_compute_force;
  gpu_bytes_fn = eam_fs_gpu_bytes;
}
