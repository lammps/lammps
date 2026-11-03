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

#include "lammps_gpu.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;


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
