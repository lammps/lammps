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
   Contributing authors: Xiaowng Zhou (Sandia)
------------------------------------------------------------------------- */

#include "pair_eam_he_gpu.h"

#include "lammps_gpu.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;

/* ---------------------------------------------------------------------- */

PairEAMHEGPU::PairEAMHEGPU(LAMMPS *lmp) : PairEAMGPU(lmp)
{
  fileformat = FS;
  one_coeff = 1;
  he_flag = 1;

  gpu_init_fn = eam_he_gpu_init;
  gpu_clear_fn = eam_he_gpu_clear;
  gpu_compute_n_fn = eam_he_gpu_compute_n;
  gpu_compute_fn = eam_he_gpu_compute;
  gpu_compute_force_fn = eam_he_gpu_compute_force;
  gpu_bytes_fn = eam_he_gpu_bytes;
}
