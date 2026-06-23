// clang-format off
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
   Contributing authors: Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "ewald_gpu.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

// external functions from gpu library

void ewald_gpu_init(const int nlocal, const int nall, FILE *screen,
                         int &success);
void ewald_gpu_setup(const int kmax, const int kcount, int *kxvecs,
                          int *kyvecs, int *kzvecs, double *ug, double **eg,
                          double **vg, double *unitk, int &success);
int ewald_gpu_structure(const int ago, const int nlocal, const int nall,
                             double **host_x, int *host_type, double *host_q,
                             double *host_sfacrl, double *host_sfacim,
                             bool &success);
void ewald_gpu_compute(double *host_sfacrl_all, double *host_sfacim_all,
                            const double qscale, const int slabflag,
                            const int eflag_atom, const int vflag_atom,
                            double *host_eatom, double **host_vatom,
                            bool &success);
void ewald_gpu_clear(const double cpu_time);
double ewald_gpu_bytes();

/* ---------------------------------------------------------------------- */

EwaldGPU::EwaldGPU(LAMMPS *lmp) : Ewald(lmp), cpu_time(0.0)
{
  triclinic_support = 0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldGPU::~EwaldGPU()
{
  ewald_gpu_clear(cpu_time);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void EwaldGPU::init()
{
  // unsupported configurations: the device path assumes a single standard
  // Verlet partition and a fixed domain decomposition

  if (strcmp(update->integrate_style, "verlet/split") == 0)
    error->all(FLERR, "Cannot use ewald/gpu with run_style verlet/split");

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "balance") == 0)
      error->all(FLERR, "Cannot currently use ewald/gpu with fix balance");

  // initialize the GPU device and atom storage first, so that the device is
  // ready when Ewald::init() -> setup() uploads the k-space coefficients

  int success = 0;
  ewald_gpu_init(atom->nlocal, atom->nlocal+atom->nghost, screen, success);
  GPU_EXTRA::check_flag(success, error, world);

  Ewald::init();
}

/* ----------------------------------------------------------------------
   called whenever the box changes: refresh the k-space coefficients on
   the device after the base class recomputes them
------------------------------------------------------------------------- */

void EwaldGPU::setup()
{
  Ewald::setup();

  int success = 0;
  ewald_gpu_setup(kmax, kcount, kxvecs, kyvecs, kzvecs, ug, eg, vg, unitk,
                       success);
  GPU_EXTRA::check_flag(success, error, world);
}

/* ----------------------------------------------------------------------
   compute the EwaldGPU long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }
  if (qsqsum == 0.0) return;

  const int nlocal = atom->nlocal;
  const int nall = atom->nlocal + atom->nghost;
  double *q = atom->q;

  // structure factors on the device

  bool success = true;
  ewald_gpu_structure(neighbor->ago, nlocal, nall, atom->x, atom->type,
                           atom->q, sfacrl, sfacim, success);
  if (!success)
    error->one(FLERR, "Insufficient memory on accelerator for ewald/gpu");

  const double t_cpu = MPI_Wtime();

  MPI_Allreduce(sfacrl, sfacrl_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim, sfacim_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  const double qscale = qqrd2e * scale;

  // per-atom field/force on the device (force queued for fix gpu to merge
  // with the pair force; raw per-atom energy/virial copied into eatom/vatom)

  ewald_gpu_compute(sfacrl_all, sfacim_all, qscale, slabflag, eflag_atom,
                         vflag_atom, eatom, vatom, success);
  if (!success)
    error->one(FLERR, "Insufficient memory on accelerator for ewald/gpu");

  // global energy across k-vectors plus self-energy and volume corrections

  if (eflag_global) {
    for (int k = 0; k < kcount; k++)
      energy += ug[k] * (sfacrl_all[k]*sfacrl_all[k] +
                         sfacim_all[k]*sfacim_all[k]);
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // global virial

  if (vflag_global) {
    double uk;
    for (int k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
      for (int j = 0; j < 6; j++) virial[j] += uk*vg[k][j];
    }
    for (int j = 0; j < 6; j++) virial[j] *= qscale;
  }

  // per-atom energy/virial: the device returned the raw k-vector sums; apply
  // the q_i factor, self-energy correction, and qscale here to match the CPU

  if (evflag_atom) {
    if (eflag_atom) {
      for (int i = 0; i < nlocal; i++) {
        eatom[i] = q[i]*eatom[i] - (g_ewald*q[i]*q[i]/MY_PIS +
          MY_PI2*q[i]*qsum / (g_ewald*g_ewald*volume));
        eatom[i] *= qscale;
      }
    }
    if (vflag_atom)
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < 6; j++) vatom[i][j] *= q[i]*qscale;
  }

  cpu_time += MPI_Wtime() - t_cpu;

  // 2d slab correction (host; adds to atom->f and the energy)

  if (slabflag == 1) slabcorr();
}

/* ---------------------------------------------------------------------- */

double EwaldGPU::memory_usage()
{
  double bytes = Ewald::memory_usage();
  bytes += ewald_gpu_bytes();
  return bytes;
}
