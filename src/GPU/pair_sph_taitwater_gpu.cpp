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
   Contributing author: Trung Dac Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "pair_sph_taitwater_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "info.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int sph_taitwater_gpu_init(const int ntypes, double **cutsq, double** host_cut,
                           double **host_viscosity, double* host_mass, double* host_rho0,
                           double* host_soundspeed, double* host_B, const int dimension,
                           double *special_lj, const int inum, const int nall,
                           const int max_nbors,  const int maxspecial,
                           const double cell_size, int &gpu_mode, FILE *screen);
void sph_taitwater_gpu_clear();
int **sph_taitwater_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time, bool &success,
                         double **host_v);
void sph_taitwater_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, tagint *tag,
                        double **host_v);
void sph_taitwater_gpu_get_extra_data(double *host_rho);
void sph_taitwater_gpu_update_drhoE(void **drhoE_ptr);
double sph_taitwater_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterGPU::PairSPHTaitwaterGPU(LAMMPS *lmp) : PairSPHTaitwater(lmp), gpu_mode(GPU_FORCE)
{
  drhoE_pinned = nullptr;
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairSPHTaitwaterGPU::~PairSPHTaitwaterGPU()
{
  sph_taitwater_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairSPHTaitwaterGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;

  double *rho = atom->rho;
  sph_taitwater_gpu_get_extra_data(rho);

  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = sph_taitwater_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
        cpu_time, success, atom->v);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    sph_taitwater_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                       eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time, success,
                       atom->tag, atom->v);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  // get the drho and dE from device

  double *drho = atom->drho;
  double *desph = atom->desph;
  sph_taitwater_gpu_update_drhoE(&drhoE_pinned);

  int nlocal = atom->nlocal;
  if (acc_float) {
    auto drhoE_ptr = (float *)drhoE_pinned;
    int idx = 0;
    for (int i = 0; i < nlocal; i++) {
      drho[i] = drhoE_ptr[idx];
      desph[i] = drhoE_ptr[idx+1];
      idx += 2;
    }

  } else {
    auto drhoE_ptr = (double *)drhoE_pinned;
    int idx = 0;
    for (int i = 0; i < nlocal; i++) {
      drho[i] = drhoE_ptr[idx];
      desph[i] = drhoE_ptr[idx+1];
      idx += 2;
    }
  }

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSPHTaitwaterGPU::init_style()
{

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double mcut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mcut = init_one(i, j);
        mcut *= mcut;
        if (mcut > maxcut) maxcut = mcut;
        cutsq[i][j] = cutsq[j][i] = mcut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      sph_taitwater_gpu_init(atom->ntypes + 1, cutsq, cut, viscosity, atom->mass,
                             rho0, soundspeed, B, domain->dimension, force->special_lj,
                             atom->nlocal, atom->nlocal + atom->nghost,
                             mnf, maxspecial, cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success, error, world);

  acc_float = Info::has_accelerator_feature("GPU", "precision", "single");

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairSPHTaitwaterGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + sph_taitwater_gpu_bytes();
}
