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
   Contributing author: Mike Brown (ORNL)
------------------------------------------------------------------------- */

#include "pair_sw_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int sw_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                const double cell_size, int &gpu_mode, FILE *screen, double **ncutsq, double **ncut,
                double **sigma, double **powerp, double **powerq, double **sigma_gamma, double **c1,
                double **c2, double **c3, double **c4, double **c5, double **c6,
                double ***lambda_epsilon, double ***costheta, const int *map, int ***e2param);
void sw_gpu_clear();
int **sw_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int &host_start, int **ilist, int **jnum,
                       const double cpu_time, bool &success);
void sw_gpu_compute(const int ago, const int nloc, const int nall, const int ln, double **host_x,
                    int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom, int &host_start,
                    const double cpu_time, bool &success);
double sw_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairSWGPU::PairSWGPU(LAMMPS *lmp) : PairSW(lmp), gpu_mode(GPU_FORCE)
{
  cpu_time = 0.0;
  reinitflag = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);

  cutghost = nullptr;
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWGPU::~PairSWGPU()
{
  sw_gpu_clear();
  if (allocated) memory->destroy(cutghost);
}

/* ---------------------------------------------------------------------- */

void PairSWGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
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
    firstneigh =
        sw_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag,
                         atom->nspecial, atom->special, eflag, vflag, eflag_atom, vflag_atom,
                         host_start, &ilist, &numneigh, cpu_time, success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    sw_gpu_compute(neighbor->ago, inum, nall, inum + list->gnum, atom->x, atom->type, ilist,
                   numneigh, firstneigh, eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                   success);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");
  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ---------------------------------------------------------------------- */

void PairSWGPU::allocate()
{
  PairSW::allocate();
  int n = atom->ntypes;

  memory->create(cutghost, n + 1, n + 1, "pair:cutghost");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWGPU::init_style()
{
  double cell_size = cutmax + neighbor->skin;

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style sw/gpu requires atom IDs");

  double **c1, **c2, **c3, **c4, **c5, **c6;
  double **ncutsq, **ncut, **sigma, **powerp, **powerq, **sigma_gamma;
  double ***lambda_epsilon, ***costheta;
  c1 = c2 = c3 = c4 = c5 = c6 = nullptr;
  ncutsq = ncut = sigma = powerp = powerq = sigma_gamma = nullptr;
  lambda_epsilon = costheta = nullptr;

  const int tp1 = atom->ntypes + 1;

  memory->create(ncutsq, tp1, tp1, "pair:ncutsq");
  memory->create(ncut, tp1, tp1, "pair:ncut");
  memory->create(sigma, tp1, tp1, "pair:sigma");
  memory->create(powerp, tp1, tp1, "pair:powerp");
  memory->create(powerq, tp1, tp1, "pair:powerq");
  memory->create(sigma_gamma, tp1, tp1, "pair:sigma_gamma");
  memory->create(c1, tp1, tp1, "pair:c1");
  memory->create(c2, tp1, tp1, "pair:c2");
  memory->create(c3, tp1, tp1, "pair:c3");
  memory->create(c4, tp1, tp1, "pair:c4");
  memory->create(c5, tp1, tp1, "pair:c5");
  memory->create(c6, tp1, tp1, "pair:c6");
  memory->create(lambda_epsilon, tp1, tp1, tp1, "pair:lambda_epsilon");
  memory->create(costheta, tp1, tp1, tp1, "pair:costheta");

  for (int ii = 1; ii < tp1; ii++) {
    int i = map[ii];
    for (int jj = 1; jj < tp1; jj++) {
      int j = map[jj];
      if (i < 0 || j < 0)
        continue;
      else {
        int ijparam = elem3param[i][j][j];
        ncutsq[ii][jj] = params[ijparam].cutsq;
        ncut[ii][jj] = params[ijparam].cut;
        sigma[ii][jj] = params[ijparam].sigma;
        powerp[ii][jj] = params[ijparam].powerp;
        powerq[ii][jj] = params[ijparam].powerq;
        sigma_gamma[ii][jj] = params[ijparam].sigma_gamma;
        c1[ii][jj] = params[ijparam].c1;
        c2[ii][jj] = params[ijparam].c2;
        c3[ii][jj] = params[ijparam].c3;
        c4[ii][jj] = params[ijparam].c4;
        c5[ii][jj] = params[ijparam].c5;
        c6[ii][jj] = params[ijparam].c6;
      }

      for (int kk = 1; kk < tp1; kk++) {
        int k = map[kk];
        if (k < 0)
          continue;
        else {
          int ijkparam = elem3param[i][j][k];
          costheta[ii][jj][kk] = params[ijkparam].costheta;
          lambda_epsilon[ii][jj][kk] = params[ijkparam].lambda_epsilon;
        }
      }
    }
  }

  int mnf = 5e-2 * neighbor->oneatom;
  int success = sw_gpu_init(tp1, atom->nlocal, atom->nlocal + atom->nghost, mnf, cell_size,
                            gpu_mode, screen, ncutsq, ncut, sigma, powerp, powerq, sigma_gamma, c1,
                            c2, c3, c4, c5, c6, lambda_epsilon, costheta, map, elem3param);

  memory->destroy(ncutsq);
  memory->destroy(ncut);
  memory->destroy(sigma);
  memory->destroy(powerp);
  memory->destroy(powerq);
  memory->destroy(sigma_gamma);
  memory->destroy(c1);
  memory->destroy(c2);
  memory->destroy(c3);
  memory->destroy(c4);
  memory->destroy(c5);
  memory->destroy(c6);
  memory->destroy(lambda_epsilon);
  memory->destroy(costheta);

  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE)
    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
  if (comm->get_comm_cutoff() < (2.0 * cutmax + neighbor->skin)) {
    comm->cutghostuser = 2.0 * cutmax + neighbor->skin;
    if (comm->me == 0)
      error->warning(FLERR, "Increasing communication cutoff to {:.8} for GPU pair style",
                     comm->cutghostuser);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSWGPU::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  cutghost[i][j] = cutmax;
  cutghost[j][i] = cutmax;

  return cutmax;
}
