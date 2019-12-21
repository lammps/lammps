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
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#include "pair_tersoff_gpu.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int tersoff_gpu_init(const int ntypes, const int inum, const int nall,
                     const int max_nbors, const double cell_size, int &gpu_mode,
                     FILE *screen, int* host_map, const int nelements,
                     int*** host_elem2param, const int nparams,
                     const double* ts_lam1, const double* ts_lam2,
                     const double* ts_lam3, const double* ts_powermint,
                     const double* ts_biga, const double* ts_bigb,
                     const double* ts_bigr, const double* ts_bigd,
                     const double* ts_c1, const double* ts_c2,
                     const double* ts_c3, const double* ts_c4,
                     const double* ts_c, const double* ts_d,
                     const double* ts_h, const double* ts_gamma,
                     const double* ts_beta, const double* ts_powern,
                     const double* ts_cutsq);
void tersoff_gpu_clear();
int ** tersoff_gpu_compute_n(const int ago, const int inum_full,
                    const int nall, double **host_x, int *host_type,
                    double *sublo, double *subhi, tagint *tag, int **nspecial,
                    tagint **special, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, int &host_start,
                    int **ilist, int **jnum, const double cpu_time,
                    bool &success);
void tersoff_gpu_compute(const int ago, const int nlocal, const int nall,
                    const int nlist, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom,
                    int &host_start, const double cpu_time, bool &success);
double tersoff_gpu_bytes();
extern double lmp_gpu_forces(double **f, double **tor, double *eatom,
                             double **vatom, double *virial, double &ecoul);

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairTersoffGPU::PairTersoffGPU(LAMMPS *lmp) : PairTersoff(lmp), gpu_mode(GPU_FORCE)
{
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);

  cutghost = NULL;
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoffGPU::~PairTersoffGPU()
{
  tersoff_gpu_clear();
  if (allocated)
    memory->destroy(cutghost);
}

/* ---------------------------------------------------------------------- */

void PairTersoffGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = tersoff_gpu_compute_n(neighbor->ago, inum, nall,
                                  atom->x, atom->type, domain->sublo,
                                  domain->subhi, atom->tag, atom->nspecial,
                                  atom->special, eflag, vflag, eflag_atom,
                                  vflag_atom, host_start,
                                  &ilist, &numneigh, cpu_time, success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    tersoff_gpu_compute(neighbor->ago, inum, nall, inum+list->gnum,
                   atom->x, atom->type, ilist, numneigh, firstneigh, eflag,
                   vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                   success);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");
}

/* ---------------------------------------------------------------------- */

void PairTersoffGPU::allocate()
{
  PairTersoff::allocate();
  int n = atom->ntypes;

  memory->create(cutghost,n+1,n+1,"pair:cutghost");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTersoffGPU::init_style()
{
  double cell_size = cutmax + neighbor->skin;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style tersoff/gpu requires atom IDs");
  if (force->newton_pair != 0)
    error->all(FLERR,"Pair style tersoff/gpu requires newton pair off");

  double *lam1, *lam2, *lam3, *powermint;
  double *biga, *bigb, *bigr, *bigd;
  double *c1, *c2, *c3, *c4;
  double *c, *d, *h, *gamma;
  double *beta, *powern, *_cutsq;
  lam1 = lam2 = lam3 = powermint = NULL;
  biga = bigb = bigr = bigd = NULL;
  c1 = c2 = c3 = c4 = NULL;
  c = d = h = gamma = NULL;
  beta = powern = _cutsq = NULL;

  memory->create(lam1,nparams,"pair:lam1");
  memory->create(lam2,nparams,"pair:lam2");
  memory->create(lam3,nparams,"pair:lam3");
  memory->create(powermint,nparams,"pair:powermint");
  memory->create(biga,nparams,"pair:biga");
  memory->create(bigb,nparams,"pair:bigb");
  memory->create(bigr,nparams,"pair:bigr");
  memory->create(bigd,nparams,"pair:bigd");
  memory->create(c1,nparams,"pair:c1");
  memory->create(c2,nparams,"pair:c2");
  memory->create(c3,nparams,"pair:c3");
  memory->create(c4,nparams,"pair:c4");
  memory->create(c,nparams,"pair:c");
  memory->create(d,nparams,"pair:d");
  memory->create(h,nparams,"pair:h");
  memory->create(gamma,nparams,"pair:gamma");
  memory->create(beta,nparams,"pair:beta");
  memory->create(powern,nparams,"pair:powern");
  memory->create(_cutsq,nparams,"pair:_cutsq");

  for (int i = 0; i < nparams; i++) {
    lam1[i] = params[i].lam1;
    lam2[i] = params[i].lam2;
    lam3[i] = params[i].lam3;
    powermint[i] = params[i].powermint;
    biga[i] = params[i].biga;
    bigb[i] = params[i].bigb;
    bigr[i] = params[i].bigr;
    bigd[i] = params[i].bigd;
    c1[i] = params[i].c1;
    c2[i] = params[i].c2;
    c3[i] = params[i].c3;
    c4[i] = params[i].c4;
    c[i] = params[i].c;
    d[i] = params[i].d;
    h[i] = params[i].h;
    gamma[i] = params[i].gamma;
    beta[i] = params[i].beta;
    powern[i] = params[i].powern;
    _cutsq[i] = params[i].cutsq;
  }

  int success = tersoff_gpu_init(atom->ntypes+1, atom->nlocal,
                                 atom->nlocal+atom->nghost, 300,
                                 cell_size, gpu_mode, screen, map, nelements,
                                 elem2param, nparams, lam1, lam2, lam3,
                                 powermint, biga, bigb, bigr, bigd,
                                 c1, c2, c3, c4, c, d, h, gamma,
                                 beta, powern, _cutsq);

  memory->destroy(lam1);
  memory->destroy(lam2);
  memory->destroy(lam3);
  memory->destroy(powermint);
  memory->destroy(biga);
  memory->destroy(bigb);
  memory->destroy(bigr);
  memory->destroy(bigd);
  memory->destroy(c1);
  memory->destroy(c2);
  memory->destroy(c3);
  memory->destroy(c4);
  memory->destroy(c);
  memory->destroy(d);
  memory->destroy(h);
  memory->destroy(gamma);
  memory->destroy(beta);
  memory->destroy(powern);
  memory->destroy(_cutsq);

  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->ghost = 1;
  }

  if (comm->cutghostuser < (2.0*cutmax + neighbor->skin)) {
    comm->cutghostuser = 2.0*cutmax + neighbor->skin;
    if (comm->me == 0)
       error->warning(FLERR,"Increasing communication cutoff for GPU style");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoffGPU::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  cutghost[i][j] = cutmax;
  cutghost[j][i] = cutmax;

  return cutmax;
}

