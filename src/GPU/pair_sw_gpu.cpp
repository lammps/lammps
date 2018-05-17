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
   Contributing author: Mike Brown (ORNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_sw_gpu.h"
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

int sw_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                const double cell_size, int &gpu_mode, FILE *screen,
                int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                const double* sw_epsilon, const double* sw_sigma,
                const double* sw_lambda, const double* sw_gamma,
                const double* sw_costheta, const double* sw_biga,
                const double* sw_bigb, const double* sw_powerp,
                const double* sw_powerq, const double* sw_cut,
                const double* sw_cutsq);
void sw_gpu_clear();
int ** sw_gpu_compute_n(const int ago, const int inum,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum,
                        const double cpu_time, bool &success);
void sw_gpu_compute(const int ago, const int nloc, const int nall, const int ln,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, int &host_start,
                    const double cpu_time, bool &success);
double sw_gpu_bytes();
extern double lmp_gpu_forces(double **f, double **tor, double *eatom,
                             double **vatom, double *virial, double &ecoul);

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWGPU::PairSWGPU(LAMMPS *lmp) : PairSW(lmp), gpu_mode(GPU_FORCE)
{
  cpu_time = 0.0;
  reinitflag = 0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);

  cutghost = NULL;
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWGPU::~PairSWGPU()
{
  sw_gpu_clear();
  if (allocated)
    memory->destroy(cutghost);
}

/* ---------------------------------------------------------------------- */

void PairSWGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = sw_gpu_compute_n(neighbor->ago, inum, nall,
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

    sw_gpu_compute(neighbor->ago, inum, nall, inum+list->gnum,
                   atom->x, atom->type, ilist, numneigh, firstneigh, eflag,
                   vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                   success);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");
}

/* ---------------------------------------------------------------------- */

void PairSWGPU::allocate()
{
  PairSW::allocate();
  int n = atom->ntypes;

  memory->create(cutghost,n+1,n+1,"pair:cutghost");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWGPU::init_style()
{
  double cell_size = cutmax + neighbor->skin;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style sw/gpu requires atom IDs");
  if (force->newton_pair != 0)
    error->all(FLERR,"Pair style sw/gpu requires newton pair off");

  double *epsilon, *sigma, *lambda, *gamma;
  double *biga, *bigb, *powerp, *powerq;
  double *_cut, *_cutsq, *costheta;
  epsilon = sigma = lambda = gamma = NULL;
  biga = bigb = powerp = powerq = NULL;
  _cut = _cutsq = costheta = NULL;

  memory->create(epsilon,nparams,"pair:epsilon");
  memory->create(sigma,nparams,"pair:sigma");
  memory->create(lambda,nparams,"pair:lambda");
  memory->create(gamma,nparams,"pair:gamma");
  memory->create(biga,nparams,"pair:biga");
  memory->create(bigb,nparams,"pair:bigb");
  memory->create(powerp,nparams,"pair:powerp");
  memory->create(powerq,nparams,"pair:powerq");
  memory->create(_cut,nparams,"pair:_cut");
  memory->create(_cutsq,nparams,"pair:_cutsq");
  memory->create(costheta,nparams,"pair:costheta");

  for (int i = 0; i < nparams; i++) {
    epsilon[i] = params[i].epsilon;
    sigma[i] = params[i].sigma;
    lambda[i] = params[i].lambda;
    gamma[i] = params[i].gamma;
    biga[i] = params[i].biga;
    bigb[i] = params[i].bigb;
    powerp[i] = params[i].powerp;
    powerq[i] = params[i].powerq;
    _cut[i] = params[i].cut;
    _cutsq[i] = params[i].cutsq;
    costheta[i] = params[i].costheta;
  }

  int success = sw_gpu_init(atom->ntypes+1, atom->nlocal, atom->nlocal+atom->nghost, 300,
                            cell_size, gpu_mode, screen, map, nelements,
                            elem2param, nparams, epsilon,
                            sigma, lambda, gamma, costheta, biga, bigb,
                            powerp, powerq, _cut, _cutsq);

  memory->destroy(epsilon);
  memory->destroy(sigma);
  memory->destroy(lambda);
  memory->destroy(gamma);
  memory->destroy(biga);
  memory->destroy(bigb);
  memory->destroy(powerp);
  memory->destroy(powerq);
  memory->destroy(_cut);
  memory->destroy(_cutsq);
  memory->destroy(costheta);

  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->ghost = 1;
  }

  if (comm->cutghostuser < (2.0*cutmax + neighbor->skin) )
    comm->cutghostuser=2.0*cutmax + neighbor->skin;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSWGPU::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  cutghost[i][j] = cutmax;
  cutghost[j][i] = cutmax;

  return cutmax;
}

