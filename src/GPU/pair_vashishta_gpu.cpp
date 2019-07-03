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
   Contributing author: Anders Hafreager (UiO)
------------------------------------------------------------------------- */

#include "pair_vashishta_gpu.h"
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

int vashishta_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                const double cell_size, int &gpu_mode, FILE *screen,
                int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                const double* cutsq, const double* r0,
                const double* gamma, const double* eta,
                const double* lam1inv, const double* lam4inv,
                const double* zizj, const double* mbigd,
                const double* dvrc, const double* big6w,
                const double* heta, const double* bigh,
                const double* bigw, const double* c0,
                const double* costheta, const double* bigb,
                const double* big2b, const double* bigc);
void vashishta_gpu_clear();
int ** vashishta_gpu_compute_n(const int ago, const int inum,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum,
                        const double cpu_time, bool &success);
void vashishta_gpu_compute(const int ago, const int nloc, const int nall, const int ln,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, int &host_start,
                    const double cpu_time, bool &success);
double vashishta_gpu_bytes();
extern double lmp_gpu_forces(double **f, double **tor, double *eatom,
                             double **vatom, double *virial, double &ecoul);

/* ---------------------------------------------------------------------- */

PairVashishtaGPU::PairVashishtaGPU(LAMMPS *lmp) : PairVashishta(lmp), gpu_mode(GPU_FORCE)
{
  cpu_time = 0.0;
  reinitflag = 0;
  gpu_allocated = false;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);

  cutghost = NULL;
  ghostneigh = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairVashishtaGPU::~PairVashishtaGPU()
{
  vashishta_gpu_clear();
  if (allocated)
    memory->destroy(cutghost);
}

/* ---------------------------------------------------------------------- */

void PairVashishtaGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = vashishta_gpu_compute_n(neighbor->ago, inum, nall,
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

    vashishta_gpu_compute(neighbor->ago, inum, nall, inum+list->gnum,
                   atom->x, atom->type, ilist, numneigh, firstneigh, eflag,
                   vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                   success);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");
}

/* ---------------------------------------------------------------------- */

void PairVashishtaGPU::allocate()
{
  if(!allocated) {
    PairVashishta::allocate();
  }
  int n = atom->ntypes;

  memory->create(cutghost,n+1,n+1,"pair:cutghost");
  gpu_allocated = true;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairVashishtaGPU::init_style()
{
  double cell_size = cutmax + neighbor->skin;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style vashishta/gpu requires atom IDs");
  if (force->newton_pair != 0)
    error->all(FLERR,"Pair style vashishta/gpu requires newton pair off");

  double *cutsq, *r0, *gamma, *eta;
  double *lam1inv, *lam4inv, *zizj, *mbigd;
  double *dvrc, *big6w, *heta, *bigh;
  double *bigw, *c0, *costheta, *bigb;
  double *big2b, *bigc;

  cutsq = r0 = gamma = eta = NULL;
  lam1inv = lam4inv = zizj = mbigd = NULL;
  dvrc = big6w = heta = bigh = NULL;
  bigw = c0 = costheta = bigb = NULL;
  big2b = bigc = NULL;

  memory->create(cutsq,nparams,"pair:cutsq");
  memory->create(r0,nparams,"pair:r0");
  memory->create(gamma,nparams,"pair:gamma");
  memory->create(eta,nparams,"pair:eta");
  memory->create(lam1inv,nparams,"pair:lam1inv");
  memory->create(lam4inv,nparams,"pair:lam4inv");
  memory->create(zizj,nparams,"pair:zizj");
  memory->create(mbigd,nparams,"pair:mbigd");
  memory->create(dvrc,nparams,"pair:dvrc");
  memory->create(big6w,nparams,"pair:big6w");
  memory->create(heta,nparams,"pair:heta");
  memory->create(bigh,nparams,"pair:bigh");
  memory->create(bigw,nparams,"pair:bigw");
  memory->create(c0,nparams,"pair:c0");
  memory->create(costheta,nparams,"pair:costheta");
  memory->create(bigb,nparams,"pair:bigb");
  memory->create(big2b,nparams,"pair:big2b");
  memory->create(bigc,nparams,"pair:bigc");

  for (int i = 0; i < nparams; i++) {
    cutsq[i] = params[i].cutsq;
    r0[i] = params[i].r0;
    gamma[i] = params[i].gamma;
    eta[i] = params[i].eta;
    lam1inv[i] = params[i].lam1inv;
    lam4inv[i] = params[i].lam4inv;
    zizj[i] = params[i].zizj;
    mbigd[i] = params[i].mbigd;
    dvrc[i] = params[i].dvrc;
    big6w[i] = params[i].big6w;
    heta[i] = params[i].heta;
    bigh[i] = params[i].bigh;
    bigw[i] = params[i].bigw;
    c0[i] = params[i].c0;
    costheta[i] = params[i].costheta;
    bigb[i] = params[i].bigb;
    big2b[i] = params[i].big2b;
    bigc[i] = params[i].bigc;
  }
  int success = vashishta_gpu_init(atom->ntypes+1, atom->nlocal, atom->nlocal+atom->nghost, 500,
                            cell_size, gpu_mode, screen, map, nelements,
                            elem2param, nparams, cutsq, r0, gamma, eta, lam1inv,
                            lam4inv, zizj, mbigd, dvrc, big6w, heta, bigh, bigw,
                            c0, costheta, bigb, big2b, bigc);
  memory->destroy(cutsq);
  memory->destroy(r0);
  memory->destroy(gamma);
  memory->destroy(eta);
  memory->destroy(lam1inv);
  memory->destroy(lam4inv);
  memory->destroy(zizj);
  memory->destroy(mbigd);
  memory->destroy(dvrc);
  memory->destroy(big6w);
  memory->destroy(heta);
  memory->destroy(bigh);
  memory->destroy(bigw);
  memory->destroy(c0);
  memory->destroy(costheta);
  memory->destroy(bigb);
  memory->destroy(big2b);
  memory->destroy(bigc);

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

double PairVashishtaGPU::init_one(int i, int j)
{
  if(!gpu_allocated) {
    allocate();
  }
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  cutghost[i][j] = cutmax;
  cutghost[j][i] = cutmax;

  return cutmax;
}

