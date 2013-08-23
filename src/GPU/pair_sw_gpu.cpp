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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
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

// External functions from cuda library for atom decomposition

int sw_gpu_init(const int nlocal, const int nall, const int max_nbors, 
                const double cell_size, int &gpu_mode, FILE *screen,
                const double, const double, const double, const double, 
                const double, const double, const double, const double, 
                const double, const double, const double);
void sw_gpu_clear();
int ** sw_gpu_compute_n(const int ago, const int inum,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, int *tag, int **nspecial,
                        int **special, const bool eflag, const bool vflag,
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

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWGPU::PairSWGPU(LAMMPS *lmp) : PairSW(lmp), gpu_mode(GPU_FORCE)
{
  cpu_time = 0.0;
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

    sw_gpu_compute(neighbor->ago, atom->nlocal, nall, inum+list->gnum,
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
  double cell_size = sqrt(params[0].cutsq) + neighbor->skin;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style sw/gpu requires atom IDs");
  if (force->newton_pair != 0)
    error->all(FLERR,"Pair style sw/gpu requires newton pair off");
  if (nparams > 1)
    error->all(FLERR,"Pair style sw/gpu is currently limited to one element.");

  int success = sw_gpu_init(atom->nlocal, atom->nlocal+atom->nghost, 300, 
                            cell_size, gpu_mode, screen,params[0].epsilon, 
                            params[0].sigma, params[0].lambda, params[0].gamma, 
                            params[0].costheta, params[0].biga, params[0].bigb, 
                            params[0].powerp, params[0].powerq, params[0].cut,
                            params[0].cutsq);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
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

