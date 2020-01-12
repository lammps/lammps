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
   Contributing author: Vsevolod Nikolskiy (HSE)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_cut_tip4p_long_gpu.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "kspace.h"
#include "angle.h"
#include "bond.h"
#include "gpu_extra.h"

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int ljtip4p_long_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
    double **host_lj2, double **host_lj3, double **host_lj4,
    double **offset, double *special_lj, const int nlocal,
    const int tH, const int tO, const double alpha, const double qdist,
    const int nall, const int max_nbors, const int maxspecial,
    const double cell_size, int &gpu_mode, FILE *screen,
    double **host_cut_ljsq, const double host_cut_coulsq,
    const double host_cut_coulsqplus, double *host_special_coul,
    const double qqrd2e, const double g_ewald,
    int map_size, int max_same);
void ljtip4p_long_gpu_clear();
int ** ljtip4p_long_gpu_compute_n(const int ago, const int inum,
    const int nall, double **host_x, int *host_type,
    double *sublo, double *subhi,
    tagint *tag, int *map_array, int map_size,
    int *sametag, int max_same,
    int **nspecial,
    tagint **special, const bool eflag, const bool vflag,
    const bool eatom, const bool vatom, int &host_start,
    int **ilist, int **jnum,
    const double cpu_time, bool &success, double *host_q,
    double *boxlo, double *prd);
void ljtip4p_long_gpu_compute(const int ago, const int inum, const int nall,
    double **host_x, int *host_type, int *ilist, int *numj,
    int **firstneigh, const bool eflag, const bool vflag,
    const bool eatom, const bool vatom, int &host_start,
    const double cpu_time,
    bool &success, double *host_q, const int nlocal,
    double *boxlo, double *prd);
double ljtip4p_long_gpu_bytes();
void ljtip4p_long_copy_molecule_data(int, tagint *, int *,
                                     int, int *, int, int);

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongGPU::PairLJCutTIP4PLongGPU(LAMMPS *lmp)
: PairLJCutTIP4PLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutTIP4PLongGPU::~PairLJCutTIP4PLongGPU()
{
  ljtip4p_long_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCutTIP4PLongGPU::compute(int eflag, int vflag)
{

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = ljtip4p_long_gpu_compute_n(neighbor->ago, inum, nall,
        atom->x, atom->type, domain->sublo,
        domain->subhi,
        atom->tag, atom->get_map_array(), atom->get_map_size(),
        atom->sametag, atom->get_max_same(),
        atom->nspecial,
        atom->special, eflag, vflag, eflag_atom,
        vflag_atom, host_start, &ilist, &numneigh,
        cpu_time, success, atom->q, domain->boxlo,
        domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    ljtip4p_long_copy_molecule_data(nall, atom->tag,
        atom->get_map_array(), atom->get_map_size(),
        atom->sametag, atom->get_max_same(), neighbor->ago);
    ljtip4p_long_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
        ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
        vflag_atom, host_start, cpu_time, success, atom->q,
        atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

//  if (host_start<inum) {
//    cpu_time = MPI_Wtime();
//    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
//    cpu_time = MPI_Wtime() - cpu_time;
//  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongGPU::init_style()
{

  cut_respa = NULL;
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style lj/cut/tip4p/long/gpu requires atom IDs");
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style lj/cut/tip4p/long/gpu requires atom attribute q");
  if (force->bond == NULL)
    error->all(FLERR,"Must use a bond style with TIP4P potential");
  if (force->angle == NULL)
    error->all(FLERR,"Must use an angle style with TIP4P potential");

  if (atom->map_style == 2)
    error->all(FLERR,"GPU-accelerated lj/cut/tip4p/long currently"
        " requires map style 'array' (atom_modify map array)");

  //PairLJCutCoulLong::init_style();
  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  // insure use of KSpace long-range solver, set g_ewald
  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables
  if (ncoultablebits) init_tables(cut_coul,cut_respa);

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;

  // set alpha parameter
  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);

  cut_coulsq = cut_coul * cut_coul;
  double cut_coulsqplus = (cut_coul+qdist+blen) * (cut_coul+qdist+blen);
  if (maxcut < cut_coulsqplus) {
    cell_size = (cut_coul+qdist+blen) + neighbor->skin;
  }
  if (comm->cutghostuser < cell_size) {
    comm->cutghostuser = cell_size;
    if (comm->me == 0)
      error->warning(FLERR,"Increasing communication cutoff for TIP4P GPU style");
  }

  int success = ljtip4p_long_gpu_init(atom->ntypes+1, cutsq, lj1, lj2, lj3, lj4,
                             offset, force->special_lj, atom->nlocal,
                             typeH, typeO, alpha, qdist,
                             atom->nlocal+atom->nghost, 300, maxspecial,
                             cell_size, gpu_mode, screen, cut_ljsq,
                             cut_coulsq, cut_coulsqplus,
                             force->special_coul, force->qqrd2e,
                             g_ewald, atom->get_map_size(),
                             atom->get_max_same());
  GPU_EXTRA::check_flag(success,error,world);
  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = cut_coul+qdist+blen + neighbor->skin;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutTIP4PLongGPU::memory_usage()
{
  double bytes = PairLJCutTIP4PLong::memory_usage();
  return bytes + ljtip4p_long_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

