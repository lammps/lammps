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
   Contributing author: Vsevolod Nikolskiy (HSE)
------------------------------------------------------------------------- */

#include "pair_lj_cut_tip4p_long_gpu.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "kspace.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int ljtip4p_long_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                          const int nlocal, const int tH, const int tO, const double alpha,
                          const double qdist, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, const double host_cut_coulsq,
                          const double host_cut_coulsqplus, double *host_special_coul,
                          const double qqrd2e, const double g_ewald, int map_size, int max_same);
void ljtip4p_long_gpu_clear();
int **ljtip4p_long_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int *map_array, int map_size, int *sametag, int max_same,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom,
                                 int &host_start, int **ilist, int **jnum, const double cpu_time,
                                 bool &success, double *host_q, double *boxlo, double *prd);
void ljtip4p_long_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, int &host_start, const double cpu_time,
                              bool &success, double *host_q, const int nlocal, double *boxlo,
                              double *prd);
double ljtip4p_long_gpu_bytes();
void ljtip4p_long_copy_molecule_data(int, tagint *, int *, int, int *, int, int);

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongGPU::PairLJCutTIP4PLongGPU(LAMMPS *lmp) :
    PairLJCutTIP4PLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
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
    firstneigh = ljtip4p_long_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag,
        atom->get_map_array(), atom->get_map_size(), atom->sametag, atom->get_max_same(),
        atom->nspecial, atom->special, eflag, vflag, eflag_atom, vflag_atom, host_start, &ilist,
        &numneigh, cpu_time, success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    ljtip4p_long_copy_molecule_data(nall, atom->tag, atom->get_map_array(), atom->get_map_size(),
                                    atom->sametag, atom->get_max_same(), neighbor->ago);
    ljtip4p_long_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh,
                             firstneigh, eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                             success, atom->q, atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");
  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongGPU::init_style()
{

  cut_respa = nullptr;
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style lj/cut/tip4p/long/gpu requires atom IDs");
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/cut/tip4p/long/gpu requires atom attribute q");
  if (force->bond == nullptr) error->all(FLERR, "Must use a bond style with TIP4P potential");
  if (force->angle == nullptr) error->all(FLERR, "Must use an angle style with TIP4P potential");

  if (atom->map_style == Atom::MAP_HASH)
    error->all(FLERR,
               "GPU-accelerated pair style lj/cut/tip4p/long currently"
               " requires an 'array' style atom map (atom_modify map array)");

  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i, j);
        cut *= cut;
        if (cut > maxcut) maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  // ensure use of KSpace long-range solver, set g_ewald
  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables
  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;

  // set alpha parameter
  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5 * theta) * blen);

  cut_coulsq = cut_coul * cut_coul;
  double cut_coulplus = cut_coul + qdist + blen;
  double cut_coulsqplus = cut_coulplus * cut_coulplus;
  if (maxcut < cut_coulsqplus) cell_size = cut_coulplus + neighbor->skin;
  if (comm->get_comm_cutoff() < cell_size) {
    if (comm->me == 0)
      error->warning(FLERR, "Increasing communication cutoff to {:.8} for TIP4P GPU style",
                     cell_size);
    comm->cutghostuser = cell_size;
  }

  int mnf = 5e-2 * neighbor->oneatom;
  int success = ljtip4p_long_gpu_init(
      atom->ntypes + 1, cutsq, lj1, lj2, lj3, lj4, offset, force->special_lj, atom->nlocal, typeH,
      typeO, alpha, qdist, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode,
      screen, cut_ljsq, cut_coulsq, cut_coulsqplus, force->special_coul, force->qqrd2e, g_ewald,
      atom->get_map_size(), atom->get_max_same());
  GPU_EXTRA::check_flag(success, error, world);
  if (gpu_mode == GPU_FORCE) {
    auto req = neighbor->add_request(this, NeighConst::REQ_FULL);
    req->set_cutoff(cut_coulplus + neighbor->skin);
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutTIP4PLongGPU::memory_usage()
{
  double bytes = PairLJCutTIP4PLong::memory_usage();
  return bytes + ljtip4p_long_gpu_bytes();
}

/* ---------------------------------------------------------------------- */
