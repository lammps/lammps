// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

   Per-atom energy/virial added by Ray Shan (Materials Design, Inc.)
   Fix reaxff/bonds and fix reaxff/species for pair_style reaxff added
   by Ray Shan (Materials Design)

   OpenMP based threading support for pair_style reaxff/omp added
   by Hasan Metin Aktulga (MSU), Chris Knight (ALCF), Paul Coffman (ALCF),
   Kurt O'Hearn (MSU), Ray Shan (Materials Design), Wei Jiang (ALCF)

   Integration of the pair_style reaxff/omp into the OPENMP package
   by Axel Kohlmeyer (Temple U.)

   Please cite the related publication:
   H. M. Aktulga, C. Knight, P. Coffman, K. A. O'Hearn, T. R. Shan,
   W. Jiang, "Optimizing the performance of reactive molecular dynamics
   simulations for multi-core architectures", International Journal of
   High Performance Computing Applications, to appear.
 ------------------------------------------------------------------------- */

#include "pair_reaxff_omp.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "fix_reaxff.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

#include "reaxff_api.h"
#include "reaxff_omp.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "omp_compat.h"
#include "suffix.h"
using namespace LAMMPS_NS;
using namespace ReaxFF;

static const char cite_pair_reaxff_omp[] =
  "pair reaxff/omp and fix qeq/reaxff/omp command: doi:10.1177/1094342017746221\n\n"
  "@Article{Aktulga17,\n"
  " author =  {H. M. Aktulga and C. Knight and P. Coffman and\n"
  "    K. A. O'Hearn and T. R. Shan and W. Jiang},\n"
  " title =   {Optimizing the Performance of Reactive Molecular Dynamics\n"
  "    Simulations for Multi-Core Architectures},\n"
  " journal = {International Journal of High Performance Computing Applications},\n"
  " year =    2018\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairReaxFFOMP::PairReaxFFOMP(LAMMPS *lmp) : PairReaxFF(lmp), ThrOMP(lmp, THR_PAIR)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_reaxff_omp);

  suffix_flag |= Suffix::OMP;
  api->system->pair_ptr = this;
  api->system->omp_active = 1;

  num_nbrs_offset = nullptr;
}

/* ---------------------------------------------------------------------- */

PairReaxFFOMP::~PairReaxFFOMP()
{
  if (setup_flag) {
    reax_list * bonds = api->lists+BONDS;
    for (int i=0; i<bonds->num_intrs; ++i)
      sfree(error, bonds->select.bond_list[i].bo_data.CdboReduction, "CdboReduction");
  }
  memory->destroy(num_nbrs_offset);
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::init_style()
{
  if (!atom->q_flag) error->all(FLERR,"Pair style reaxff/omp requires atom attribute q");

  auto acks2_fixes = modify->get_fix_by_style("^acks2/reax");
  int have_qeq = modify->get_fix_by_style("^qeq/reax").size()
    + modify->get_fix_by_style("^qeq/shielded").size() + acks2_fixes.size();

  if (qeqflag && (have_qeq != 1))
    error->all(FLERR,"Pair style reaxff/omp requires use of exactly one of the "
               "fix qeq/reaxff or fix qeq/shielded or fix acks2/reaxff commands");

  api->system->acks2_flag = acks2_fixes.size();
  if (api->system->acks2_flag)
    error->all(FLERR,"Cannot (yet) use ACKS2 with OPENMP ReaxFF");

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts
  api->system->wsize = comm->nprocs;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style reaxff/omp requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style reaxff/omp requires newton pair on");

  // need a half neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  neighbor->add_request(this, NeighConst::REQ_GHOST | NeighConst::REQ_NEWTON_OFF);

  cutmax = MAX3(api->control->nonb_cut, api->control->hbond_cut, api->control->bond_cut);
  if ((cutmax < 2.0*api->control->bond_cut) && (comm->me == 0))
    error->warning(FLERR,"Total cutoff < 2*bond cutoff. May need to use an "
                   "increased neighbor list skin.");

  if (fix_reaxff == nullptr)
    fix_reaxff = dynamic_cast<FixReaxFF *>(modify->add_fix(fmt::format("{} all REAXFF",fix_id)));

  api->control->nthreads = comm->nthreads;
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::setup()
{
  int oldN;
  int mincap = api->system->mincap;
  double safezone = api->system->safezone;

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts
  oldN = api->system->N;

  if (api->system->N > nmax) {
    memory->destroy(num_nbrs_offset);
    // Don't update nmax here. It is updated at end of compute().
    memory->create(num_nbrs_offset, api->system->N, "pair:num_nbrs_offset");
  }

  if (setup_flag == 0) {

    setup_flag = 1;

    int *num_bonds = fix_reaxff->num_bonds;
    int *num_hbonds = fix_reaxff->num_hbonds;

    // determine the local and total capacity

    api->system->local_cap = MAX((int)(api->system->n * safezone), mincap);
    api->system->total_cap = MAX((int)(api->system->N * safezone), mincap);

    // initialize my data structures

    PreAllocate_Space(api->system, api->workspace);
    write_reax_atoms();

    api->system->wsize = comm->nprocs;

    int num_nbrs = estimate_reax_lists();
    if (num_nbrs < 0)
      error->all(FLERR,"Too many neighbors for pair style reaxff");

    Make_List(api->system->total_cap,num_nbrs,TYP_FAR_NEIGHBOR,api->lists+FAR_NBRS);
    (api->lists+FAR_NBRS)->error_ptr=error;

    write_reax_lists();

    InitializeOMP(api->system,api->control,api->data,api->workspace,&api->lists,world);
    for (int k = 0; k < api->system->N; ++k) {
      num_bonds[k] = api->system->my_atoms[k].num_bonds;
      num_hbonds[k] = api->system->my_atoms[k].num_hbonds;
    }

  } else {

    // fill in reax datastructures

    write_reax_atoms();

    // reset the bond list info for new atoms

    for (int k = oldN; k < api->system->N; ++k)
      Set_End_Index(k, Start_Index(k, api->lists+BONDS), api->lists+BONDS);

    // estimate far neighbor list size
    // Not present in MPI-only version
    api->workspace->realloc.num_far = estimate_reax_lists();

    // check if I need to shrink/extend my data-structs

    ReAllocate(api->system, api->control, api->data, api->workspace, &api->lists);
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::compute(int eflag, int vflag)
{

  // communicate num_bonds once every reneighboring
  // 2 num arrays stored by fix, grab ptr to them

  if (neighbor->ago == 0) comm->forward_comm(fix_reaxff);
  int *num_bonds = fix_reaxff->num_bonds;
  int *num_hbonds = fix_reaxff->num_hbonds;

  ev_init(eflag,vflag);

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts
  const int nall = api->system->N;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
#if defined(_OPENMP)
     int tid = omp_get_thread_num();
#else
     int tid = 0;
#endif
     ThrData *thr = fix->get_thr(tid);
     thr->timer(Timer::START);
     ev_setup_thr(eflag, vflag, api->system->N, eatom, vatom, nullptr, thr);
  }
  // setup data structures

  setup();

  Reset(api->system, api->control, api->data, api->workspace, &api->lists);

  // Why not update workspace like in MPI-only code?
  // Using the MPI-only way messes up the hb energy
  //workspace->realloc.num_far = write_reax_lists();
  write_reax_lists();

  // forces

  Compute_ForcesOMP(api->system,api->control,api->data,api->workspace,&api->lists);
  read_reax_forces(vflag);

  const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(vflag)
#endif
  {
#if defined(_OPENMP)
     int tid = omp_get_thread_num();
#else
     int tid = 0;
#endif
     ThrData *thr = fix->get_thr(tid);
     thr->timer(Timer::PAIR);

     // the pair style reduces energy and forces directly. so only reduce virial/
     // per-atom virial and per-atom centroid virial are the same for two-body
     // many-body pair styles not yet implemented
      if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
        data_reduce_thr(&(vatom[0][0]), nall , nthreads, 6, tid);
      }
  }

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < api->system->N; ++k) {
    num_bonds[k] = api->system->my_atoms[k].num_bonds;
    num_hbonds[k] = api->system->my_atoms[k].num_hbonds;
  }

  // energies and pressure

  if (eflag_global) {

    // Store the different parts of the energy
    // in a list for output by compute pair command

    pvector[0] = api->data->my_en.e_bond;
    pvector[1] = api->data->my_en.e_ov + api->data->my_en.e_un;
    pvector[2] = api->data->my_en.e_lp;
    pvector[3] = 0.0;
    pvector[4] = api->data->my_en.e_ang;
    pvector[5] = api->data->my_en.e_pen;
    pvector[6] = api->data->my_en.e_coa;
    pvector[7] = api->data->my_en.e_hb;
    pvector[8] = api->data->my_en.e_tor;
    pvector[9] = api->data->my_en.e_con;
    pvector[10] = api->data->my_en.e_vdW;
    pvector[11] = api->data->my_en.e_ele;
    pvector[12] = 0.0;
    pvector[13] = api->data->my_en.e_pol;
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // Set internal timestep counter to that of LAMMPS

  api->data->step = update->ntimestep;

  // populate tmpid and tmpbo arrays for fix reaxff/species

  if (fixspecies_flag) {
    if (api->system->N > nmax) {
      memory->destroy(tmpid);
      memory->destroy(tmpbo);
      nmax = api->system->N;
      memory->create(tmpid,nmax,MAXSPECBOND,"pair:tmpid");
      memory->create(tmpbo,nmax,MAXSPECBOND,"pair:tmpbo");
    }

#if defined(_OPENMP)
#pragma omp parallel for collapse(2) schedule(static) default(shared)
#endif
    for (int i = 0; i < api->system->N; i++)
      for (int j = 0; j < MAXSPECBOND; j++) {
        tmpbo[i][j] = 0.0;
        tmpid[i][j] = 0;
      }

    FindBond();
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::write_reax_atoms()
{
  int *num_bonds = fix_reaxff->num_bonds;
  int *num_hbonds = fix_reaxff->num_hbonds;

  if (api->system->N > api->system->total_cap)
    error->all(FLERR,"Too many ghost atoms");

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (int i = 0; i < api->system->N; ++i) {
    api->system->my_atoms[i].orig_id = atom->tag[i];
    api->system->my_atoms[i].type = map[atom->type[i]];
    api->system->my_atoms[i].x[0] = atom->x[i][0];
    api->system->my_atoms[i].x[1] = atom->x[i][1];
    api->system->my_atoms[i].x[2] = atom->x[i][2];
    api->system->my_atoms[i].q = atom->q[i];
    api->system->my_atoms[i].num_bonds = num_bonds[i];
    api->system->my_atoms[i].num_hbonds = num_hbonds[i];
  }
}

/* ---------------------------------------------------------------------- */

int PairReaxFFOMP::estimate_reax_lists()
{
  int i;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int numall = list->inum + list->gnum;
  int mincap = api->system->mincap;

  // for good performance in the OpenMP implementation, each thread needs
  // to know where to place the neighbors of the atoms it is responsible for.
  // The sumscan values for the list->numneigh will be used to determine the
  // neighbor offset of each atom. Note that this may cause some significant
  // memory overhead if delayed neighboring is used - so it may be desirable
  // to work on this part to reduce the memory footprint of the far_nbrs list.

  int num_nbrs = 0;

  for (int itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    num_nbrs += numneigh[i];
  }

  int new_estimate = MAX(num_nbrs, mincap*REAX_MIN_NBRS);

  return new_estimate;
}

/* ---------------------------------------------------------------------- */

int PairReaxFFOMP::write_reax_lists()
{
  int num_mynbrs;
  double d_sqr, dist, cutoff_sqr;
  rvec dvec;

  double **x = atom->x;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  reax_list *far_nbrs = api->lists + FAR_NBRS;
  far_neighbor_data *far_list = far_nbrs->select.far_nbr_list;

  int num_nbrs = 0;
  int inum = list->inum;
  int gnum = list->gnum;
  int numall = inum + gnum;

  // sumscan of the number of neighbors per atom to determine the offsets
  // most likely, we are overallocating. desirable to work on this part
  // to reduce the memory footprint of the far_nbrs list.

  num_nbrs = 0;

  for (int itr_i = 0; itr_i < numall; ++itr_i) {
    int i = ilist[itr_i];
    num_nbrs_offset[i] = num_nbrs;
    num_nbrs += numneigh[i];
  }

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) default(shared)           \
  private(cutoff_sqr, num_mynbrs, d_sqr, dvec, dist)
#endif
  for (int itr_i = 0; itr_i < numall; ++itr_i) {
    int i = ilist[itr_i];
    auto jlist = firstneigh[i];
    Set_Start_Index(i, num_nbrs_offset[i], far_nbrs);

    if (i < inum)
      cutoff_sqr = SQR(api->control->nonb_cut);
    else
      cutoff_sqr = SQR(api->control->bond_cut);

    num_mynbrs = 0;

    for (int itr_j = 0; itr_j < numneigh[i]; ++itr_j) {
      int j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance(x[j], x[i], &d_sqr, &dvec);

      if (d_sqr <= cutoff_sqr) {
        dist = sqrt(d_sqr);
        set_far_nbr(&far_list[num_nbrs_offset[i] + num_mynbrs], j, dist, dvec);
        ++num_mynbrs;
      }
    }
    Set_End_Index(i, num_nbrs_offset[i] + num_mynbrs, far_nbrs);
  }

  return num_nbrs;
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::read_reax_forces(int /* vflag */)
{
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (int i = 0; i < api->system->N; ++i) {
    api->system->my_atoms[i].f[0] = api->workspace->f[i][0];
    api->system->my_atoms[i].f[1] = api->workspace->f[i][1];
    api->system->my_atoms[i].f[2] = api->workspace->f[i][2];

    atom->f[i][0] = -api->workspace->f[i][0];
    atom->f[i][1] = -api->workspace->f[i][1];
    atom->f[i][2] = -api->workspace->f[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFFOMP::FindBond()
{
  const double bo_cut = 0.10;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (int i = 0; i < api->system->n; i++) {
    int j, pj, nj;
    double bo_tmp;
    bond_data *bo_ij;

    nj = 0;
    for (pj = Start_Index(i, api->lists); pj < End_Index(i, api->lists); ++pj) {
      bo_ij = &(api->lists->select.bond_list[pj]);
      j = bo_ij->nbr;
      if (j < i) continue;

      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp >= bo_cut) {
        tmpid[i][nj] = j;
        tmpbo[i][nj] = bo_tmp;
        nj ++;
        if (nj > MAXSPECBOND) error->all(FLERR,"Increase MAXSPECBOND in fix_reaxff_species.h");
      }
    }
  }
}
