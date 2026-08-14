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
   Author: Brian Dandurand (Queen's University Belfast)
   Based on VerletSplit as originally developed by:
        Yuxing Peng and Chris Knight (U Chicago)
   This class VerletSplitRK corresponds to the enhanced baseline of:
        Brian Dandurand, Hans Vandierendonck, and Bronis de Supinski.
        "Improving Parallel Scalability for Molecular Dynamics Simulations in the Exascale Era".
        in Proceedings of the IPDPS Conference. 2025.
   The enhanced baseline in turn was inspired by the earlier contribution of
           D. F. Richards, J. N. Glosli, B. Chan, M. R. Dorr, E. W. Draeger, J.-
        L. Fattebert, W. D. Krauss, T. Spelce, F. H. Streitz, M. P. Surh, and
        J. A. Gunnels,
        "Beyond homogeneous decomposition: scaling long-range forces
        on massively parallel systems,"
        in Proceedings of the Conference on High Performance Computing Networking,
        Storage and Analysis, ser. SC '09. New York, NY, USA:
        Association for Computing Machinery, 2009.
------------------------------------------------------------------------- */

#include "verlet_split_rk.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

using namespace LAMMPS_NS;

static const char cite_verlet_split_rk[] =
    "run_style verlet/split/rk command: "
    "https://doi.org/10.1109/IPDPS64566.2025.00077\n\n"
    "@inproceedings{11078563,\n"
    "author={Dandurand, Brian and Vandierendonck, Hans and De Supinski, Bronis R.},\n"
    "booktitle={2025 IEEE International Parallel and Distributed Processing Symposium (IPDPS)},\n"
    "title={Improving Parallel Scalability for Molecular Dynamics Simulations in the Exascale Era},\n"
    "year={2025},\n"
    "volume={},\n"
    "number={},\n"
    "pages={813-823},\n"
    "keywords={Distributed processing;Accuracy;Scalability;Computational modeling;Parallel processing;"
    "Trajectory;Timing;Maintenance;History;Optimization;molecular dynamics simulation;parallelism;"
    "asynchrony;scalability;LAMMPS;MPI},\n"
    "doi={10.1109/IPDPS64566.2025.00077}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

VerletSplitRK::VerletSplitRK(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg)
{
  rk_flag = 1;
  if (lmp->citeme) {
      lmp->citeme->add(cite_verlet_split_rk);
  }
}

/* ---------------------------------------------------------------------- */

VerletSplitRK::~VerletSplitRK()
{
  ;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletSplitRK::init()
{
  if (comm->style != Comm::BRICK)
    error->universe_all(FLERR,"Verlet/split/rk can only currently be used with comm_style brick");
  if (!force->kspace)
    error->universe_all(FLERR,"A KSpace style must be defined with verlet/split/rk");

  // error for as-yet unsupported verlet/split KSpace options

  int errflag = 0;
  if (!force->kspace->rk_flag) //Currently requires the use of pppm/rk
    error->all(FLERR,"The kspace style {} does not support an rk decomposition. Verlet/split/rk requires a kspace style (such as pppm/rk) that supports an rk decomposition.", force->kspace_style);
  if (force->kspace->tip4pflag) errflag = 1;
  if (force->kspace->dipoleflag) errflag = 1;
  if (force->kspace->spinflag) errflag = 1;

  if (errflag)
    error->all(FLERR,"Verlet/split/rk cannot (yet) be used with kspace style {}", force->kspace_style);


  // invoke parent Verlet init

  Verlet::init();

}

/* ----------------------------------------------------------------------
   setup before run
   kproc partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplitRK::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fputs("Setting up Verlet/split/rk run ...\n",screen);
    if (flag) {
      utils::print(screen,"  Unit style    : {}\n"
                   "  Current step  : {}\n"
                   "  Time step     : {}\n",
                   update->unit_style,update->ntimestep,update->dt);
      timer->print_timeout(screen);
    }
  }

  if (lmp->kokkos)
    error->all(FLERR,"KOKKOS package requires run_style verlet/kk");

  rproc = (universe->iworld == 0) ? 1 : 0;
  if (rproc) {
    update->setupflag = 1;

    // setup domain, communication and neighboring
    // acquire ghosts
    // build neighbor lists

    atom->setup();
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    if (atom->sortfreq > 0) atom->sort();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;

    // compute all forces

    force->setup();
    ev_set(update->ntimestep);
    force_clear();
    modify->setup_pre_force(vflag);
  }//rproc

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag)
      force->kspace->r2k_comm(eflag,vflag);
  }
  if (rproc) {
    if (pair_compute_flag) force->pair->compute(eflag,vflag);
    else if (force->pair) force->pair->compute_dummy(eflag,vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
    }
  }//rproc

  if (force->kspace) {
    if (kspace_compute_flag) {
      if(!rproc) force->kspace->compute_grid_potentials(eflag,vflag);
      force->kspace->k2r_comm(eflag,vflag);
    }
    else force->kspace->compute_dummy(eflag,vflag);
  }

  if (rproc) {
    modify->setup_pre_reverse(eflag,vflag);
    if (force->newton) comm->reverse_comm();

    modify->setup(vflag);
    output->setup(flag);
    update->setupflag = 0;
  }//rproc
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
   kproc partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplitRK::setup_minimal(int flag)
{
  rproc = (universe->iworld == 0) ? 1 : 0;
  if (rproc) {
    update->setupflag = 1;

    // setup domain, communication and neighboring
    // acquire ghosts
    // build neighbor lists

    if (flag) {
      modify->setup_pre_exchange();
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
      comm->exchange();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      domain->image_check();
      domain->box_too_small_check();
      modify->setup_pre_neighbor();
      neighbor->build(1);
      modify->setup_post_neighbor();
      neighbor->ncalls = 0;
    }

    // compute all forces

    ev_set(update->ntimestep);
    force_clear();
    modify->setup_pre_force(vflag);
  }//rproc

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag)
      force->kspace->r2k_comm(eflag,vflag);
  }

  if (rproc) {
    if (pair_compute_flag) force->pair->compute(eflag,vflag);
    else if (force->pair) force->pair->compute_dummy(eflag,vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
    }
  }//rproc

  if (force->kspace) {
    if (kspace_compute_flag) {
      if(!rproc) force->kspace->compute_grid_potentials(eflag,vflag);
      force->kspace->k2r_comm(eflag,vflag);
    }
    else force->kspace->compute_dummy(eflag,vflag);
  }

  if (rproc) {
    modify->setup_pre_reverse(eflag,vflag);
    if (force->newton) comm->reverse_comm();
    modify->setup(vflag);
    update->setupflag = 0;
  }//rproc
}

/* ----------------------------------------------------------------------
   run for N steps
   rproc partition does everything but Kspace
   kproc partition does just Kspace
   communicate back and forth every step:
     atom coords from rproc -> kproc
     kspace forces from kproc -> rproc
     also box bounds from rproc -> kproc if necessary
------------------------------------------------------------------------- */

void VerletSplitRK::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;
  // sync both partitions before start timer

  MPI_Barrier(universe->uworld);
  timer->init();
  timer->barrier_start();

  // setup initial Rspace <-> Kspace comm params
  //rk_setup();

  // flags for timestepping iterations

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force_any;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    if (rproc) {
      timer->stamp();
      modify->initial_integrate(vflag);
      if (n_post_integrate) modify->post_integrate();
      timer->stamp(Timer::MODIFY);

      /*Atom-specific information is maintained on R-processes,
        reneighboring only occurs for R-processes*/
      timer->stamp(Timer::COMM);
      nflag = neighbor->decide();
      timer->stamp(Timer::NEIGH);
      if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
      } else {
        if (n_pre_exchange) modify->pre_exchange();
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
          domain->reset_box();
          comm->setup();
          if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (sortflag && ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (n_pre_neighbor) modify->pre_neighbor();
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
      }
    }

    // if reneighboring occurred, re-setup Rspace <-> Kspace comm params
    // comm Rspace atom coords to Kspace procs

    //Call to VerletSplit::rk_setup() no longer necessary
    timer->stamp();
    force->kspace->r2k_comm(eflag,vflag);
    timer->stamp(Timer::KSPACE);

    // force computations
    if (rproc) {
      timer->stamp();
      force_clear(); //Only needs to be called by R-processes
      if (n_pre_force) modify->pre_force(vflag);
      timer->stamp(Timer::MODIFY);


      timer->stamp();
      if (force->pair) {
        force->pair->compute(eflag,vflag);
        timer->stamp(Timer::PAIR);
      }

      if (atom->molecular != Atom::ATOMIC) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
      }

    } else {
      timer->stamp();
      force->kspace->compute_grid_potentials(eflag,vflag);
      timer->stamp(Timer::KSPACE);
    }

    // comm and sum Kspace forces back to Rspace procs

    timer->stamp();
    force->kspace->k2r_comm(eflag,vflag);
    timer->stamp(Timer::KSPACE);
    if(rproc) {
      /*pre_reverse and reverse_comm only need to be addressed for R-processes
        after kspace forces have been accumulated.*/
      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }
      if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }
      // force modifications, final time integration, diagnostics
      // all output
      timer->stamp();
      if (n_post_force) modify->post_force(vflag);
      modify->final_integrate();
      if (n_end_of_step) modify->end_of_step();
      timer->stamp(Timer::MODIFY);

      if (ntimestep == output->next) {
        timer->stamp();
        output->write(ntimestep);
        timer->stamp(Timer::OUTPUT);
      }
    } //rproc
  } //for i
}


