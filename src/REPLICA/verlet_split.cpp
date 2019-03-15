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
   Contributing authors: Yuxing Peng and Chris Knight (U Chicago)
------------------------------------------------------------------------- */

#include <cstring>
#include "verlet_split.h"
#include "universe.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "fix.h"
#include "modify.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletSplit::VerletSplit(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg), qsize(NULL), qdisp(NULL), xsize(NULL), xdisp(NULL), f_kspace(NULL)
{
  // error checks on partitions

  if (universe->nworlds != 2)
    error->universe_all(FLERR,"Verlet/split requires 2 partitions");
  if (universe->procs_per_world[0] % universe->procs_per_world[1])
    error->universe_all(FLERR,"Verlet/split requires Rspace partition "
                        "size be multiple of Kspace partition size");
  if (comm->style != 0)
    error->universe_all(FLERR,"Verlet/split can only currently be used with "
                        "comm_style brick");

  // master = 1 for Rspace procs, 0 for Kspace procs

  if (universe->iworld == 0) master = 1;
  else master = 0;

  ratio = universe->procs_per_world[0] / universe->procs_per_world[1];

  // Kspace root proc broadcasts info about Kspace proc layout to Rspace procs

  int kspace_procgrid[3];

  if (universe->me == universe->root_proc[1]) {
    kspace_procgrid[0] = comm->procgrid[0];
    kspace_procgrid[1] = comm->procgrid[1];
    kspace_procgrid[2] = comm->procgrid[2];
  }
  MPI_Bcast(kspace_procgrid,3,MPI_INT,universe->root_proc[1],universe->uworld);

  int ***kspace_grid2proc;
  memory->create(kspace_grid2proc,kspace_procgrid[0],
                 kspace_procgrid[1],kspace_procgrid[2],
                 "verlet/split:kspace_grid2proc");

  if (universe->me == universe->root_proc[1]) {
    for (int i = 0; i < comm->procgrid[0]; i++)
      for (int j = 0; j < comm->procgrid[1]; j++)
        for (int k = 0; k < comm->procgrid[2]; k++)
          kspace_grid2proc[i][j][k] = comm->grid2proc[i][j][k];
  }
  MPI_Bcast(&kspace_grid2proc[0][0][0],
            kspace_procgrid[0]*kspace_procgrid[1]*kspace_procgrid[2],MPI_INT,
            universe->root_proc[1],universe->uworld);

  // Rspace partition must be multiple of Kspace partition in each dim
  // so atoms of one Kspace proc coincide with atoms of several Rspace procs

  if (master) {
    int flag = 0;
    if (comm->procgrid[0] % kspace_procgrid[0]) flag = 1;
    if (comm->procgrid[1] % kspace_procgrid[1]) flag = 1;
    if (comm->procgrid[2] % kspace_procgrid[2]) flag = 1;
    if (flag)
      error->one(FLERR,
                 "Verlet/split requires Rspace partition layout be "
                 "multiple of Kspace partition layout in each dim");
  }

  // block = 1 Kspace proc with set of Rspace procs it overlays
  // me_block = 0 for Kspace proc
  // me_block = 1 to ratio for Rspace procs
  // block = MPI communicator for that set of procs

  int iblock,key;

  if (!master) {
    iblock = comm->me;
    key = 0;
  } else {
    int kpx = comm->myloc[0] / (comm->procgrid[0]/kspace_procgrid[0]);
    int kpy = comm->myloc[1] / (comm->procgrid[1]/kspace_procgrid[1]);
    int kpz = comm->myloc[2] / (comm->procgrid[2]/kspace_procgrid[2]);
    iblock = kspace_grid2proc[kpx][kpy][kpz];
    key = 1;
  }

  MPI_Comm_split(universe->uworld,iblock,key,&block);
  MPI_Comm_rank(block,&me_block);

  // output block groupings to universe screen/logfile
  // bmap is ordered by block and then by proc within block

  int *bmap = new int[universe->nprocs];
  for (int i = 0; i < universe->nprocs; i++) bmap[i] = -1;
  bmap[iblock*(ratio+1)+me_block] = universe->me;

  int *bmapall = new int[universe->nprocs];
  MPI_Allreduce(bmap,bmapall,universe->nprocs,MPI_INT,MPI_MAX,universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,
              "Per-block Rspace/Kspace proc IDs (original proc IDs):\n");
      int m = 0;
      for (int i = 0; i < universe->nprocs/(ratio+1); i++) {
        fprintf(universe->uscreen,"  block %d:",i);
        int kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++)
          fprintf(universe->uscreen," %d",bmapall[m+j]);
        fprintf(universe->uscreen," %d",kspace_proc);
        kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++) {
          if (j == 1) fprintf(universe->uscreen," (");
          else fprintf(universe->uscreen," ");
          fprintf(universe->uscreen,"%d",
                  universe->uni2orig[bmapall[m+j]]);
        }
        fprintf(universe->uscreen," %d)\n",universe->uni2orig[kspace_proc]);
        m += ratio + 1;
      }
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,
              "Per-block Rspace/Kspace proc IDs (original proc IDs):\n");
      int m = 0;
      for (int i = 0; i < universe->nprocs/(ratio+1); i++) {
        fprintf(universe->ulogfile,"  block %d:",i);
        int kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++)
          fprintf(universe->ulogfile," %d",bmapall[m+j]);

        fprintf(universe->ulogfile," %d",kspace_proc);
        kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++) {
          if (j == 1) fprintf(universe->ulogfile," (");
          else fprintf(universe->ulogfile," ");
          fprintf(universe->ulogfile,"%d",
                  universe->uni2orig[bmapall[m+j]]);
        }
        fprintf(universe->ulogfile," %d)\n",universe->uni2orig[kspace_proc]);
        m += ratio + 1;
      }
    }
  }

  memory->destroy(kspace_grid2proc);
  delete [] bmap;
  delete [] bmapall;

  // size/disp = vectors for MPI gather/scatter within block

  qsize = new int[ratio+1];
  qdisp = new int[ratio+1];
  xsize = new int[ratio+1];
  xdisp = new int[ratio+1];

  // f_kspace = Rspace copy of Kspace forces
  // allocate dummy version for Kspace partition

  maxatom = 0;
  f_kspace = NULL;
  if (!master) memory->create(f_kspace,1,1,"verlet/split:f_kspace");
}

/* ---------------------------------------------------------------------- */

VerletSplit::~VerletSplit()
{
  delete [] qsize;
  delete [] qdisp;
  delete [] xsize;
  delete [] xdisp;
  memory->destroy(f_kspace);
  MPI_Comm_free(&block);
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletSplit::init()
{
  if (comm->style != 0)
    error->universe_all(FLERR,"Verlet/split can only currently be used with "
                        "comm_style brick");
  if (!force->kspace && comm->me == 0)
    error->warning(FLERR,"No Kspace calculation with verlet/split");

  if (force->kspace_match("/tip4p",0)) tip4p_flag = 1;
  else tip4p_flag = 0;

  // currently TIP4P does not work with verlet/split, so generate error
  // see Axel email on this, also other TIP4P notes below

  if (tip4p_flag) error->all(FLERR,"Verlet/split does not yet support TIP4P");

  Verlet::init();
}

/* ----------------------------------------------------------------------
   setup before run
   servant partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplit::setup(int flag)
{
  if (comm->me == 0 && screen)
    fprintf(screen,"Setting up Verlet/split run ...\n");

  if (!master) force->kspace->setup();
  else Verlet::setup(flag);
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
   servant partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplit::setup_minimal(int flag)
{
  if (!master) force->kspace->setup();
  else Verlet::setup_minimal(flag);
}

/* ----------------------------------------------------------------------
   run for N steps
   master partition does everything but Kspace
   servant partition does just Kspace
   communicate back and forth every step:
     atom coords from master -> servant
     kspace forces from servant -> master
     also box bounds from master -> servant if necessary
------------------------------------------------------------------------- */

void VerletSplit::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  // sync both partitions before start timer

  MPI_Barrier(universe->uworld);
  timer->init();
  timer->barrier_start();

  // setup initial Rspace <-> Kspace comm params

  rk_setup();

  // check if OpenMP support fix defined

  Fix *fix_omp;
  int ifix = modify->find_fix("package_omp");
  if (ifix < 0) fix_omp = NULL;
  else fix_omp = modify->fix[ifix];

  // flags for timestepping iterations

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    if (master) {
      modify->initial_integrate(vflag);
      if (n_post_integrate) modify->post_integrate();
    }

    // regular communication vs neighbor list rebuild

    if (master) nflag = neighbor->decide();
    MPI_Bcast(&nflag,1,MPI_INT,1,block);

    if (master) {
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

    if (nflag) rk_setup();
    r2k_comm();

    // force computations

    force_clear();

    if (master) {
      if (n_pre_force) modify->pre_force(vflag);

      timer->stamp();
      if (force->pair) {
        force->pair->compute(eflag,vflag);
        timer->stamp(Timer::PAIR);
      }

      if (atom->molecular) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
      }

      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }
      if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }

    } else {

      // run FixOMP as sole pre_force fix, if defined

      if (fix_omp) fix_omp->pre_force(vflag);

      if (force->kspace) {
        timer->stamp();
        force->kspace->compute(eflag,vflag);
        timer->stamp(Timer::KSPACE);
      }

      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }

      // TIP4P PPPM puts forces on ghost atoms, so must reverse_comm()

      if (tip4p_flag && force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }
    }

    // comm and sum Kspace forces back to Rspace procs

    k2r_comm();

    // force modifications, final time integration, diagnostics
    // all output

    if (master) {
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
    }
  }
}

/* ----------------------------------------------------------------------
   setup params for Rspace <-> Kspace communication
   called initially and after every reneighbor
   also communcicate atom charges from Rspace to KSpace since static
------------------------------------------------------------------------- */

void VerletSplit::rk_setup()
{
  // grow f_kspace array on master procs if necessary

  if (master) {
    if (atom->nmax > maxatom) {
      memory->destroy(f_kspace);
      maxatom = atom->nmax;
      memory->create(f_kspace,maxatom,3,"verlet/split:f_kspace");
    }
  }

  // qsize = # of atoms owned by each master proc in block

  int n = 0;
  if (master) n = atom->nlocal;
  MPI_Gather(&n,1,MPI_INT,qsize,1,MPI_INT,0,block);

  // setup qdisp, xsize, xdisp based on qsize
  // only needed by Kspace proc
  // set Kspace nlocal to sum of Rspace nlocals
  // insure Kspace atom arrays are large enough

  if (!master) {
    qsize[0] = qdisp[0] = xsize[0] = xdisp[0] = 0;
    for (int i = 1; i <= ratio; i++) {
      qdisp[i] = qdisp[i-1]+qsize[i-1];
      xsize[i] = 3*qsize[i];
      xdisp[i] = xdisp[i-1]+xsize[i-1];
    }

    atom->nlocal = qdisp[ratio] + qsize[ratio];
    while (atom->nmax <= atom->nlocal) atom->avec->grow(0);
    atom->nghost = 0;
  }

  // one-time gather of Rspace atom charges to Kspace proc

  MPI_Gatherv(atom->q,n,MPI_DOUBLE,atom->q,qsize,qdisp,MPI_DOUBLE,0,block);

  // for TIP4P also need to send atom type and tag
  // KSpace procs need to acquire ghost atoms and map all their atoms
  // map_clear() call is in lieu of comm->exchange() which performs map_clear
  // borders() call acquires ghost atoms and maps them
  // NOTE: do atom coords need to be communicated here before borders() call?
  //   could do this by calling r2k_comm() here and not again from run()
  //   except that forward_comm() in r2k_comm() is wrong

  if (tip4p_flag) {
    //r2k_comm();
    MPI_Gatherv(atom->type,n,MPI_INT,atom->type,qsize,qdisp,MPI_INT,0,block);
    MPI_Gatherv(atom->tag,n,MPI_LMP_TAGINT,
                atom->tag,qsize,qdisp,MPI_LMP_TAGINT,0,block);
    if (!master) {
      if (triclinic) domain->x2lamda(atom->nlocal);
      if (domain->box_change) comm->setup();
      timer->stamp();
      atom->map_clear();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      timer->stamp(Timer::COMM);
    }
  }
}

/* ----------------------------------------------------------------------
   communicate Rspace atom coords to Kspace
   also eflag,vflag and box bounds if needed
------------------------------------------------------------------------- */

void VerletSplit::r2k_comm()
{
  int n = 0;
  if (master) n = atom->nlocal;
  MPI_Gatherv(atom->x[0],n*3,MPI_DOUBLE,atom->x[0],xsize,xdisp,
              MPI_DOUBLE,0,block);

  // send eflag,vflag from Rspace to Kspace

  if (me_block == 1) {
    int flags[2];
    flags[0] = eflag; flags[1] = vflag;
    MPI_Send(flags,2,MPI_INT,0,0,block);
  } else if (!master) {
    int flags[2];
    MPI_Recv(flags,2,MPI_INT,1,0,block,MPI_STATUS_IGNORE);
    eflag = flags[0]; vflag = flags[1];
  }

  // send box bounds from Rspace to Kspace if simulation box is dynamic

  if (domain->box_change) {
    if (me_block == 1) {
      MPI_Send(domain->boxlo,3,MPI_DOUBLE,0,0,block);
      MPI_Send(domain->boxhi,3,MPI_DOUBLE,0,0,block);
    } else if (!master) {
      MPI_Recv(domain->boxlo,3,MPI_DOUBLE,1,0,block,MPI_STATUS_IGNORE);
      MPI_Recv(domain->boxhi,3,MPI_DOUBLE,1,0,block,MPI_STATUS_IGNORE);
      domain->set_global_box();
      domain->set_local_box();
      force->kspace->setup();
    }
  }

  // for TIP4P, Kspace partition needs to update its ghost atoms

  if (tip4p_flag && !master) {
    timer->stamp();
    comm->forward_comm();
    timer->stamp(Timer::COMM);
  }
}

/* ----------------------------------------------------------------------
   communicate and sum Kspace atom forces back to Rspace
------------------------------------------------------------------------- */

void VerletSplit::k2r_comm()
{
  if (eflag) MPI_Bcast(&force->kspace->energy,1,MPI_DOUBLE,0,block);
  if (vflag) MPI_Bcast(force->kspace->virial,6,MPI_DOUBLE,0,block);

  int n = 0;
  if (master) n = atom->nlocal;
  MPI_Scatterv(atom->f[0],xsize,xdisp,MPI_DOUBLE,
               f_kspace[0],n*3,MPI_DOUBLE,0,block);

  if (master) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += f_kspace[i][0];
      f[i][1] += f_kspace[i][1];
      f[i][2] += f_kspace[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of Kspace force array on master procs
------------------------------------------------------------------------- */

bigint VerletSplit::memory_usage()
{
  bigint bytes = maxatom*3 * sizeof(double);
  return bytes;
}
