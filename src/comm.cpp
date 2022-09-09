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

#include "comm.h"

#include "accelerator_kokkos.h"
#include "atom.h"               // IWYU pragma: keep
#include "atom_vec.h"
#include "bond.h"
#include "compute.h"
#include "domain.h"             // IWYU pragma: keep
#include "dump.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "irregular.h"
#include "memory.h"             // IWYU pragma: keep
#include "modify.h"
#include "neighbor.h"           // IWYU pragma: keep
#include "output.h"
#include "pair.h"
#include "procmap.h"
#include "universe.h"
#include "update.h"

#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define BUFEXTRA 1024

enum{ONELEVEL,TWOLEVEL,NUMA,CUSTOM};
enum{CART,CARTREORDER,XYZ};

/* ---------------------------------------------------------------------- */

Comm::Comm(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  mode = 0;
  bordergroup = 0;
  cutghostuser = 0.0;
  cutusermulti = nullptr;
  cutusermultiold = nullptr;
  ncollections = 0;
  ncollections_cutoff = 0;
  ghost_velocity = 0;

  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  coregrid[0] = coregrid[1] = coregrid[2] = 1;
  gridflag = ONELEVEL;
  mapflag = CART;
  customfile = nullptr;
  outfile = nullptr;
  recv_from_partition = send_to_partition = -1;
  otherflag = 0;

  maxexchange = maxexchange_atom = maxexchange_fix = 0;
  maxexchange_fix_dynamic = 0;
  bufextra = BUFEXTRA;

  grid2proc = nullptr;
  xsplit = ysplit = zsplit = nullptr;
  rcbnew = 0;
  multi_reduce = 0;

  // use of OpenMP threads
  // query OpenMP for number of threads/process set by user at run-time
  // if the OMP_NUM_THREADS environment variable is not set, we default
  // to using 1 thread. This follows the principle of the least surprise,
  // while practically all OpenMP implementations violate it by using
  // as many threads as there are (virtual) CPU cores by default.

  nthreads = 1;
#ifdef _OPENMP
  if (lmp->kokkos) {
    nthreads = lmp->kokkos->nthreads * lmp->kokkos->numa;
  } else if (getenv("OMP_NUM_THREADS") == nullptr) {
    nthreads = 1;
    if (me == 0)
      error->message(FLERR,"OMP_NUM_THREADS environment is not set. "
                           "Defaulting to 1 thread.");
  } else {
    nthreads = omp_get_max_threads();
  }

  // enforce consistent number of threads across all MPI tasks

  MPI_Bcast(&nthreads,1,MPI_INT,0,world);
  if (!lmp->kokkos) omp_set_num_threads(nthreads);

  if (me == 0)
    utils::logmesg(lmp,"  using {} OpenMP thread(s) per MPI task\n",nthreads);
#endif

}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  memory->destroy(grid2proc);
  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);
  memory->destroy(cutusermulti);
  memory->destroy(cutusermultiold);
  delete [] customfile;
  delete [] outfile;
}

/* ----------------------------------------------------------------------
   deep copy of arrays from old Comm class to new one
   all public/protected vectors/arrays in parent Comm class must be copied
   called from alternate constructor of child classes
   when new comm style is created from Input
------------------------------------------------------------------------- */

void Comm::copy_arrays(Comm *oldcomm)
{
  if (oldcomm->grid2proc) {
    memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
                   "comm:grid2proc");
    memcpy(&grid2proc[0][0][0],&oldcomm->grid2proc[0][0][0],
           (procgrid[0]*procgrid[1]*procgrid[2])*sizeof(int));

    memory->create(xsplit,procgrid[0]+1,"comm:xsplit");
    memory->create(ysplit,procgrid[1]+1,"comm:ysplit");
    memory->create(zsplit,procgrid[2]+1,"comm:zsplit");
    memcpy(xsplit,oldcomm->xsplit,(procgrid[0]+1)*sizeof(double));
    memcpy(ysplit,oldcomm->ysplit,(procgrid[1]+1)*sizeof(double));
    memcpy(zsplit,oldcomm->zsplit,(procgrid[2]+1)*sizeof(double));
  }

  ncollections = oldcomm->ncollections;
  ncollections_cutoff = oldcomm->ncollections_cutoff;
  if (oldcomm->cutusermulti) {
    memory->create(cutusermulti,ncollections_cutoff,"comm:cutusermulti");
    memcpy(cutusermulti,oldcomm->cutusermulti,ncollections_cutoff);
  }

  if (oldcomm->cutusermultiold) {
    memory->create(cutusermultiold,atom->ntypes+1,"comm:cutusermultiold");
    memcpy(cutusermultiold,oldcomm->cutusermultiold,atom->ntypes+1);
  }

  if (customfile)
    customfile = utils::strdup(oldcomm->customfile);

  if (outfile)
    outfile = utils::strdup(oldcomm->outfile);
}

/* ----------------------------------------------------------------------
   common to all Comm styles
------------------------------------------------------------------------- */

void Comm::init()
{
  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // check warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()
  // really should check every exchange() in case box size is shrinking
  //   but seems overkill to do that (fix balance does perform this check)

  domain->subbox_too_small_check(neighbor->skin);

  // comm_only = 1 if only x,f are exchanged in forward/reverse comm
  // comm_x_only = 0 if ghost_velocity since velocities are added

  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;
  if (ghost_velocity) comm_x_only = 0;

  // set per-atom sizes for forward/reverse/border comm
  // augment by velocity and fix quantities if needed

  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;

  if (ghost_velocity) size_forward += atom->avec->size_velocity;
  if (ghost_velocity) size_border += atom->avec->size_velocity;

  const auto &fix_list = modify->get_fix_list();
  for (const auto &fix : fix_list)
    size_border += fix->comm_border;

  // per-atom limits for communication
  // maxexchange = max # of datums in exchange comm, set in exchange()
  // maxforward = # of datums in largest forward comm
  // maxreverse = # of datums in largest reverse comm
  // query pair,fix,compute,dump for their requirements
  // pair style can force reverse comm even if newton off

  maxforward = MAX(size_forward,size_border);
  maxreverse = size_reverse;

  if (force->pair) maxforward = MAX(maxforward,force->pair->comm_forward);
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse);

  if (force->bond) maxforward = MAX(maxforward,force->bond->comm_forward);
  if (force->bond) maxreverse = MAX(maxreverse,force->bond->comm_reverse);

  for (const auto &fix : fix_list) {
    maxforward = MAX(maxforward, fix->comm_forward);
    maxreverse = MAX(maxreverse, fix->comm_reverse);
  }

  for (const auto &compute : modify->get_compute_list()) {
    maxforward = MAX(maxforward,compute->comm_forward);
    maxreverse = MAX(maxreverse,compute->comm_reverse);
  }

  for (const auto &dump: output->get_dump_list()) {
    maxforward = MAX(maxforward,dump->comm_forward);
    maxreverse = MAX(maxreverse,dump->comm_reverse);
  }

  if (force->newton == 0) maxreverse = 0;
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse_off);
  if (force->bond) maxreverse = MAX(maxreverse,force->bond->comm_reverse_off);

  // maxexchange_atom = size of an exchanged atom, set by AtomVec
  //   only needs to be set if size > BUFEXTRA
  // maxexchange_fix_dynamic = 1 if any fix sets its maxexchange dynamically

  maxexchange_atom = atom->avec->maxexchange;

  maxexchange_fix_dynamic = 0;
  for (const auto &fix : fix_list) if (fix->maxexchange_dynamic) maxexchange_fix_dynamic = 1;

  if ((mode == Comm::MULTI) && (neighbor->style != Neighbor::MULTI))
    error->all(FLERR,"Cannot use comm mode multi without multi-style neighbor lists");

  if (multi_reduce) {
    if (force->newton == 0)
      error->all(FLERR,"Cannot use multi/reduce communication with Newton off");
    if (neighbor->any_full())
      error->all(FLERR,"Cannot use multi/reduce communication with a full neighbor list");
    if (mode != Comm::MULTI)
      error->all(FLERR,"Cannot use multi/reduce communication without mode multi");
  }
}

/* ----------------------------------------------------------------------
   set maxexchange based on AtomVec and fixes
------------------------------------------------------------------------- */

void Comm::init_exchange()
{
  maxexchange_fix = 0;
  for (const auto &fix : modify->get_fix_list()) maxexchange_fix += fix->maxexchange;

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
}

/* ----------------------------------------------------------------------
   modify communication params
   invoked from input script by comm_modify command
------------------------------------------------------------------------- */

void Comm::modify_params(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "comm_modify", error);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "comm_modify mode", error);
      if (strcmp(arg[iarg+1],"single") == 0) {
        // need to reset cutghostuser when switching comm mode
        if (mode == Comm::MULTI) cutghostuser = 0.0;
        if (mode == Comm::MULTIOLD) cutghostuser = 0.0;
        memory->destroy(cutusermulti);
        memory->destroy(cutusermultiold);
        mode = Comm::SINGLE;
      } else if (strcmp(arg[iarg+1],"multi") == 0) {
        if (neighbor->style != Neighbor::MULTI)
          error->all(FLERR,"Cannot use comm mode 'multi' without 'multi' style neighbor lists");
        // need to reset cutghostuser when switching comm mode
        if (mode == Comm::SINGLE) cutghostuser = 0.0;
        if (mode == Comm::MULTIOLD) cutghostuser = 0.0;
        memory->destroy(cutusermultiold);
        mode = Comm::MULTI;
      } else if (strcmp(arg[iarg+1],"multi/old") == 0) {
        if (neighbor->style == Neighbor::MULTI)
          error->all(FLERR,"Cannot use comm mode 'multi/old' with 'multi' style neighbor lists");
        // need to reset cutghostuser when switching comm mode
        if (mode == Comm::SINGLE) cutghostuser = 0.0;
        if (mode == Comm::MULTI) cutghostuser = 0.0;
        memory->destroy(cutusermulti);
        mode = Comm::MULTIOLD;
      } else error->all(FLERR,"Unknown comm_modify mode argument: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "comm_modify group", error);
      bordergroup = group->find(arg[iarg+1]);
      if (bordergroup < 0)
        error->all(FLERR, "Invalid comm_modify keyword: group {} not found", arg[iarg+1]);
      if (bordergroup && ((atom->firstgroupname == nullptr) || strcmp(arg[iarg+1],atom->firstgroupname) != 0))
        error->all(FLERR, "Comm_modify group != atom_modify first group: {}", atom->firstgroupname);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "comm_modify cutoff", error);
      if (mode == Comm::MULTI)
        error->all(FLERR, "Use cutoff/multi keyword to set cutoff in multi mode");
      if (mode == Comm::MULTIOLD)
        error->all(FLERR, "Use cutoff/multi/old keyword to set cutoff in multi mode");
      cutghostuser = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutghostuser < 0.0)
        error->all(FLERR,"Invalid cutoff {} in comm_modify command", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff/multi") == 0) {
      int i,nlo,nhi;
      double cut;
      if (mode == Comm::SINGLE)
        error->all(FLERR,"Use cutoff keyword to set cutoff in single mode");
      if (mode == Comm::MULTIOLD)
        error->all(FLERR,"Use cutoff/multi/old keyword to set cutoff in multi/old mode");
      if (domain->box_exist == 0)
        error->all(FLERR, "Cannot set cutoff/multi before simulation box is defined");

      // Check if # of collections has changed, if so erase any previously defined cutoffs
      // Neighbor will reset ncollections if collections are redefined
      if (! cutusermulti || ncollections_cutoff != neighbor->ncollections) {
        ncollections_cutoff = neighbor->ncollections;
        memory->destroy(cutusermulti);
        memory->create(cutusermulti,ncollections_cutoff,"comm:cutusermulti");
        for (i=0; i < ncollections_cutoff; ++i)
          cutusermulti[i] = -1.0;
      }
      utils::bounds(FLERR,arg[iarg+1],1,ncollections_cutoff,nlo,nhi,error);
      cut = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      cutghostuser = MAX(cutghostuser,cut);
      if (cut < 0.0)
        error->all(FLERR,"Invalid cutoff {} in comm_modify command", arg[iarg+2]);
      // collections use 1-based indexing externally and 0-based indexing internally
      for (i=nlo; i<=nhi; ++i)
        cutusermulti[i-1] = cut;
      iarg += 3;
    }  else if (strcmp(arg[iarg],"cutoff/multi/old") == 0) {
      int i,nlo,nhi;
      double cut;
      if (mode == Comm::SINGLE)
        error->all(FLERR,"Use cutoff keyword to set cutoff in single mode");
      if (mode == Comm::MULTI)
        error->all(FLERR,"Use cutoff/multi keyword to set cutoff in multi mode");
      if (domain->box_exist == 0)
        error->all(FLERR, "Cannot set cutoff/multi before simulation box is defined");
      const int ntypes = atom->ntypes;
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "comm_modify cutoff/multi/old", error);
      if (cutusermultiold == nullptr) {
        memory->create(cutusermultiold,ntypes+1,"comm:cutusermultiold");
        for (i=0; i < ntypes+1; ++i)
          cutusermultiold[i] = -1.0;
      }
      utils::bounds(FLERR,arg[iarg+1],1,ntypes,nlo,nhi,error);
      cut = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      cutghostuser = MAX(cutghostuser,cut);
      if (cut < 0.0)
        error->all(FLERR,"Invalid cutoff {} in comm_modify command", arg[iarg+2]);
      for (i=nlo; i<=nhi; ++i)
        cutusermultiold[i] = cut;
      iarg += 3;
    } else if (strcmp(arg[iarg],"reduce/multi") == 0) {
      if (mode == Comm::SINGLE)
        error->all(FLERR,"Use reduce/multi in mode multi only");
      multi_reduce = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "comm_modify vel", error);
      ghost_velocity = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown comm_modify keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   set dimensions for 3d grid of processors, and associated flags
   invoked from input script by processors command
------------------------------------------------------------------------- */

void Comm::set_processors(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal processors command");

  if (strcmp(arg[0],"*") == 0) user_procgrid[0] = 0;
  else user_procgrid[0] = utils::inumeric(FLERR,arg[0],false,lmp);
  if (strcmp(arg[1],"*") == 0) user_procgrid[1] = 0;
  else user_procgrid[1] = utils::inumeric(FLERR,arg[1],false,lmp);
  if (strcmp(arg[2],"*") == 0) user_procgrid[2] = 0;
  else user_procgrid[2] = utils::inumeric(FLERR,arg[2],false,lmp);

  if (user_procgrid[0] < 0 || user_procgrid[1] < 0 || user_procgrid[2] < 0)
    error->all(FLERR,"Illegal processors command");

  int p = user_procgrid[0]*user_procgrid[1]*user_procgrid[2];
  if (p && p != nprocs)
    error->all(FLERR,"Specified processors != physical processors");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");

      if (strcmp(arg[iarg+1],"onelevel") == 0) {
        gridflag = ONELEVEL;

      } else if (strcmp(arg[iarg+1],"twolevel") == 0) {
        if (iarg+6 > narg) error->all(FLERR,"Illegal processors command");
        gridflag = TWOLEVEL;

        ncores = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (strcmp(arg[iarg+3],"*") == 0) user_coregrid[0] = 0;
        else user_coregrid[0] = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
        if (strcmp(arg[iarg+4],"*") == 0) user_coregrid[1] = 0;
        else user_coregrid[1] = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
        if (strcmp(arg[iarg+5],"*") == 0) user_coregrid[2] = 0;
        else user_coregrid[2] = utils::inumeric(FLERR,arg[iarg+5],false,lmp);

        if (ncores <= 0 || user_coregrid[0] < 0 ||
            user_coregrid[1] < 0 || user_coregrid[2] < 0)
          error->all(FLERR,"Illegal processors command");
        iarg += 4;

      } else if (strcmp(arg[iarg+1],"numa") == 0) {
        gridflag = NUMA;

      } else if (strcmp(arg[iarg+1],"custom") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal processors command");
        gridflag = CUSTOM;
        delete [] customfile;
        customfile = utils::strdup(arg[iarg+2]);
        iarg += 1;

      } else error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"map") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      if (strcmp(arg[iarg+1],"cart") == 0) mapflag = CART;
      else if (strcmp(arg[iarg+1],"cart/reorder") == 0) mapflag = CARTREORDER;
      else if (strcmp(arg[iarg+1],"xyz") == 0 ||
               strcmp(arg[iarg+1],"xzy") == 0 ||
               strcmp(arg[iarg+1],"yxz") == 0 ||
               strcmp(arg[iarg+1],"yzx") == 0 ||
               strcmp(arg[iarg+1],"zxy") == 0 ||
               strcmp(arg[iarg+1],"zyx") == 0) {
        mapflag = XYZ;
        strncpy(xyz,arg[iarg+1],3);
      } else error->all(FLERR,"Illegal processors command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"part") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal processors command");
      if (universe->nworlds == 1)
        error->all(FLERR,
                   "Cannot use processors part command "
                   "without using partitions");
      int isend = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      int irecv = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (isend < 1 || isend > universe->nworlds ||
          irecv < 1 || irecv > universe->nworlds || isend == irecv)
        error->all(FLERR,"Invalid partitions in processors part command");
      if (isend-1 == universe->iworld) {
        if (send_to_partition >= 0)
          error->all(FLERR,
                     "Sending partition in processors part command "
                     "is already a sender");
        send_to_partition = irecv-1;
      }
      if (irecv-1 == universe->iworld) {
        if (recv_from_partition >= 0)
          error->all(FLERR,
                     "Receiving partition in processors part command "
                     "is already a receiver");
        recv_from_partition = isend-1;
      }

      // only receiver has otherflag dependency

      if (strcmp(arg[iarg+3],"multiple") == 0) {
        if (universe->iworld == irecv-1) {
          otherflag = 1;
          other_style = Comm::MULTIPLE;
        }
      } else error->all(FLERR,"Illegal processors command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal processors command");
      delete [] outfile;
      outfile = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal processors command");
  }

  // error checks

  if (gridflag == NUMA && mapflag != CART)
    error->all(FLERR,"Processors grid numa and map style are incompatible");
  if (otherflag && (gridflag == NUMA || gridflag == CUSTOM))
    error->all(FLERR,
               "Processors part option and grid style are incompatible");
}

/* ----------------------------------------------------------------------
   create a 3d grid of procs based on Nprocs and box size & shape
   map processors to grid, setup xyz split for a uniform grid
------------------------------------------------------------------------- */

void Comm::set_proc_grid(int outflag)
{
  // recv 3d proc grid of another partition if my 3d grid depends on it

  if (recv_from_partition >= 0) {
    if (me == 0) {
      MPI_Recv(other_procgrid,3,MPI_INT,
               universe->root_proc[recv_from_partition],0,
               universe->uworld,MPI_STATUS_IGNORE);
      MPI_Recv(other_coregrid,3,MPI_INT,
               universe->root_proc[recv_from_partition],0,
               universe->uworld,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(other_procgrid,3,MPI_INT,0,world);
    MPI_Bcast(other_coregrid,3,MPI_INT,0,world);
  }

  // create ProcMap class to create 3d grid and map procs to it

  auto pmap = new ProcMap(lmp);

  // create 3d grid of processors
  // produces procgrid and coregrid (if relevant)

  if (gridflag == ONELEVEL) {
    pmap->onelevel_grid(nprocs,user_procgrid,procgrid,
                        otherflag,other_style,other_procgrid,other_coregrid);

  } else if (gridflag == TWOLEVEL) {
    pmap->twolevel_grid(nprocs,user_procgrid,procgrid,
                        ncores,user_coregrid,coregrid,
                        otherflag,other_style,other_procgrid,other_coregrid);

  } else if (gridflag == NUMA) {
    pmap->numa_grid(nprocs,user_procgrid,procgrid,coregrid);

  } else if (gridflag == CUSTOM) {
    pmap->custom_grid(customfile,nprocs,user_procgrid,procgrid);
  }

  // error check on procgrid
  // should not be necessary due to ProcMap

  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs)
    error->all(FLERR,"Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all(FLERR,"Processor count in z must be 1 for 2d simulation");

  // grid2proc[i][j][k] = proc that owns i,j,k location in 3d grid

  if (grid2proc) memory->destroy(grid2proc);
  memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
                 "comm:grid2proc");

  // map processor IDs to 3d processor grid
  // produces myloc, procneigh, grid2proc

  if (gridflag == ONELEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0,procgrid,myloc,procneigh,grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1,procgrid,myloc,procneigh,grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz,procgrid,myloc,procneigh,grid2proc);

  } else if (gridflag == TWOLEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz,procgrid,ncores,coregrid,myloc,procneigh,grid2proc);

  } else if (gridflag == NUMA) {
    pmap->numa_map(0,coregrid,myloc,procneigh,grid2proc);

  } else if (gridflag == CUSTOM) {
    pmap->custom_map(procgrid,myloc,procneigh,grid2proc);
  }

  // print 3d grid info to screen and logfile

  if (outflag && me == 0) {
    auto mesg = fmt::format("  {} by {} by {} MPI processor grid\n",
                            procgrid[0],procgrid[1],procgrid[2]);
    if (gridflag == NUMA || gridflag == TWOLEVEL)
      mesg += fmt::format("  {} by {} by {} core grid within node\n",
                          coregrid[0],coregrid[1],coregrid[2]);
    utils::logmesg(lmp,mesg);
  }

  // print 3d grid details to outfile

  if (outfile) pmap->output(outfile,procgrid,grid2proc);

  // free ProcMap class

  delete pmap;

  // set xsplit,ysplit,zsplit for uniform spacings

  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);

  memory->create(xsplit,procgrid[0]+1,"comm:xsplit");
  memory->create(ysplit,procgrid[1]+1,"comm:ysplit");
  memory->create(zsplit,procgrid[2]+1,"comm:zsplit");

  for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
  for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
  for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];

  xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;

  // set lamda box params after procs are assigned
  // only set once unless load-balancing occurs

  if (domain->triclinic) domain->set_lamda_box();

  // send my 3d proc grid to another partition if requested

  if (send_to_partition >= 0) {
    if (me == 0) {
      MPI_Send(procgrid,3,MPI_INT,
               universe->root_proc[send_to_partition],0,
               universe->uworld);
      MPI_Send(coregrid,3,MPI_INT,
               universe->root_proc[send_to_partition],0,
               universe->uworld);
    }
  }
}

/* ----------------------------------------------------------------------
   determine suitable communication cutoff.
   this uses three inputs: 1) maximum neighborlist cutoff, 2) an estimate
   based on bond lengths and bonded interaction styles present, and 3) a
   user supplied communication cutoff.
   the neighbor list cutoff (1) is *always* used, since it is a requirement
   for neighborlists working correctly. the bond length based cutoff is
   *only* used, if no pair style is defined and no user cutoff is provided.
   otherwise, a warning is printed. if the bond length based estimate is
   larger than what is used.
   print a warning, if a user specified communication cutoff is overridden.
------------------------------------------------------------------------- */

double Comm::get_comm_cutoff()
{
  double maxcommcutoff, maxbondcutoff = 0.0;

  if (force->bond) {
    int n = atom->nbondtypes;
    for (int i = 1; i <= n; ++i)
      maxbondcutoff = MAX(maxbondcutoff,force->bond->equilibrium_distance(i));

    // apply bond length based heuristics.

    if (force->newton_bond) {
      if (force->dihedral || force->improper) {
        maxbondcutoff *= 2.25;
      } else {
        maxbondcutoff *=1.5;
      }
    } else {
      if (force->dihedral || force->improper) {
        maxbondcutoff *= 3.125;
      } else if (force->angle) {
        maxbondcutoff *= 2.25;
      } else {
        maxbondcutoff *=1.5;
      }
    }
    maxbondcutoff += neighbor->skin;
  }

  // always take the larger of max neighbor list and user specified cutoff

  maxcommcutoff = MAX(cutghostuser,neighbor->cutneighmax);

  // use cutoff estimate from bond length only if no user specified
  // cutoff was given and no pair style present. Otherwise print a
  // warning, if the estimated bond based cutoff is larger than what
  // is currently used.

  if (!force->pair && (cutghostuser == 0.0)) {
    maxcommcutoff = MAX(maxcommcutoff,maxbondcutoff);
  } else {
    if ((me == 0) && (maxbondcutoff > maxcommcutoff))
      error->warning(FLERR,"Communication cutoff {} is shorter than a bond "
                     "length based estimate of {}. This may lead to errors.",
                     maxcommcutoff,maxbondcutoff);
  }

  // print warning if neighborlist cutoff overrides user cutoff

  if ((me == 0) && (update->setupflag == 1)) {
    if ((cutghostuser > 0.0) && (maxcommcutoff > cutghostuser))
      error->warning(FLERR,"Communication cutoff adjusted to {}",maxcommcutoff);
  }

  // check maximum interval size for neighbor multi

  if (neighbor->interval_collection_flag) {
    for (int i = 0; i < neighbor->ncollections; i++){
      maxcommcutoff = MAX(maxcommcutoff, neighbor->collection2cut[i]);
    }
  }

  return maxcommcutoff;
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with coord x[3] based on current decomp
   x will be in box (orthogonal) or lamda coords (triclinic)
   if layout = UNIFORM, calculate owning proc directly
   if layout = NONUNIFORM, iteratively find owning proc via binary search
   if layout = TILED, CommTiled has its own method
   return owning proc ID via grid2proc
   return igx,igy,igz = logical grid loc of owing proc within 3d grid of procs
------------------------------------------------------------------------- */

int Comm::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  double *prd = domain->prd;
  double *boxlo = domain->boxlo;

  // initialize triclinic b/c coord2proc can be called before Comm::init()
  // via Irregular::migrate_atoms()

  triclinic = domain->triclinic;

  if (layout == Comm::LAYOUT_UNIFORM) {
    if (triclinic == 0) {
      igx = static_cast<int> (procgrid[0] * (x[0]-boxlo[0]) / prd[0]);
      igy = static_cast<int> (procgrid[1] * (x[1]-boxlo[1]) / prd[1]);
      igz = static_cast<int> (procgrid[2] * (x[2]-boxlo[2]) / prd[2]);
    } else {
      igx = static_cast<int> (procgrid[0] * x[0]);
      igy = static_cast<int> (procgrid[1] * x[1]);
      igz = static_cast<int> (procgrid[2] * x[2]);
    }

  } else if (layout == Comm::LAYOUT_NONUNIFORM) {
    if (triclinic == 0) {
      igx = utils::binary_search((x[0]-boxlo[0])/prd[0],procgrid[0],xsplit);
      igy = utils::binary_search((x[1]-boxlo[1])/prd[1],procgrid[1],ysplit);
      igz = utils::binary_search((x[2]-boxlo[2])/prd[2],procgrid[2],zsplit);
    } else {
      igx = utils::binary_search(x[0],procgrid[0],xsplit);
      igy = utils::binary_search(x[1],procgrid[1],ysplit);
      igz = utils::binary_search(x[2],procgrid[2],zsplit);
    }
  }

  if (igx < 0) igx = 0;
  if (igx >= procgrid[0]) igx = procgrid[0] - 1;
  if (igy < 0) igy = 0;
  if (igy >= procgrid[1]) igy = procgrid[1] - 1;
  if (igz < 0) igz = 0;
  if (igz >= procgrid[2]) igz = procgrid[2] - 1;

  return grid2proc[igx][igy][igz];
}

/* ----------------------------------------------------------------------
   partition a global regular grid into one brick-shaped sub-grid per proc
   if grid point is inside my sub-domain I own it,
     this includes sub-domain lo boundary but excludes hi boundary
   nx,ny,nz = extent of global grid
     indices into the global grid range from 0 to N-1 in each dim
   zfactor = 0.0 if the grid exactly covers the simulation box
   zfactor > 1.0 if the grid extends beyond the +z boundary by this factor
     used by 2d slab-mode PPPM
     this effectively maps proc sub-grids to a smaller subset of the grid
   nxyz lo/hi = inclusive lo/hi bounds of global grid sub-brick I own
   if proc owns no grid cells in a dim, then nlo > nhi
   special case: 2 procs share boundary which a grid point is exactly on
     2 equality if tests insure a consistent decision as to which proc owns it
------------------------------------------------------------------------- */

void Comm::partition_grid(int nx, int ny, int nz, double zfactor,
                          int &nxlo, int &nxhi, int &nylo, int &nyhi,
                          int &nzlo, int &nzhi)
{
  double xfraclo,xfrachi,yfraclo,yfrachi,zfraclo,zfrachi;

  if (layout != LAYOUT_TILED) {
    xfraclo = xsplit[myloc[0]];
    xfrachi = xsplit[myloc[0]+1];
    yfraclo = ysplit[myloc[1]];
    yfrachi = ysplit[myloc[1]+1];
    zfraclo = zsplit[myloc[2]];
    zfrachi = zsplit[myloc[2]+1];
  } else {
    xfraclo = mysplit[0][0];
    xfrachi = mysplit[0][1];
    yfraclo = mysplit[1][0];
    yfrachi = mysplit[1][1];
    zfraclo = mysplit[2][0];
    zfrachi = mysplit[2][1];
  }

  nxlo = static_cast<int> (xfraclo * nx);
  if (1.0*nxlo != xfraclo*nx) nxlo++;
  nxhi = static_cast<int> (xfrachi * nx);
  if (1.0*nxhi == xfrachi*nx) nxhi--;

  nylo = static_cast<int> (yfraclo * ny);
  if (1.0*nylo != yfraclo*ny) nylo++;
  nyhi = static_cast<int> (yfrachi * ny);
  if (1.0*nyhi == yfrachi*ny) nyhi--;

  if (zfactor == 0.0) {
    nzlo = static_cast<int> (zfraclo * nz);
    if (1.0*nzlo != zfraclo*nz) nzlo++;
    nzhi = static_cast<int> (zfrachi * nz);
    if (1.0*nzhi == zfrachi*nz) nzhi--;
  } else {
    nzlo = static_cast<int> (zfraclo * nz/zfactor);
    if (1.0*nzlo != zfraclo*nz) nzlo++;
    nzhi = static_cast<int> (zfrachi * nz/zfactor);
    if (1.0*nzhi == zfrachi*nz) nzhi--;
  }

  // OLD code
  // could sometimes map grid points slightly outside a proc to the proc

  /*
  if (layout != LAYOUT_TILED) {
    nxlo = static_cast<int> (xsplit[myloc[0]] * nx);
    nxhi = static_cast<int> (xsplit[myloc[0]+1] * nx) - 1;

    nylo = static_cast<int> (ysplit[myloc[1]] * ny);
    nyhi = static_cast<int> (ysplit[myloc[1]+1] * ny) - 1;

    if (zfactor == 0.0) {
      nzlo = static_cast<int> (zsplit[myloc[2]] * nz);
      nzhi = static_cast<int> (zsplit[myloc[2]+1] * nz) - 1;
    } else {
      nzlo = static_cast<int> (zsplit[myloc[2]] * nz/zfactor);
      nzhi = static_cast<int> (zsplit[myloc[2]+1] * nz/zfactor) - 1;
    }

  } else {
    nxlo = static_cast<int> (mysplit[0][0] * nx);
    nxhi = static_cast<int> (mysplit[0][1] * nx) - 1;

    nylo = static_cast<int> (mysplit[1][0] * ny);
    nyhi = static_cast<int> (mysplit[1][1] * ny) - 1;

    if (zfactor == 0.0) {
      nzlo = static_cast<int> (mysplit[2][0] * nz);
      nzhi = static_cast<int> (mysplit[2][1] * nz) - 1;
    } else {
      nzlo = static_cast<int> (mysplit[2][0] * nz/zfactor);
      nzhi = static_cast<int> (mysplit[2][1] * nz/zfactor) - 1;
    }
  }
  */
}

/* ----------------------------------------------------------------------
   communicate inbuf around full ring of processors with messtag
   nbytes = size of inbuf = n datums * nper bytes
   callback() is invoked to allow caller to process/update each proc's inbuf
   if self=1 (default), then callback() is invoked on final iteration
     using original inbuf, which may have been updated
   for non-nullptr outbuf, final updated inbuf is copied to it
     ok to specify outbuf = inbuf
   the ptr argument is a pointer to the instance of calling class
------------------------------------------------------------------------- */

void Comm::ring(int n, int nper, void *inbuf, int messtag,
                void (*callback)(int, char *, void *),
                void *outbuf, void *ptr, int self)
{
  MPI_Request request;
  MPI_Status status;

  int nbytes = n*nper;
  int maxbytes;
  MPI_Allreduce(&nbytes,&maxbytes,1,MPI_INT,MPI_MAX,world);

  // no need to communicate without data

  if (maxbytes == 0) return;

  // sanity check

  if ((nbytes > 0) && inbuf == nullptr)
    error->one(FLERR,"Cannot put data on ring from NULL pointer");

  char *buf,*bufcopy;
  memory->create(buf,maxbytes,"comm:buf");
  memory->create(bufcopy,maxbytes,"comm:bufcopy");
  if (nbytes && inbuf) memcpy(buf,inbuf,nbytes);

  int next = me + 1;
  int prev = me - 1;
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  for (int loop = 0; loop < nprocs; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxbytes,MPI_CHAR,prev,messtag,world,&request);
      MPI_Send(buf,nbytes,MPI_CHAR,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&nbytes);
      if (nbytes) memcpy(buf,bufcopy,nbytes);
    }
    if (self || loop < nprocs-1) callback(nbytes/nper,buf,ptr);
  }

  if (nbytes && outbuf) memcpy(outbuf,buf,nbytes);

  memory->destroy(buf);
  memory->destroy(bufcopy);
}

/* ----------------------------------------------------------------------
   rendezvous communication operation
   three stages:
     first comm sends inbuf from caller decomp to rvous decomp
     callback operates on data in rendezvous decomp
     second comm sends outbuf from rvous decomp back to caller decomp
   inputs:
     which = perform (0) irregular or (1) MPI_All2allv communication
     n = # of datums in inbuf
     inbuf = vector of input datums
     insize = byte size of each input datum
     inorder = 0 for inbuf in random proc order, 1 for datums ordered by proc
     procs: inorder 0 = proc to send each datum to, 1 = # of datums/proc,
     callback = caller function to invoke in rendezvous decomposition
                takes input datums, returns output datums
     outorder = same as inorder, but for datums returned by callback()
     ptr = pointer to caller class, passed to callback()
   outputs:
     nout = # of output datums (function return)
     outbuf = vector of output datums
     outsize = byte size of each output datum
   callback inputs:
     nrvous = # of rvous decomp datums in inbuf_rvous
     inbuf_rvous = vector of rvous decomp input datums
     ptr = pointer to caller class
   callback outputs:
     nrvous_out = # of rvous decomp output datums (function return)
     flag = 0 for no second comm, 1 for outbuf_rvous = inbuf_rvous,
            2 for second comm with new outbuf_rvous
     procs_rvous = outorder 0 = proc to send each datum to, 1 = # of datums/proc
                   allocated
     outbuf_rvous = vector of rvous decomp output datums
   NOTE: could use MPI_INT or MPI_DOUBLE insead of MPI_CHAR
         to avoid checked-for overflow in MPI_Alltoallv?
------------------------------------------------------------------------- */

int Comm::
rendezvous(int which, int n, char *inbuf, int insize,
           int inorder, int *procs,
           int (*callback)(int, char *, int &, int *&, char *&, void *),
           int outorder, char *&outbuf, int outsize, void *ptr, int statflag)
{
  if (which == 0)
    return rendezvous_irregular(n,inbuf,insize,inorder,procs,callback,
                                outorder,outbuf,outsize,ptr,statflag);
  else
    return rendezvous_all2all(n,inbuf,insize,inorder,procs,callback,
                              outorder,outbuf,outsize,ptr,statflag);
}

/* ---------------------------------------------------------------------- */

int Comm::
rendezvous_irregular(int n, char *inbuf, int insize, int inorder, int *procs,
                     int (*callback)(int, char *, int &, int *&, char *&, void *),
                     int outorder, char *&outbuf,
                     int outsize, void *ptr, int statflag)
{
  // irregular comm of inbuf from caller decomp to rendezvous decomp

  auto irregular = new Irregular(lmp);

  int nrvous;
  if (inorder) nrvous = irregular->create_data_grouped(n,procs);
  else nrvous = irregular->create_data(n,procs);

  // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

  auto inbuf_rvous = (char *) memory->smalloc((bigint) nrvous*insize+1, "rendezvous:inbuf");
  irregular->exchange_data(inbuf,insize,inbuf_rvous);

  bigint irregular1_bytes = irregular->memory_usage();
  irregular->destroy_data();
  delete irregular;

  // peform rendezvous computation via callback()
  // callback() allocates/populates proclist_rvous and outbuf_rvous

  int flag;
  int *procs_rvous;
  char *outbuf_rvous;
  int nrvous_out = callback(nrvous,inbuf_rvous,flag, procs_rvous,outbuf_rvous,ptr);

  if (flag != 1) memory->sfree(inbuf_rvous);  // outbuf_rvous = inbuf_vous
  if (flag == 0) {
    if (statflag) rendezvous_stats(n,0,nrvous,nrvous_out,insize,outsize,
                                   (bigint) nrvous_out*sizeof(int) + irregular1_bytes);
    return 0;    // all nout_rvous are 0, no 2nd comm stage
  }

  // irregular comm of outbuf from rendezvous decomp back to caller decomp
  // caller will free outbuf

  irregular = new Irregular(lmp);

  int nout;
  if (outorder) nout = irregular->create_data_grouped(nrvous_out,procs_rvous);
  else nout = irregular->create_data(nrvous_out,procs_rvous);

  // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

  outbuf = (char *) memory->smalloc((bigint) nout*outsize+1, "rendezvous:outbuf");
  irregular->exchange_data(outbuf_rvous,outsize,outbuf);

  bigint irregular2_bytes = irregular->memory_usage();
  irregular->destroy_data();
  delete irregular;

  memory->destroy(procs_rvous);
  memory->sfree(outbuf_rvous);

  // return number of output datums
  // last arg to stats() = memory for procs_rvous + irregular comm

  if (statflag) rendezvous_stats(n,nout,nrvous,nrvous_out,insize,outsize,
                                 (bigint) nrvous_out*sizeof(int) +
                                 MAX(irregular1_bytes,irregular2_bytes));
  return nout;
}

/* ---------------------------------------------------------------------- */

int Comm::
rendezvous_all2all(int n, char *inbuf, int insize, int inorder, int *procs,
                   int (*callback)(int, char *, int &, int *&, char *&, void *),
                   int outorder, char *&outbuf, int outsize, void *ptr,
                   int statflag)
{
  int iproc;
  bigint all2all1_bytes,all2all2_bytes;
  int *sendcount,*sdispls,*recvcount,*rdispls;
  int *procs_a2a;
  bigint *offsets;
  char *inbuf_a2a,*outbuf_a2a;

  // create procs and inbuf for All2all if necessary

  if (!inorder) {
    memory->create(procs_a2a,nprocs,"rendezvous:procs");

    // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

    inbuf_a2a = (char *) memory->smalloc((bigint) n*insize+1,
                                         "rendezvous:inbuf");
    memset(inbuf_a2a,0,(bigint)n*insize*sizeof(char));
    memory->create(offsets,nprocs,"rendezvous:offsets");

    for (int i = 0; i < nprocs; i++) procs_a2a[i] = 0;
    for (int i = 0; i < n; i++) procs_a2a[procs[i]]++;

    offsets[0] = 0;
    for (int i = 1; i < nprocs; i++)
      offsets[i] = offsets[i-1] + (bigint)insize*procs_a2a[i-1];

    bigint offset = 0;
    for (int i = 0; i < n; i++) {
      iproc = procs[i];
      memcpy(&inbuf_a2a[offsets[iproc]],&inbuf[offset],insize);
      offsets[iproc] += insize;
      offset += insize;
    }

    all2all1_bytes = nprocs*sizeof(int) + nprocs*sizeof(bigint)
                     + (bigint)n*insize;

  } else {
    procs_a2a = procs;
    inbuf_a2a = inbuf;
    all2all1_bytes = 0;
  }

  // create args for MPI_Alltoallv() on input data

  memory->create(sendcount,nprocs,"rendezvous:sendcount");
  memcpy(sendcount,procs_a2a,nprocs*sizeof(int));

  memory->create(recvcount,nprocs,"rendezvous:recvcount");
  MPI_Alltoall(sendcount,1,MPI_INT,recvcount,1,MPI_INT,world);

  memory->create(sdispls,nprocs,"rendezvous:sdispls");
  memory->create(rdispls,nprocs,"rendezvous:rdispls");
  sdispls[0] = rdispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    sdispls[i] = sdispls[i-1] + sendcount[i-1];
    rdispls[i] = rdispls[i-1] + recvcount[i-1];
  }
  int nrvous = rdispls[nprocs-1] + recvcount[nprocs-1];

  // test for overflow of input data due to imbalance or insize
  // means that individual sdispls or rdispls values overflow

  int overflow = 0;
  if ((bigint) n*insize > MAXSMALLINT) overflow = 1;
  if ((bigint) nrvous*insize > MAXSMALLINT) overflow = 1;
  int overflowall;
  MPI_Allreduce(&overflow,&overflowall,1,MPI_INT,MPI_MAX,world);
  if (overflowall) error->all(FLERR,"Overflow input size in rendezvous_a2a");

  for (int i = 0; i < nprocs; i++) {
    sendcount[i] *= insize;
    sdispls[i] *= insize;
    recvcount[i] *= insize;
    rdispls[i] *= insize;
  }

  // all2all comm of inbuf from caller decomp to rendezvous decomp
  // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

  auto inbuf_rvous = (char *) memory->smalloc((bigint) nrvous*insize+1, "rendezvous:inbuf");
  memset(inbuf_rvous,0,(bigint) nrvous*insize*sizeof(char));

  MPI_Alltoallv(inbuf_a2a,sendcount,sdispls,MPI_CHAR,
                inbuf_rvous,recvcount,rdispls,MPI_CHAR,world);

  if (!inorder) {
    memory->destroy(procs_a2a);
    memory->sfree(inbuf_a2a);
    memory->destroy(offsets);
  }

  // peform rendezvous computation via callback()
  // callback() allocates/populates proclist_rvous and outbuf_rvous

  int flag;
  int *procs_rvous;
  char *outbuf_rvous;

  int nrvous_out = callback(nrvous,inbuf_rvous,flag, procs_rvous,outbuf_rvous,ptr);

  if (flag != 1) memory->sfree(inbuf_rvous);  // outbuf_rvous = inbuf_vous
  if (flag == 0) {
    memory->destroy(sendcount);
    memory->destroy(recvcount);
    memory->destroy(sdispls);
    memory->destroy(rdispls);
    if (statflag) rendezvous_stats(n,0,nrvous,nrvous_out,insize,outsize,
                                   (bigint) nrvous_out*sizeof(int) +
                                   4*nprocs*sizeof(int) + all2all1_bytes);
    return 0;    // all nout_rvous are 0, no 2nd irregular
  }

  // create procs and outbuf for All2all if necessary

  if (!outorder) {
    memory->create(procs_a2a,nprocs,"rendezvous_a2a:procs");

    // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

    outbuf_a2a = (char *) memory->smalloc((bigint) nrvous_out*outsize+1, "rendezvous:outbuf");
    memory->create(offsets,nprocs,"rendezvous:offsets");

    for (int i = 0; i < nprocs; i++) procs_a2a[i] = 0;
    for (int i = 0; i < nrvous_out; i++) procs_a2a[procs_rvous[i]]++;

    offsets[0] = 0;
    for (int i = 1; i < nprocs; i++)
      offsets[i] = offsets[i-1] + (bigint)outsize*procs_a2a[i-1];

    bigint offset = 0;
    for (int i = 0; i < nrvous_out; i++) {
      iproc = procs_rvous[i];
      memcpy(&outbuf_a2a[offsets[iproc]],&outbuf_rvous[offset],outsize);
      offsets[iproc] += outsize;
      offset += outsize;
    }

    all2all2_bytes = nprocs*sizeof(int) + nprocs*sizeof(bigint) + (bigint)nrvous_out*outsize;

  } else {
    procs_a2a = procs_rvous;
    outbuf_a2a = outbuf_rvous;
    all2all2_bytes = 0;
  }

  // comm outbuf from rendezvous decomposition back to caller

  memcpy(sendcount,procs_a2a,nprocs*sizeof(int));

  MPI_Alltoall(sendcount,1,MPI_INT,recvcount,1,MPI_INT,world);

  sdispls[0] = rdispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    sdispls[i] = sdispls[i-1] + sendcount[i-1];
    rdispls[i] = rdispls[i-1] + recvcount[i-1];
  }
  int nout = rdispls[nprocs-1] + recvcount[nprocs-1];

  // test for overflow of outbuf due to imbalance or outsize
  // means that individual sdispls or rdispls values overflow

  overflow = 0;
  if ((bigint) nrvous*outsize > MAXSMALLINT) overflow = 1;
  if ((bigint) nout*outsize > MAXSMALLINT) overflow = 1;
  MPI_Allreduce(&overflow,&overflowall,1,MPI_INT,MPI_MAX,world);
  if (overflowall) error->all(FLERR,"Overflow output in rendezvous_a2a");

  for (int i = 0; i < nprocs; i++) {
    sendcount[i] *= outsize;
    sdispls[i] *= outsize;
    recvcount[i] *= outsize;
    rdispls[i] *= outsize;
  }

  // all2all comm of outbuf from rendezvous decomp back to caller decomp
  // caller will free outbuf
  // add 1 item to the allocated buffer size, so the returned pointer is not a null pointer

  outbuf = (char *) memory->smalloc((bigint) nout*outsize+1,"rendezvous:outbuf");

  MPI_Alltoallv(outbuf_a2a,sendcount,sdispls,MPI_CHAR,
                outbuf,recvcount,rdispls,MPI_CHAR,world);

  memory->destroy(procs_rvous);
  memory->sfree(outbuf_rvous);

  if (!outorder) {
    memory->destroy(procs_a2a);
    memory->sfree(outbuf_a2a);
    memory->destroy(offsets);
  }

  // clean up

  memory->destroy(sendcount);
  memory->destroy(recvcount);
  memory->destroy(sdispls);
  memory->destroy(rdispls);

  // return number of output datums
  // last arg to stats() = mem for procs_rvous + per-proc vecs + reordering ops

  if (statflag) rendezvous_stats(n,nout,nrvous,nrvous_out,insize,outsize,
                                 (bigint) nrvous_out*sizeof(int) +
                                 4*nprocs*sizeof(int) +
                                 MAX(all2all1_bytes,all2all2_bytes));
  return nout;
}

/* ----------------------------------------------------------------------
   print balance and memory info for rendezvous operation
   useful for debugging
------------------------------------------------------------------------- */

void Comm::rendezvous_stats(int n, int nout, int nrvous, int nrvous_out,
                            int insize, int outsize, bigint commsize)
{
  bigint size_in_all,size_in_max,size_in_min;
  bigint size_out_all,size_out_max,size_out_min;
  bigint size_inrvous_all,size_inrvous_max,size_inrvous_min;
  bigint size_outrvous_all,size_outrvous_max,size_outrvous_min;
  bigint size_comm_all,size_comm_max,size_comm_min;

  bigint size = (bigint) n*insize;
  MPI_Allreduce(&size,&size_in_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_in_max,1,MPI_LMP_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_in_min,1,MPI_LMP_BIGINT,MPI_MIN,world);

  size = (bigint) nout*outsize;
  MPI_Allreduce(&size,&size_out_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_out_max,1,MPI_LMP_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_out_min,1,MPI_LMP_BIGINT,MPI_MIN,world);

  size = (bigint) nrvous*insize;
  MPI_Allreduce(&size,&size_inrvous_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_inrvous_max,1,MPI_LMP_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_inrvous_min,1,MPI_LMP_BIGINT,MPI_MIN,world);

  size = (bigint) nrvous_out*insize;
  MPI_Allreduce(&size,&size_outrvous_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_outrvous_max,1,MPI_LMP_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_outrvous_min,1,MPI_LMP_BIGINT,MPI_MIN,world);

  size = commsize;
  MPI_Allreduce(&size,&size_comm_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_comm_max,1,MPI_LMP_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_comm_min,1,MPI_LMP_BIGINT,MPI_MIN,world);

  int mbytes = 1024*1024;

  if (me == 0) {
    std::string mesg = "Rendezvous balance and memory info: (tot,ave,max,min) \n";
    mesg += fmt::format("  input datum count: {} {} {} {}\n",
                        size_in_all/insize,1.0*size_in_all/nprocs/insize,
                        size_in_max/insize,size_in_min/insize);
    mesg += fmt::format("  input data (MB): {:.6} {:.6} {:.6} {:.6}\n",
                        1.0*size_in_all/mbytes,1.0*size_in_all/nprocs/mbytes,
                        1.0*size_in_max/mbytes,1.0*size_in_min/mbytes);
    if (outsize)
      mesg += fmt::format("  output datum count: {} {} {} {}\n",
                          size_out_all/outsize,1.0*size_out_all/nprocs/outsize,
                          size_out_max/outsize,size_out_min/outsize);
    else
      mesg += fmt::format("  output datum count: {} {:.6} {} {}\n",0,0.0,0,0);

    mesg += fmt::format("  output data (MB): {:.6} {:.6} {:.6} {:.6}\n",
                        1.0*size_out_all/mbytes,1.0*size_out_all/nprocs/mbytes,
                        1.0*size_out_max/mbytes,1.0*size_out_min/mbytes);
    mesg += fmt::format("  input rvous datum count: {} {} {} {}\n",
                        size_inrvous_all/insize,1.0*size_inrvous_all/nprocs/insize,
                        size_inrvous_max/insize,size_inrvous_min/insize);
    mesg += fmt::format("  input rvous data (MB): {:.6} {:.6} {:.6} {:.6}\n",
                        1.0*size_inrvous_all/mbytes,1.0*size_inrvous_all/nprocs/mbytes,
                        1.0*size_inrvous_max/mbytes,1.0*size_inrvous_min/mbytes);
    if (outsize)
      mesg += fmt::format("  output rvous datum count: {} {} {} {}\n",
                          size_outrvous_all/outsize,1.0*size_outrvous_all/nprocs/outsize,
                          size_outrvous_max/outsize,size_outrvous_min/outsize);
    else
      mesg += fmt::format("  output rvous datum count: {} {:.6} {} {}\n",0,0.0,0,0);
    mesg += fmt::format("  output rvous data (MB): {:.6} {:.6} {:.6} {:.6}\n",
                        1.0*size_outrvous_all/mbytes,1.0*size_outrvous_all/nprocs/mbytes,
                        1.0*size_outrvous_max/mbytes,1.0*size_outrvous_min/mbytes);
    mesg += fmt::format("  rvous comm (MB): {:.6} {:.6} {:.6} {:.6}\n",
                        1.0*size_comm_all/mbytes,1.0*size_comm_all/nprocs/mbytes,
                        1.0*size_comm_max/mbytes,1.0*size_comm_min/mbytes);
    utils::logmesg(lmp,mesg);
  }
}
