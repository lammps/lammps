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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "universe.h"
#include "version.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
------------------------------------------------------------------------- */

Universe::Universe(LAMMPS *lmp, MPI_Comm communicator) : Pointers(lmp)
{
  version = (char *) LAMMPS_VERSION;

  uworld = original = communicator;
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  uscreen = stdout;
  ulogfile = NULL;

  existflag = 0;
  nworlds = 0;
  procs_per_world = NULL;
  root_proc = NULL;

  memory->create(proc2original,nprocs,"universe:proc2original");
  for (int i = 0; i < nprocs; i++) proc2original[i] = i;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  if (uworld != original) MPI_Comm_free(&uworld);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(proc2original);
}

/* ----------------------------------------------------------------------
   placeholder routine, not yet operational
   permute the mapping of universe procs in uworld to procs in original
------------------------------------------------------------------------- */

void Universe::reorder(int key)
{
  if (uworld != original) MPI_Comm_free(&uworld);

  MPI_Comm_split(original,0,key,&uworld);
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  int ome;
  MPI_Comm_rank(original,&ome);
  MPI_Allgather(&ome,1,MPI_INT,proc2original,1,MPI_INT,uworld);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == NULL -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
------------------------------------------------------------------------- */

void Universe::add_world(char *str)
{
  int n,nper;
  char *ptr;

  if (str == NULL) {
    n = 1;
    nper = nprocs;
  } else if ((ptr = strchr(str,'x')) != NULL) {
    *ptr = '\0';
    n = atoi(str);
    nper = atoi(ptr+1);
  } else {
    n = 1;
    nper = atoi(str);
  }

  memory->grow(procs_per_world,nworlds+n,"universe:procs_per_world");
  memory->grow(root_proc,(nworlds+n),"universe:root_proc");

  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0) root_proc[nworlds] = 0;
    else
      root_proc[nworlds] = root_proc[nworlds-1] + procs_per_world[nworlds-1];
    if (me >= root_proc[nworlds]) iworld = nworlds;
    nworlds++;
  }
}

/* ----------------------------------------------------------------------
   check if total procs in all worlds = procs in universe
------------------------------------------------------------------------- */

int Universe::consistent()
{
  int n = 0;
  for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
  if (n == nprocs) return 1;
  else return 0;
}
