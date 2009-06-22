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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_evaporate.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "group.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixEvaporate::FixEvaporate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all("Illegal fix evaporate command");

  scalar_flag = 1;
  scalar_vector_freq = 1;
  extscalar = 0;

  nevery = atoi(arg[3]);
  nflux = atoi(arg[4]);
  iregion = domain->find_region(arg[5]);
  int seed = atoi(arg[6]);

  if (nevery <= 0 || nflux <= 0) error->all("Illegal fix evaporate command");
  if (iregion == -1) error->all("Fix evaporate region ID does not exist");
  if (seed <= 0) error->all("Illegal fix evaporate command");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
  ndeleted = 0;

  nmax = 0;
  list = NULL;
  mark = NULL;
}

/* ---------------------------------------------------------------------- */

FixEvaporate::~FixEvaporate()
{
  delete random;
  memory->sfree(list);
  memory->sfree(mark);
}

/* ---------------------------------------------------------------------- */

int FixEvaporate::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEvaporate::init()
{
  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all("Cannot evaporate atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   perform particle deletion
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixEvaporate::pre_exchange()
{
  int i,iwhichglobal,iwhichlocal;

  if (update->ntimestep != next_reneighbor) return;

  // grow list and mark arrays if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(list);
    memory->sfree(mark);
    nmax = atom->nmax;
    list = (int *) memory->smalloc(nmax*sizeof(int),"evaporate:list");
    mark = (int *) memory->smalloc(nmax*sizeof(int),"evaporate:mark");
  }

  // nall = # of deletable atoms in region
  // nbefore = # on procs before me

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ncount = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	list[ncount++] = i;

  int nall,nbefore;
  MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,world);
  nbefore -= ncount;
  
  // nwhack = total number of atoms to delete
  // choose atoms randomly across all procs and mark them for deletion
  // shrink local list of candidates as my atoms get marked for deletion

  int nwhack = MIN(nflux,nall);

  for (i = 0; i < nlocal; i++) mark[i] = 0;
  
  for (i = 0; i < nwhack; i++) {
    iwhichglobal = static_cast<int> (nall*random->uniform());
    if (iwhichglobal < nbefore) nbefore--;
    else if (iwhichglobal < nbefore + ncount) {
      iwhichlocal = iwhichglobal - nbefore;
      mark[list[iwhichlocal]] = 1;
      list[iwhichlocal] = list[ncount-1];
      ncount--;
    }
    nall--;
  }
  
  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms
  
  AtomVec *avec = atom->avec;
  
  for (i = nlocal-1; i >= 0; i--) {
    if (mark[i]) {
      avec->copy(atom->nlocal-1,i);
      atom->nlocal--;
    }
  }

  // reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  atom->natoms -= nwhack;
  if (nwhack && atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // statistics
  
  ndeleted += nwhack;
  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
   return number of deleted particles
------------------------------------------------------------------------- */

double FixEvaporate::compute_scalar()
{
  return 1.0*ndeleted;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixEvaporate::memory_usage()
{
  double bytes = 2*nmax * sizeof(int);
  return bytes;
}
