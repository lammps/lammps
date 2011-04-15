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
  if (narg < 7) error->all("Illegal fix evaporate command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  nevery = atoi(arg[3]);
  nflux = atoi(arg[4]);
  iregion = domain->find_region(arg[5]);
  int n = strlen(arg[5]) + 1;
  idregion = new char[n];
  strcpy(idregion,arg[5]);
  int seed = atoi(arg[6]);

  if (nevery <= 0 || nflux <= 0) error->all("Illegal fix evaporate command");
  if (iregion == -1) error->all("Region ID for fix evaporate does not exist");
  if (seed <= 0) error->all("Illegal fix evaporate command");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // optional args

  molflag = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix evaporate command");
      if (strcmp(arg[iarg+1],"no") == 0) molflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) molflag = 1;
      else error->all("Illegal fix evaporate command");
      iarg += 2;
    } else error->all("Illegal fix evaporate command");
  }

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
  delete [] idregion;
  delete random;
  memory->destroy(list);
  memory->destroy(mark);
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
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all("Region ID for fix evaporate does not exist");

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

  // if molflag not set, warn if any deletable atom has a mol ID

  if (molflag == 0 && atom->molecule_flag) {
    int *molecule = atom->molecule;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	if (molecule[i]) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->
	warning("Fix evaporate may delete atom with non-zero molecule ID");
  }

  if (molflag && atom->molecule_flag == 0)
      error->all("Fix evaporate molecule requires atom attribute molecule");
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
    memory->destroy(list);
    memory->destroy(mark);
    nmax = atom->nmax;
    memory->create(list,nmax,"evaporate:list");
    memory->create(mark,nmax,"evaporate:mark");
  }

  // ncount = # of deletable atoms in region that I own
  // nall = # on all procs
  // nbefore = # on procs before me
  // list[ncount] = list of local indices

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

  // ndel = total # of atom deletions, in or out of region
  // mark[] = 1 if deleted

  int ndel = 0;
  for (i = 0; i < nlocal; i++) mark[i] = 0;

  // atomic deletions
  // choose atoms randomly across all procs and mark them for deletion
  // shrink eligible list as my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each atom deletion

  if (molflag == 0) {
    while (nall && ndel < nflux) {
      iwhichglobal = static_cast<int> (nall*random->uniform());
      if (iwhichglobal < nbefore) nbefore--;
      else if (iwhichglobal < nbefore + ncount) {
	iwhichlocal = iwhichglobal - nbefore;
	mark[list[iwhichlocal]] = 1;
	list[iwhichlocal] = list[ncount-1];
	ncount--;
      }
      ndel++;
      nall--;
    }

  // molecule deletions
  // choose one atom in one molecule randomly across all procs
  // bcast mol ID and delete all atoms in that molecule on any proc
  // update deletion count by total # of atoms in molecule
  // shrink list of eligible candidates as any of my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each molecule deletion

  } else {
    int me,proc,iatom,imolecule,ndelone,ndelall;
    int *molecule = atom->molecule;

    while (nall && ndel < nflux) {

      // pick an iatom,imolecule on proc me to delete

      iwhichglobal = static_cast<int> (nall*random->uniform());
      if (iwhichglobal >= nbefore && iwhichglobal < nbefore + ncount) {
	iwhichlocal = iwhichglobal - nbefore;
	iatom = list[iwhichlocal];
	imolecule = molecule[iatom];
	me = comm->me;
      } else me = -1;

      // bcast mol ID to delete all atoms from
      // delete single iatom if mol ID = 0

      MPI_Allreduce(&me,&proc,1,MPI_INT,MPI_MAX,world);
      MPI_Bcast(&imolecule,1,MPI_INT,proc,world);
      ndelone = 0;
      for (i = 0; i < nlocal; i++) {
	if (imolecule && molecule[i] == imolecule) {
	  mark[i] = 1;
	  ndelone++;
	} else if (me == proc && i == iatom) {
	  mark[i] = 1;
	  ndelone++;
	}
      }

      // remove any atoms marked for deletion from my eligible list

      i = 0;
      while (i < ncount) {
	if (mark[list[i]]) {
	  list[i] = list[ncount-1];
	  ncount--;
	} else i++;
      }

      // update ndel,ncount,nall,nbefore

      MPI_Allreduce(&ndelone,&ndelall,1,MPI_INT,MPI_SUM,world);
      ndel += ndelall;
      MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,world);
      MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,world);
      nbefore -= ncount;
    }
  }

  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms
  
  AtomVec *avec = atom->avec;
  
  for (i = nlocal-1; i >= 0; i--) {
    if (mark[i]) {
      avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
    }
  }

  // reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  atom->natoms -= ndel;
  if (ndel && atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // statistics
  
  ndeleted += ndel;
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
