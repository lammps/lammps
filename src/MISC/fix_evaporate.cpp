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
#include "force.h"
#include "group.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEvaporate::FixEvaporate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix evaporate command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  nevery = force->inumeric(FLERR,arg[3]);
  nflux = force->inumeric(FLERR,arg[4]);
  iregion = domain->find_region(arg[5]);
  int n = strlen(arg[5]) + 1;
  idregion = new char[n];
  strcpy(idregion,arg[5]);
  int seed = force->inumeric(FLERR,arg[6]);

  if (nevery <= 0 || nflux <= 0)
    error->all(FLERR,"Illegal fix evaporate command");
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix evaporate does not exist");
  if (seed <= 0) error->all(FLERR,"Illegal fix evaporate command");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // optional args

  molflag = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix evaporate command");
      if (strcmp(arg[iarg+1],"no") == 0) molflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) molflag = 1;
      else error->all(FLERR,"Illegal fix evaporate command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix evaporate command");
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
    error->all(FLERR,"Region ID for fix evaporate does not exist");

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
      error->all(FLERR,"Cannot evaporate atoms in atom_modify first group");
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
      error->warning(FLERR,
                     "Fix evaporate may delete atom with non-zero molecule ID");
  }

  if (molflag && atom->molecule_flag == 0)
      error->all(FLERR,
                 "Fix evaporate molecule requires atom attribute molecule");
}

/* ----------------------------------------------------------------------
   perform particle deletion
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixEvaporate::pre_exchange()
{
  int i,j,m,iwhichglobal,iwhichlocal;
  int ndel,ndeltopo[4];

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
  // list[ncount] = list of local indices of atoms I can delete

  double **x = atom->x;
  int *mask = atom->mask;
  int *tag = atom->tag;
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
  // ndeltopo[1,2,3,4] = ditto for bonds, angles, dihedrals, impropers
  // mark[] = 1 if deleted

  ndel = 0;
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
  // keep ndel,ndeltopo,ncount,nall,nbefore current after each mol deletion

  } else {
    int me,proc,iatom,imolecule,ndelone,ndelall;
    int *molecule = atom->molecule;

    ndeltopo[0] = ndeltopo[1] = ndeltopo[2] = ndeltopo[3] = 0;

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
      // if mol ID > 0, delete any atom in molecule and decrement counters
      // if mol ID == 0, delete single iatom
      // be careful to delete correct # of bond,angle,etc for newton on or off

      MPI_Allreduce(&me,&proc,1,MPI_INT,MPI_MAX,world);
      MPI_Bcast(&imolecule,1,MPI_INT,proc,world);
      ndelone = 0;
      for (i = 0; i < nlocal; i++) {
        if (imolecule && molecule[i] == imolecule) {
          mark[i] = 1;
          ndelone++;

          if (atom->avec->bonds_allow) {
            if (force->newton_bond) ndeltopo[0] += atom->num_bond[i];
            else {
              for (j = 0; j < atom->num_bond[i]; j++) {
                if (tag[i] < atom->bond_atom[i][j]) ndeltopo[0]++;
              }
            }
          }
          if (atom->avec->angles_allow) {
            if (force->newton_bond) ndeltopo[1] += atom->num_angle[i];
            else {
              for (j = 0; j < atom->num_angle[i]; j++) {
                m = atom->map(atom->angle_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[1]++;
              }
            }
          }
          if (atom->avec->dihedrals_allow) {
            if (force->newton_bond) ndeltopo[2] += atom->num_dihedral[i];
            else {
              for (j = 0; j < atom->num_dihedral[i]; j++) {
                m = atom->map(atom->dihedral_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[2]++;
              }
            }
          }
          if (atom->avec->impropers_allow) {
            if (force->newton_bond) ndeltopo[3] += atom->num_improper[i];
            else {
              for (j = 0; j < atom->num_improper[i]; j++) {
                m = atom->map(atom->improper_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[3]++;
              }
            }
          }

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
      // ndelall is total atoms deleted on this iteration
      // ncount is already correct, so resum to get nall and nbefore

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

  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  atom->natoms -= ndel;
  if (molflag) {
    int all[4];
    MPI_Allreduce(ndeltopo,all,4,MPI_INT,MPI_SUM,world);
    atom->nbonds -= all[0];
    atom->nangles -= all[1];
    atom->ndihedrals -= all[2];
    atom->nimpropers -= all[3];
  }

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
