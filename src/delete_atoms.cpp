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

#include "stdlib.h"
#include "string.h"
#include "delete_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "region.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace LAMMPS_NS;

// allocate space for static class variable

DeleteAtoms *DeleteAtoms::cptr;

/* ---------------------------------------------------------------------- */

DeleteAtoms::DeleteAtoms(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DeleteAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Delete_atoms command before simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal delete_atoms command");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use delete_atoms unless atoms have IDs");

  // store state before delete

  bigint natoms_previous = atom->natoms;

  // delete the atoms

  if (strcmp(arg[0],"group") == 0) delete_group(narg,arg);
  else if (strcmp(arg[0],"region") == 0) delete_region(narg,arg);
  else if (strcmp(arg[0],"overlap") == 0) delete_overlap(narg,arg);
  else if (strcmp(arg[0],"porosity") == 0) delete_porosity(narg,arg);
  else error->all(FLERR,"Illegal delete_atoms command");

  // delete local atoms flagged in dlist
  // reset nlocal

  AtomVec *avec = atom->avec;
  int nlocal = atom->nlocal;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  // if non-molecular system and compress flag set,
  // reset atom tags to be contiguous
  // set all atom IDs to 0, call tag_extend()

  if (atom->molecular == 0 && compress_flag) {
    int *tag = atom->tag;
    for (i = 0; i < nlocal; i++) tag[i] = 0;
    atom->tag_extend();
  }

  // reset atom->natoms
  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // print before and after atom count

  bigint ndelete = natoms_previous - atom->natoms;

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Deleted " BIGINT_FORMAT
                        " atoms, new total = " BIGINT_FORMAT "\n",
                        ndelete,atom->natoms);
    if (logfile) fprintf(logfile,"Deleted " BIGINT_FORMAT
                         " atoms, new total = " BIGINT_FORMAT "\n",
                         ndelete,atom->natoms);
  }
}

/* ----------------------------------------------------------------------
   delete all atoms in group
   group will still exist
------------------------------------------------------------------------- */

void DeleteAtoms::delete_group(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal delete_atoms command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find delete_atoms group ID");
  options(narg-2,&arg[2]);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete all atoms in region
   if mol_flag is set, also delete atoms in molecules with any deletions
------------------------------------------------------------------------- */

void DeleteAtoms::delete_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal delete_atoms command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Could not find delete_atoms region ID");
  options(narg-2,&arg[2]);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++)
    if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) dlist[i] = 1;

  if (mol_flag == 0) return;

  // delete entire molecules if any atom in molecule was deleted
  // store list of molecule IDs I delete atoms from in list
  // pass list from proc to proc via ring communication

  hash = new std::map<int,int>();

  int *molecule = atom->molecule;
  for (int i = 0; i < nlocal; i++)
    if (dlist[i] && hash->find(molecule[i]) == hash->end())
      (*hash)[molecule[i]] = 1;

  int n = hash->size();
  int *list;
  memory->create(list,n,"delete_atoms:list");

  n = 0;
  std::map<int,int>::iterator pos;
  for (pos = hash->begin(); pos != hash->end(); ++pos) list[n++] = pos->first;

  cptr = this;
  comm->ring(n,sizeof(int),list,1,molring,NULL);

  delete hash;
  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   callback from comm->ring()
   cbuf = list of N molecule IDs, put them in hash
   loop over my atoms, if matches moleculed ID in hash, delete that atom
------------------------------------------------------------------------- */

void DeleteAtoms::molring(int n, char *cbuf)
{
  int *list = (int *) cbuf;
  int *dlist = cptr->dlist;
  std::map<int,int> *hash = cptr->hash;
  int nlocal = cptr->atom->nlocal;
  int *molecule = cptr->atom->molecule;

  hash->clear();
  for (int i = 0; i < n; i++) (*hash)[list[i]] = 1;

  for (int i = 0; i < nlocal; i++)
    if (hash->find(molecule[i]) != hash->end()) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete atoms so there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

void DeleteAtoms::delete_overlap(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal delete_atoms command");

  // read args

  double cut = force->numeric(FLERR,arg[1]);
  double cutsq = cut*cut;

  int igroup1 = group->find(arg[2]);
  int igroup2 = group->find(arg[3]);
  if (igroup1 < 0 || igroup2 < 0)
    error->all(FLERR,"Could not find delete_atoms group ID");
  options(narg-4,&arg[4]);

  int group1bit = group->bitmask[igroup1];
  int group2bit = group->bitmask[igroup2];

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for delete_atoms ...\n");

  // request a full neighbor list for use by this command

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  lmp->init();

  // error check on cutoff
  // if no pair style, neighbor list will be empty

  if (force->pair == NULL)
    error->all(FLERR,"Delete_atoms requires a pair style be defined");
  if (cut > neighbor->cutneighmax)
    error->all(FLERR,"Delete_atoms cutoff > neighbor cutoff");

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor list based on earlier request

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  NeighList *list = neighbor->lists[irequest];
  neighbor->build_one(irequest);

  // allocate and initialize deletion list
  // must be after exchange potentially changes nlocal

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // double loop over owned atoms and their full neighbor list
  // at end of loop, there are no more overlaps
  // only ever delete owned atom I, never J even if owned

  int *tag = atom->tag;
  int *mask = atom->mask;
  double **x = atom->x;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      // only consider deletion if I,J distance < cutoff

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;

      // only consider deletion if I,J are in groups 1,2 respectively
      // true whether J is owned or ghost atom

      if (!(mask[i] & group1bit)) continue;
      if (!(mask[j] & group2bit)) continue;

      // J is owned atom:
      //   delete atom I if atom J has not already been deleted
      // J is ghost atom:
      //   delete atom I if J,I is not a candidate deletion pair
      //     due to being in groups 1,2 respectively
      //   if they are candidate pair, then either:
      //      another proc owns J and could delete J
      //      J is a ghost of another of my owned atoms, and I could delete J
      //   test on tags of I,J insures that only I or J is deleted

      if (j < nlocal) {
        if (dlist[j]) continue;
      } else if ((mask[i] & group2bit) && (mask[j] & group1bit)) {
        if (tag[i] > tag[j]) continue;
      }

      dlist[i] = 1;
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   create porosity by deleting atoms in a specified region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_porosity(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal delete_atoms command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Could not find delete_atoms region ID");

  double porosity_fraction = force->numeric(FLERR,arg[2]);
  int seed = force->inumeric(FLERR,arg[3]);
  options(narg-4,&arg[4]);

  RanMars *random = new RanMars(lmp,seed + comm->me);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++)
    if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
      if (random->uniform() <= porosity_fraction) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   process command options
------------------------------------------------------------------------- */

void DeleteAtoms::options(int narg, char **arg)
{
  compress_flag = 1;
  mol_flag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal delete_atoms command");
      if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
      else error->all(FLERR,"Illegal delete_atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal delete_atoms command");
      if (strcmp(arg[iarg+1],"yes") == 0) mol_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) mol_flag = 0;
      else error->all(FLERR,"Illegal delete_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal delete_atoms command");
  }
}
