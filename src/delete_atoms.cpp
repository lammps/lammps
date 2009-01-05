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

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DeleteAtoms::DeleteAtoms(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DeleteAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Delete_atoms command before simulation box is defined");
  if (narg < 1) error->all("Illegal delete_atoms command");

  // store state before delete

  if (atom->tag_enable == 0)
    error->all("Cannot use delete_atoms unless atoms have IDs");

  int natoms_previous = static_cast<int> (atom->natoms);

  // delete the atoms

  if (strcmp(arg[0],"group") == 0) delete_group(narg,arg);
  else if (strcmp(arg[0],"region") == 0) delete_region(narg,arg);
  else if (strcmp(arg[0],"overlap") == 0) delete_overlap(narg,arg);
  else if (strcmp(arg[0],"porosity") == 0) delete_porosity(narg,arg);
  else error->all("Illegal delete_atoms command");

  // delete local atoms flagged in dlist
  // reset nlocal

  AtomVec *avec = atom->avec;
  int nlocal = atom->nlocal;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->sfree(dlist);

  // if non-molecular system, reset atom tags to be contiguous
  // set all atom IDs to 0, call tag_extend()

  if (atom->molecular == 0) {
    int *tag = atom->tag;
    for (i = 0; i < nlocal; i++) tag[i] = 0;
    atom->tag_extend();
  }

  // reset atom->natoms
  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped

  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // print before and after atom count

  int ndelete = static_cast<int> (natoms_previous - atom->natoms);

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Deleted %d atoms, new total = %.15g\n",
			ndelete,atom->natoms);
    if (logfile) fprintf(logfile,"Deleted %d atoms, new total = %.15g\n",
			 ndelete,atom->natoms);
  }
}

/* ----------------------------------------------------------------------
   delete all atoms in group
   group will still exist
------------------------------------------------------------------------- */

void DeleteAtoms::delete_group(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal delete_atoms command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find delete_atoms group ID");

  // allocate and initialize deletion list
  
  int nlocal = atom->nlocal;
  dlist = (int *) memory->smalloc(nlocal*sizeof(int),"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete all atoms in region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_region(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal delete_atoms command");
  
  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all("Could not find delete_atoms region ID");

  // allocate and initialize deletion list
  
  int nlocal = atom->nlocal;
  dlist = (int *) memory->smalloc(nlocal*sizeof(int),"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++)
    if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete atoms so there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

void DeleteAtoms::delete_overlap(int narg, char **arg)
{
  if (narg < 4) error->all("Illegal delete_atoms command");
    
  // read args

  double cut = atof(arg[1]);
  double cutsq = cut*cut;

  int igroup1 = group->find(arg[2]);
  int igroup2 = group->find(arg[3]);
  if (igroup1 < 0 || igroup2 < 0)
    error->all("Could not find delete_atoms group ID");

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
    error->all("Delete_atoms requires a pair style be defined");
  if (cut > neighbor->cutneighmax) 
    error->all("Delete_atoms cutoff > neighbor cutoff");

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
  dlist = (int *) memory->smalloc(nlocal*sizeof(int),"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // double loop over owned atoms and their full neighbor list
  // at end of loop, there are no more overlaps
  // only ever delete owned atom I, never J even if owned

  int *tag = atom->tag;
  int *mask = atom->mask;
  double **x = atom->x;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

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

      // if weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (j >= nall) {
	if (special_coul[j/nall] == 0.0 && special_lj[j/nall] == 0.0) continue;
	j %= nall;
      }

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
  if (narg != 4) error->all("Illegal delete_atoms command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all("Could not find delete_atoms region ID");

  double porosity_fraction = atof(arg[2]);
  int seed = atoi(arg[3]);

  RanMars *random = new RanMars(lmp,seed + comm->me);

  // allocate and initialize deletion list
 
  int nlocal = atom->nlocal;
  dlist = (int *) memory->smalloc(nlocal*sizeof(int),"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++)
    if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
      if (random->uniform() <= porosity_fraction) dlist[i] = 1;
}


