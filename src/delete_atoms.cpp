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
#include "comm.h"
#include "force.h"
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

  // allocate and initialize deletion list
  
  int nlocal = atom->nlocal;
  int *dlist = new int[nlocal];

  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // delete the atoms

  if (strcmp(arg[0],"group") == 0) delete_group(narg,arg,dlist);
  else if (strcmp(arg[0],"region") == 0) delete_region(narg,arg,dlist);
  else if (strcmp(arg[0],"overlap") == 0) delete_overlap(narg,arg,dlist);
  else error->all("Illegal delete_atoms command");

  // delete local atoms in list
  // reset nlocal

  AtomVec *avec = atom->avec;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }
  atom->nlocal = nlocal;
  delete [] dlist;

  // if non-molecular system, reset atom tags to be contiguous
  // set all atom IDs to 0, call tag_extend()

  if (atom->molecular == 0) {
    int *tag = atom->tag;
    int nlocal = atom->nlocal;
 
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

void DeleteAtoms::delete_group(int narg, char **arg, int *dlist)
{
  if (narg != 2) error->all("Illegal delete_atoms command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find delete_atoms group ID");

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete all atoms in region
------------------------------------------------------------------------- */

void DeleteAtoms::delete_region(int narg, char **arg, int *dlist)
{
  if (narg != 2) error->all("Illegal delete_atoms command");
  
  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all("Could not find delete_atoms region ID");

  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) dlist[i] = 1;
}

/* ----------------------------------------------------------------------
   delete atoms so there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

void DeleteAtoms::delete_overlap(int narg, char **arg, int *dlist)
{
  if (narg < 2) error->all("Illegal delete_atoms command");
    
  // read args including optional type info

  double cut = atof(arg[1]);
  double cutsq = cut*cut;

  int typeflag,type1,type2;
  if (narg == 2) typeflag = 0;
  else if (narg == 3) {
    typeflag = 1;
    type1 = atoi(arg[2]);
  } else if (narg == 4) {
    typeflag = 2;
    type1 = atoi(arg[2]);
    type2 = atoi(arg[3]);
  } else error->all("Illegal delete_atoms command");

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for delete_atoms ...\n");

  // request a half neighbor list for use by this command

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->occasional = 1;

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  
  lmp->init();

  // error check on cutoff

  if (cut > neighbor->cutghost) 
    error->all("Delete_atoms cutoff > ghost cutoff");

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

  // create an atom map if one doesn't exist already

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  // double loop over owned atoms and their neighbors
  // at end of loop, there are no more overlaps
  // criteria for i,j to undergo a deletion event:
  //   weighting factor != 0.0 for this pair
  //     could be 0 and still be in neigh list for long-range Coulombics
  //   local copy of j (map[tag[j]]) has not already been deleted
  //   distance between i,j is less than cutoff
  //   i,j are of valid types
  // if all criteria met, delete i and skip to next i in outer loop
  //   unless j is ghost and newton_pair off and tag[j] < tag[i]
  //   then rely on other proc to delete

  int *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  int i,j,ii,jj,m,inum,jnum,itype,jtype;
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
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j >= nall) {
	if (special_coul[j/nall] == 0.0 && special_lj[j/nall] == 0.0) continue;
	j %= nall;
      }

      if (j < nlocal) {
	if (dlist[j]) continue;
      } else {
	m = atom->map(tag[j]);
	if (m < nlocal && dlist[m]) continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;

      if (typeflag) {
	jtype = type[j];
	if (typeflag == 1 && itype != type1 && jtype != type1) continue;
	if (typeflag == 2 && !(itype == type1 && jtype == type2) && 
	    !(itype == type2 && jtype == type1)) continue;
      }

      if (j >= nlocal && newton_pair == 0 && tag[j] < tag[i]) continue;

      dlist[i] = 1;
      break;
    }
  }

  // delete temporary atom map

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }
}
