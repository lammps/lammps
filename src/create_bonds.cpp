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
#include "create_bonds.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "comm.h"
#include "group.h"
#include "special.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateBonds::CreateBonds(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBonds::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_bonds command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use create_bonds unless atoms have IDs");
  if (atom->molecular != 1) 
    error->all(FLERR,"Cannot use create_bonds with non-molecular system");

  if (narg != 5) error->all(FLERR,"Illegal create_bonds command");

  // parse args

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Cannot find create_bonds group ID");
  int group1bit = group->bitmask[igroup];
  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Cannot find create_bonds group ID");
  int group2bit = group->bitmask[igroup];

  int btype = force->inumeric(FLERR,arg[2]);
  double rmin = force->numeric(FLERR,arg[3]);
  double rmax = force->numeric(FLERR,arg[4]);

  if (btype <= 0 || btype > atom->nbondtypes) 
    error->all(FLERR,"Invalid bond type in create_bonds command");
  if (rmin > rmax) error->all(FLERR,"Illegal create_bonds command");

  double rminsq = rmin*rmin;
  double rmaxsq = rmax*rmax;

  // store state before bond creation

  bigint nbonds_previous = atom->nbonds;

  // request a full neighbor list for use by this command

  int irequest = neighbor->request(this);
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
    error->all(FLERR,"Create_bonds requires a pair style be defined");
  if (rmax > neighbor->cutneighmax)
    error->all(FLERR,"Create_bonds max distance > neighbor cutoff");
  if (rmax > neighbor->cutneighmin && comm->me == 0)
    error->warning(FLERR,"Create_bonds max distance > minimum neighbor cutoff");

  // require special_bonds 1-2 weights = 0.0 and KSpace = NULL
  // so that already bonded atom pairs do not appear in neighbor list
  // otherwise with newton_bond = 1,
  //   would be hard to check if I-J bond already existed
  // note that with KSpace, pair with weight = 0 could still be in neigh list

  if (force->special_lj[1] != 0.0 || force->special_coul[1] != 0.0)
    error->all(FLERR,"Create_bonds command requires "
               "special_bonds 1-2 weights be 0.0");
  if (force->kspace_style)
    error->all(FLERR,"Create_bonds command requires "
               "no kspace_style be defined");

  // setup domain, communication and neighboring
  // acquire ghosts and build standard neighbor lists

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  neighbor->build();

  // build neighbor list this command needs based on earlier request

  NeighList *list = neighbor->lists[irequest];
  neighbor->build_one(list);

  // loop over all neighs of each atom
  // compute distance between two atoms consistently on both procs
  // add bond if group and distance criteria are met
  // check that bond list does not overflow

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  double **x = atom->x;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  double newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;
  
  int i,j,ii,jj,inum,jnum,flag;
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
      j &= NEIGHMASK;

      // only consider bond creation if I,J distance between 2 cutoffs
      // compute rsq identically on both I,J loop iterations
      // if I,J tags equal, do not bond atom to itself
      
      if (tag[i] < tag[j]) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
      } else if (tag[i] > tag[j]) {
        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;
      } else continue;
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rminsq || rsq > rmaxsq) continue;

      // only consider bond creation if igroup and jgroup match I,J atoms

      flag = 0;
      if ((mask[i] & group1bit) && (mask[j] & group2bit)) flag = 1;
      if ((mask[i] & group2bit) && (mask[j] & group1bit)) flag = 1;
      if (!flag) continue;

      // create bond, check for overflow
      // on I,J loop iterations, store with 1 or 2 atoms based on newton_bond

      if (!newton_bond || tag[i] < tag[j]) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR,
                     "New bond exceeded bonds per atom in create_bonds");
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tag[j];
        num_bond[i]++;
      }
    }
  }

  // recount bonds

  bigint nbonds = 0;
  for (int i = 0; i < nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (!force->newton_bond) atom->nbonds /= 2;

  // print new bond count

  bigint nadd_bonds = atom->nbonds - nbonds_previous;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Added " BIGINT_FORMAT
              " bonds, new total = " BIGINT_FORMAT "\n",
              nadd_bonds,atom->nbonds);
    }

    if (logfile) {
      fprintf(logfile,"Added " BIGINT_FORMAT
              " bonds, new total = " BIGINT_FORMAT "\n",
              nadd_bonds,atom->nbonds);
    }
  }

  // re-trigger special list build

  Special special(lmp);
  special.build();
}
