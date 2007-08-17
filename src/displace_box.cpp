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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "displace_box.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,FINAL,DELTA,SCALE,VOLUME};
enum{ONE_FROM_ONE,ONE_FROM_TWO,TWO_FROM_ONE};
enum{NO_REMAP,X_REMAP};

/* ---------------------------------------------------------------------- */

DisplaceBox::DisplaceBox(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DisplaceBox::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0) 
    error->all("Displace_box command before simulation box is defined");
  if (narg < 2) error->all("Illegal displace_box command");

  // init entire system since comm->irregular is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for displace_box ...\n");
  lmp->init();

  if (comm->me == 0 && screen) fprintf(screen,"Displacing box ...\n");

  // group

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all("Could not find displace_box group ID");
  int groupbit = group->bitmask[igroup];

  // set defaults

  set = new Set[6];
  set[0].style = set[1].style = set[2].style = 
    set[3].style = set[4].style = set[5].style = NONE;

  // parse arguments

  int triclinic = domain->triclinic;

  int index;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0 || strcmp(arg[iarg],"y") == 0 ||
	strcmp(arg[iarg],"z") == 0) {
      if (strcmp(arg[iarg],"x") == 0) index = 0;
      else if (strcmp(arg[iarg],"y") == 0) index = 1;
      else if (strcmp(arg[iarg],"z") == 0) index = 2;

      if (iarg+2 > narg) error->all("Illegal displace_box command");
      if (strcmp(arg[iarg+1],"final") == 0) {
	if (iarg+4 > narg) error->all("Illegal displace_box command");
	set[index].style = FINAL;
	set[index].flo = atof(arg[iarg+2]);
	set[index].fhi = atof(arg[iarg+3]);
	iarg += 4;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
	if (iarg+4 > narg) error->all("Illegal displace_box command");
	set[index].style = DELTA;
	set[index].dlo = atof(arg[iarg+2]);
	set[index].dhi = atof(arg[iarg+3]);
	iarg += 4;
      } else if (strcmp(arg[iarg+1],"scale") == 0) {
	if (iarg+3 > narg) error->all("Illegal displace_box command");
	set[index].style = SCALE;
	set[index].scale = atof(arg[iarg+2]);
	iarg += 3;
      } else if (strcmp(arg[iarg+1],"volume") == 0) {
	set[index].style = VOLUME;
	iarg += 2;
      } else error->all("Illegal displace_box command");

    } else if (strcmp(arg[iarg],"xy") == 0 || strcmp(arg[iarg],"xz") == 0 ||
	strcmp(arg[iarg],"yz") == 0) {
      if (triclinic == 0)
	error->all("Displace_box tilt factors require triclinic box");
      if (strcmp(arg[iarg],"xy") == 0) index = 5;
      else if (strcmp(arg[iarg],"xz") == 0) index = 4;
      else if (strcmp(arg[iarg],"yz") == 0) index = 3;
      if (iarg+2 > narg) error->all("Illegal displace_box command");
      if (strcmp(arg[iarg+1],"final") == 0) {
	if (iarg+3 > narg) error->all("Illegal displace_box command");
	set[index].style = FINAL;
	set[index].ftilt = atof(arg[iarg+2]);
	iarg += 3;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
	if (iarg+3 > narg) error->all("Illegal displace_box command");
	set[index].style = DELTA;
	set[index].dtilt = atof(arg[iarg+2]);
	iarg += 3;
      } else error->all("Illegal displace_box command");

    } else break;
  }

  // read options from end of input line

  options(narg-iarg,&arg[iarg]);

  // check periodicity

  if ((set[0].style && domain->xperiodic == 0) ||
      (set[1].style && domain->yperiodic == 0) ||
      (set[2].style && domain->zperiodic == 0))
    error->all("Cannot displace_box on a non-periodic boundary");

  if (set[3].style && (domain->yperiodic == 0 || domain->zperiodic == 0))
    error->all("Cannot displace_box on a non-periodic boundary");
  if (set[4].style && (domain->xperiodic == 0 || domain->zperiodic == 0))
    error->all("Cannot displace_box on a non-periodic boundary");
  if (set[5].style && (domain->xperiodic == 0 || domain->yperiodic == 0))
    error->all("Cannot displace_box on a non-periodic boundary");

  // apply scaling to FINAL,DELTA since they have distance units

  int flag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style == FINAL || set[i].style == DELTA) flag = 1;

  if (flag && scaleflag && domain->lattice == NULL)
    error->all("Use of displace_box with undefined lattice");

  double xscale,yscale,zscale;
  if (flag && scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // for 3,4,5 scaling is in 1st dimension, e.g. x for xz

  double map[6];
  map[0] = xscale; map[1] = yscale; map[2] = zscale;
  map[3] = yscale; map[4] = xscale; map[5] = xscale;

  for (int i = 0; i < 3; i++) {
    if (set[i].style == FINAL) {
      set[i].flo *= map[i];
      set[i].fhi *= map[i];
    } else if (set[i].style == DELTA) {
      set[i].dlo *= map[i];
      set[i].dhi *= map[i];
    }
  }

  for (int i = 3; i < 6; i++) {
    if (set[i].style == FINAL) set[i].ftilt *= map[i];
    else if (set[i].style == DELTA) set[i].dtilt *= map[i];
  }

  // set initial/final values for box size and shape
  // final = initial if no setting

  for (int i = 0; i < 3; i++) {
    set[i].lo_stop = set[i].lo_start = domain->boxlo[i];
    set[i].hi_stop = set[i].hi_start = domain->boxhi[i];

    if (set[i].style == FINAL) {
      set[i].lo_stop = set[i].flo;
      set[i].hi_stop = set[i].fhi;
    } else if (set[i].style == DELTA) {
      set[i].lo_stop = set[i].lo_start + set[i].dlo;
      set[i].hi_stop = set[i].hi_start + set[i].dhi;
    } else if (set[i].style == SCALE) {
      set[i].lo_stop = 0.5*(set[i].lo_start+set[i].hi_start) - 
	0.5*set[i].scale*(set[i].hi_start-set[i].lo_start);
      set[i].hi_stop = 0.5*(set[i].lo_start+set[i].hi_start) + 
	0.5*set[i].scale*(set[i].hi_start-set[i].lo_start);
    }
  }

  for (int i = 3; i < 6; i++) {
    if (i == 5) set[i].tilt_start = domain->xy;
    else if (i == 4) set[i].tilt_start = domain->xz;
    else if (i == 3) set[i].tilt_start = domain->yz;
    set[i].tilt_stop = set[i].tilt_start;

    if (set[i].style == FINAL) {
      set[i].tilt_stop = set[i].ftilt;
    } else if (set[i].style == DELTA) {
      set[i].tilt_stop = set[i].tilt_start + set[i].dtilt;
    }
  }

  // for VOLUME, setup links to other dims
  // fixed, dynamic1,2, vol_start

  for (int i = 0; i < 3; i++) {
    set[i].vol_start = domain->xprd * domain->yprd * domain->zprd;

    if (set[i].style != VOLUME) continue;
    int other1 = (i+1) % 3;
    int other2 = (i+2) % 3;

    if (set[other1].style == NONE) {
      if (set[other2].style == NONE || set[other2].style == VOLUME)
	error->all("Fix deform volume setting is invalid");
      set[i].substyle = ONE_FROM_ONE;
      set[i].fixed = other1;
      set[i].dynamic1 = other2;
    } else if (set[other2].style == NONE) {
      if (set[other1].style == NONE || set[other1].style == VOLUME)
	error->all("Fix deform volume setting is invalid");
      set[i].substyle = ONE_FROM_ONE;
      set[i].fixed = other2;
      set[i].dynamic1 = other1;
    } else if (set[other1].style == VOLUME) {
      if (set[other2].style == NONE || set[other2].style == VOLUME)
	error->all("Fix deform volume setting is invalid");
      set[i].substyle = TWO_FROM_ONE;
      set[i].fixed = other1;
      set[i].dynamic1 = other2;
    } else if (set[other2].style == VOLUME) {
      if (set[other1].style == NONE || set[other1].style == VOLUME)
	error->all("Fix deform volume setting is invalid");
      set[i].substyle = TWO_FROM_ONE;
      set[i].fixed = other2;
      set[i].dynamic1 = other1;
    } else {
      set[i].substyle = ONE_FROM_TWO;
      set[i].dynamic2 = other1;
      set[i].dynamic2 = other2;
    }
  }

  // set new box size for VOLUME dims that are linked to other dims

  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;

    if (set[i].substyle == ONE_FROM_ONE) {
      set[i].lo_stop = 0.5*(set[i].lo_start+set[i].hi_start) -
	0.5*(set[i].vol_start /
	     (set[set[i].dynamic1].hi_stop -
	      set[set[i].dynamic1].lo_stop) /
	     (set[set[i].fixed].hi_start-set[set[i].fixed].lo_start));
      set[i].hi_stop = 0.5*(set[i].lo_start+set[i].hi_start) +
	0.5*(set[i].vol_start /
	     (set[set[i].dynamic1].hi_stop -
	      set[set[i].dynamic1].lo_stop) /
	     (set[set[i].fixed].hi_start-set[set[i].fixed].lo_start));

    } else if (set[i].substyle == ONE_FROM_TWO) {
      set[i].lo_stop = 0.5*(set[i].lo_start+set[i].hi_start) -
	0.5*(set[i].vol_start /
	     (set[set[i].dynamic1].hi_stop - 
	      set[set[i].dynamic1].lo_stop) /
	     (set[set[i].dynamic2].hi_stop - 
	      set[set[i].dynamic2].lo_stop));
      set[i].hi_stop = 0.5*(set[i].lo_start+set[i].hi_start) +
	0.5*(set[i].vol_start /
	     (set[set[i].dynamic1].hi_stop -
	      set[set[i].dynamic1].lo_stop) /
	     (set[set[i].dynamic2].hi_stop - 
	      set[set[i].dynamic2].lo_stop));
      
    } else if (set[i].substyle == TWO_FROM_ONE) {
      set[i].lo_stop = 0.5*(set[i].lo_start+set[i].hi_start) -
	0.5*sqrt(set[i].vol_start /
		 (set[set[i].dynamic1].hi_stop - 
		  set[set[i].dynamic1].lo_stop) /
		 (set[set[i].fixed].hi_start - 
		  set[set[i].fixed].lo_start) *
		 (set[i].hi_start - set[i].lo_start));
      set[i].hi_stop = 0.5*(set[i].lo_start+set[i].hi_start) +
	0.5*sqrt(set[i].vol_start /
		 (set[set[i].dynamic1].hi_stop - 
		  set[set[i].dynamic1].lo_stop) /
		 (set[set[i].fixed].hi_start - 
		  set[set[i].fixed].lo_start) *
		 (set[i].hi_start - set[i].lo_start));
    }
  }

  // check that final tilt is not illegal value

  double xprd_stop = set[0].hi_stop - set[0].lo_stop;
  double yprd_stop = set[0].hi_stop - set[0].lo_stop;

  if (set[3].tilt_stop < -0.5*yprd_stop || set[3].tilt_stop > 0.5*yprd_stop ||
      set[4].tilt_stop < -0.5*xprd_stop || set[4].tilt_stop > 0.5*xprd_stop ||
      set[5].tilt_stop < -0.5*xprd_stop || set[5].tilt_stop > 0.5*xprd_stop)
    error->all("Induced tilt by displace_box is too large");

  // convert atoms to lamda coords

  if (remapflag == X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	domain->x2lamda(x[i],x[i]);
 }

  // reset global and local box to new size/shape

  domain->boxlo[0] = set[0].lo_stop;
  domain->boxlo[1] = set[1].lo_stop;
  domain->boxlo[2] = set[2].lo_stop;
  domain->boxhi[0] = set[0].hi_stop;
  domain->boxhi[1] = set[1].hi_stop;
  domain->boxhi[2] = set[2].hi_stop;

  if (triclinic) {
    domain->yz = set[3].tilt_stop;
    domain->xz = set[4].tilt_stop;
    domain->xy = set[5].tilt_stop;
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms back to box coords

  if (remapflag == X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	domain->lamda2x(x[i],x[i]);
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case box
  //   moved a long distance relative to atoms
  // use irregular() instead of exchange() in case box
  //   moved a long distance relative to atoms

  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  comm->setup();
  comm->irregular();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // clean up

  delete [] set;

  // check if any atoms were lost

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (natoms != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via displace_box: original %.15g current %.15g",
	    atom->natoms,natoms);
    error->all(str);
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_box input line 
------------------------------------------------------------------------- */

void DisplaceBox::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal displace_box command");

  remapflag = X_REMAP;
  scaleflag = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all("Illegal displace_box command");
      if (strcmp(arg[iarg+1],"x") == 0) remapflag = X_REMAP;
      else if (strcmp(arg[iarg+1],"none") == 0) remapflag = NO_REMAP;
      else error->all("Illegal displace_box command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal displace_box command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal displace_box command");
      iarg += 2;
    } else error->all("Illegal displace_box command");
  }
}
