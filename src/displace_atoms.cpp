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
#include "displace_atoms.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "group.h"
#include "random_park.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{MOVE,RAMP,RANDOM};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DisplaceAtoms::DisplaceAtoms(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DisplaceAtoms::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0) 
    error->all("Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all("Illegal displace_atoms command");

  // init entire system since comm->irregular is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for displace_atoms ...\n");
  lmp->init();

  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");

  // group and style

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all("Could not find displace_atoms group ID");
  int groupbit = group->bitmask[igroup];

  int style;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  else error->all("Illegal displace_atoms command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == RANDOM) options(narg-6,&arg[6]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of displace_atoms with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector

  if (style == MOVE) {

    double delx = xscale*atof(arg[2]);
    double dely = yscale*atof(arg[3]);
    double delz = zscale*atof(arg[4]);

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	x[i][0] += delx;
	x[i][1] += dely;
	x[i][2] += delz;
      }
    }
  }

  // move atoms in ramped fashion
    
  if (style == RAMP) {

    int d_dim;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all("Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*atof(arg[3]);
      d_hi = xscale*atof(arg[4]);
    } else if (d_dim == 1) {
      d_lo = yscale*atof(arg[3]);
      d_hi = yscale*atof(arg[4]);
    } else if (d_dim == 2) {
      d_lo = zscale*atof(arg[3]);
      d_hi = zscale*atof(arg[4]);
    }

    int coord_dim;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all("Illegal displace_atoms ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*atof(arg[6]);
      coord_hi = xscale*atof(arg[7]);
    } else if (coord_dim == 1) {
      coord_lo = yscale*atof(arg[6]);
      coord_hi = yscale*atof(arg[7]);
    } else if (coord_dim == 2) {
      coord_lo = zscale*atof(arg[6]);
      coord_hi = zscale*atof(arg[7]);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double fraction,dramp;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
	fraction = MAX(fraction,0.0);
	fraction = MIN(fraction,1.0);
	dramp = d_lo + fraction*(d_hi - d_lo);
	x[i][d_dim] += dramp;
      }
    }
  }

  // move atoms randomly
  // makes atom result independent of what proc owns it via random->reset()

    
  if (style == RANDOM) {
    RanPark *random = new RanPark(lmp,1);

    double dx = xscale*atof(arg[2]);
    double dy = yscale*atof(arg[3]);
    double dz = zscale*atof(arg[4]);
    int seed = atoi(arg[5]);
    if (seed <= 0) error->all("Illegal displace_atoms random command");

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	random->reset(seed,x[i]);
	x[i][0] += dx * 2.0*(random->uniform()-0.5);
	x[i][1] += dy * 2.0*(random->uniform()-0.5);
	x[i][2] += dz * 2.0*(random->uniform()-0.5);
      }
    }

    delete random;
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() instead of exchange() in case atoms moved a long distance

  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  comm->setup();
  comm->irregular();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (natoms != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via displace_atoms: original %.15g current %.15g",
	    atom->natoms,natoms);
    error->all(str);
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_atoms input line 
------------------------------------------------------------------------- */

void DisplaceAtoms::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal displace_atoms command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal displace_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal displace_atoms command");
      iarg += 2;
    } else error->all("Illegal displace_atoms command");
  }
}
