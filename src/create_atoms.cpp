/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_atoms.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "error.h"

#define MAXATOMS 0x7FFFFFFF
#define BIG      1.0e30
#define EPSILON  1.0e-6

/* ---------------------------------------------------------------------- */

void CreateAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Create_atoms command before simulation box is defined");
  if (domain->lattice == NULL)
    error->all("Cannot create atoms with undefined lattice");

  // parse arguments

  int nbasis = domain->lattice->nbasis;
  int basistype[nbasis];

  if (narg < 1) error->all("Illegal create_atoms command");
  int itype = atoi(arg[0]);
  if (itype <= 0 || itype > atom->ntypes) 
    error->all("Invalid atom type in create_atoms command");
  for (int i = 0; i < nbasis; i++) basistype[i] = itype;

  regionflag = -1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal create_atoms command");
      int iregion;
      for (iregion = 0; iregion < domain->nregion; iregion++)
	if (strcmp(arg[iarg+1],domain->regions[iregion]->id) == 0) break;
      if (iregion == domain->nregion)
	error->all("Create_atoms region ID does not exist");
      regionflag = iregion;
      iarg += 2;
    } else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all("Illegal create_atoms command");
      int ibasis = atoi(arg[iarg+1]);
      itype = atoi(arg[iarg+2]);
      if (ibasis <= 0 || ibasis > nbasis || 
	  itype <= 0 || itype > atom->ntypes) 
	error->all("Illegal create_atoms command");
      basistype[ibasis-1] = itype;
      iarg += 3;
    } else error->all("Illegal create_atoms command");
  }

  // convert 8 corners of my sub-box from box coords to lattice coords
  // min to max = bounding box around the pts in lattice space

  subxlo = domain->subxlo;
  subxhi = domain->subxhi;
  subylo = domain->subylo;
  subyhi = domain->subyhi;
  subzlo = domain->subzlo;
  subzhi = domain->subzhi;

  double xmin,ymin,zmin,xmax,ymax,zmax;
  xmin = ymin = zmin = BIG;
  xmax = ymax = zmax = -BIG;

  domain->lattice->bbox(1,subxlo,subylo,subzlo,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxhi,subylo,subzlo,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxlo,subyhi,subzlo,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxhi,subyhi,subzlo,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxlo,subylo,subzhi,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxhi,subylo,subzhi,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxlo,subyhi,subzhi,xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,subxhi,subyhi,subzhi,xmin,ymin,zmin,xmax,ymax,zmax);

  // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my sub-box
  // overlap = any part of a unit cell (face,edge,pt) in common with my sub-box
  // in lattice space, sub-box is a tilted box
  // but bbox of sub-box is aligned with lattice axes
  // so ilo:khi unit cells should completely tile bounding box
  // decrement lo values if min < 0, since static_cast(-1.5) = -1

  int ilo,ihi,jlo,jhi,klo,khi;
  ilo = static_cast<int> (xmin);
  jlo = static_cast<int> (ymin);
  klo = static_cast<int> (zmin);
  ihi = static_cast<int> (xmax);
  jhi = static_cast<int> (ymax);
  khi = static_cast<int> (zmax);

  if (xmin < 0.0) ilo--;
  if (ymin < 0.0) jlo--;
  if (zmin < 0.0) klo--;

  // iterate on 3d periodic lattice using loop bounds
  // invoke add_atom for nbasis atoms in each unit cell
  // add_atom converts lattice coords to box coords, checks if in my sub-box

  boxxhi = domain->boxxhi;
  boxyhi = domain->boxyhi;
  boxzhi = domain->boxzhi;

  double natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  double **basis = domain->lattice->basis;

  int i,j,k,m;
  for (k = klo; k <= khi; k++)
    for (j = jlo; j <= jhi; j++)
      for (i = ilo; i <= ihi; i++)
	for (m = 0; m < nbasis; m++)
	  add_atom(basistype[m],i+basis[m][0],j+basis[m][1],k+basis[m][2]);

  // new total # of atoms

  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);

  // print status

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Created %.15g atoms\n",atom->natoms-natoms_previous);
    if (logfile)
      fprintf(logfile,"Created %.15g atoms\n",atom->natoms-natoms_previous);
  }

  // reset simulation now that more atoms are defined
  // add tags for newly created atoms if possible
  // if global map exists, reset it
  // if a molecular system, set nspecial to 0 for new atoms

  if (atom->natoms > MAXATOMS) atom->tag_enable = 0;
  if (atom->natoms <= MAXATOMS) atom->tag_extend();

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
  if (atom->molecular) {
    int **nspecial = atom->nspecial;
    for (i = nlocal_previous; i < atom->nlocal; i++) {
      nspecial[i][0] = 0;
      nspecial[i][1] = 0;
      nspecial[i][2] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   add an atom of type at lattice coords x,y,z if it meets all criteria 
------------------------------------------------------------------------- */

void CreateAtoms::add_atom(int type, double x, double y, double z)
{
  // convert from lattice coords to box coords

  domain->lattice->lattice2box(x,y,z);

  // if a region was specified, test if atom is in it

  if (regionflag >= 0)
    if (!domain->regions[regionflag]->match(x,y,z)) return;

  // test if atom is in my subbox

  if (x < subxlo || x >= subxhi || 
      y < subylo || y >= subyhi || 
      z < subzlo || z >= subzhi) return;

  // don't put atoms within EPSILON of upper periodic boundary
  // else may overlap image atom at lower boundary

  if (domain->xperiodic && fabs(x-boxxhi) < EPSILON) return;
  if (domain->yperiodic && fabs(y-boxyhi) < EPSILON) return;
  if (domain->zperiodic && fabs(z-boxzhi) < EPSILON) return;

  // add the atom to my list of atoms

  atom->create_one(type,x,y,z);
}
