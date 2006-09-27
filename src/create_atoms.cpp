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
#include "force.h"
#include "domain.h"
#include "update.h"
#include "region.h"
#include "error.h"

#define SC  1
#define BCC 2
#define FCC 3
#define SQ  4
#define SQ2 5
#define HEX 6
#define DIAMOND 7

#define MINATOMS 1000
#define MAXATOMS 0x7FFFFFFF
#define EPSILON  1.0e-6

/* ---------------------------------------------------------------------- */

void CreateAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0) 
    error->all("Create_atoms command before simulation box is defined");

  if (narg != 1 && narg != 2) error->all("Illegal create_atoms command");

  create_type = atoi(arg[0]);
  if (create_type > atom->ntypes) 
    error->all("Too large an atom type in create_atoms command");

  if (strcmp(domain->lattice_style,"none") == 0)
    error->all("Cannot create atoms with undefined lattice");
  if (!domain->orthogonality()) error->all("Non-orthogonal lattice vectors");
  if (!domain->right_handed())
    error->all("Orientation vectors are not right-handed");

  // iregion = specified region (-1 if not specified)

  iregion = -1;
  if (narg == 2) {
    for (iregion = 0; iregion < domain->nregion; iregion++)
      if (strcmp(arg[1],domain->regions[iregion]->id) == 0) break;
    if (iregion == domain->nregion)
      error->all("Create_atoms region ID does not exist");
  }

  // local copies of domain properties

  subxlo = domain->subxlo;
  subxhi = domain->subxhi;
  subylo = domain->subylo;
  subyhi = domain->subyhi;
  subzlo = domain->subzlo;
  subzhi = domain->subzhi;

  boxxhi = domain->boxxhi;
  boxyhi = domain->boxyhi;
  boxzhi = domain->boxzhi;

  // ilo:ihi,jlo:jhi,klo:khi = loop bounds of simple cubic lattice
  //   that entirely overlaps my proc's sub-box

  int ilo,ihi,jlo,jhi,klo,khi;

  loop_bounds(0,&ilo,&ihi);
  loop_bounds(1,&jlo,&jhi);
  loop_bounds(2,&klo,&khi);

  // initialize 3d periodic lattice using overlapping loop bounds
  // lattice style determines how many atoms in cubic unit cell
  // sc = 1, bcc = 2, fcc = 4, sq = 1, sq2 = 2, hex = 2, diamond = 8

  double natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  int style;
  if (strcmp(domain->lattice_style,"sc") == 0) style = SC;
  else if (strcmp(domain->lattice_style,"bcc") == 0) style = BCC;
  else if (strcmp(domain->lattice_style,"fcc") == 0) style = FCC;
  else if (strcmp(domain->lattice_style,"sq") == 0) style = SQ;
  else if (strcmp(domain->lattice_style,"sq2") == 0) style = SQ2;
  else if (strcmp(domain->lattice_style,"hex") == 0) style = HEX;
  else if (strcmp(domain->lattice_style,"diamond") == 0) style = DIAMOND;

  double ifull,ihalf,jfull,jhalf,kfull,khalf;
  double iquart,i3quart,jquart,j3quart,kquart,k3quart;
  int i,j,k;

  for (k = klo; k <= khi; k++) {
    kfull = (double) k;
    khalf = k + 0.5;
    kquart = k + 0.25;
    k3quart = k + 0.75;
    for (j = jlo; j <= jhi; j++) {
      jfull = (double) j;
      jhalf = j + 0.5;
      jquart = j + 0.25;
      j3quart = j + 0.75;
      for (i = ilo; i <= ihi; i++) {
	ifull = (double) i;
	ihalf = i + 0.5;
	iquart = i + 0.25;
	i3quart = i + 0.75;

	if (style == SC)
	  add_atom(ifull,jfull,kfull);
	else if (style == BCC) {
	  add_atom(ifull,jfull,kfull);
	  add_atom(ihalf,jhalf,khalf);
	} else if (style == FCC) {
	  add_atom(ifull,jfull,kfull);
	  add_atom(ihalf,jhalf,kfull);
	  add_atom(ihalf,jfull,khalf);
	  add_atom(ifull,jhalf,khalf);
	} else if (style == SQ) {
	  add_atom(ifull,jfull,kfull);
	} else if (style == SQ2) {
	  add_atom(ifull,jfull,kfull);
	  add_atom(ihalf,jhalf,kfull);
	} else if (style == HEX) {
	  add_atom(ifull,jfull,kfull);
	  add_atom(ihalf,jhalf,kfull);
	} else if (style == DIAMOND) {
	  add_atom(ifull,jfull,kfull);
	  add_atom(ifull,jhalf,khalf);
	  add_atom(ihalf,jfull,khalf);
	  add_atom(ihalf,jhalf,kfull);
	  add_atom(iquart,jquart,kquart);
	  add_atom(iquart,j3quart,k3quart);
	  add_atom(i3quart,jquart,k3quart);
	  add_atom(i3quart,j3quart,kquart);
	}
      }
    }
  }

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
   add an atom at x,y,z in lattice coords if it meets all criteria 
------------------------------------------------------------------------- */

void CreateAtoms::add_atom(double x, double y, double z)
{
  // convert from lattice coords to box coords

  domain->lattice2box(&x,&y,&z);

  // if a region was specified, test if atom is in it

  if (iregion >= 0)
    if (!domain->regions[iregion]->match(x,y,z)) return;

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

  atom->create_one(create_type,x,y,z);
}

/* ----------------------------------------------------------------------
   search for 2 bounding lattice planes that completely enclose my sub-box
   do this by testing if all corner points of my sub-box lie on correct side
     of a lattice plane via same_side function
   dim = 0,1,2 for x,y,z directions
   lo,hi = returned indices of 2 bounding lattice planes
   a lattice plane is defined by:
     point = point in lattice space through which the plane passes
     normal = vector normal to lattice plane
------------------------------------------------------------------------- */

void CreateAtoms::loop_bounds(int dim, int *lo, int *hi)
{
  int normal[3],point[3];

  // start search at origin

  point[0] = point[1] = point[2] = 0;

  // set lattice plane direction along positive lattice axis

  normal[0] = normal[1] = normal[2] = 0;
  normal[dim] = 1;

  // step down (if needed) until entire box is above the plane
  // step up until 1st time entire box is not above the plane

  while (!same_side(point,normal)) point[dim]--;
  while (same_side(point,normal)) point[dim]++;

  // lower loop bound = current loc minus 1 (subtract 1 more for safety)

  *lo = point[dim] - 2;

  // flip plane direction
  // step up until entire box is below the plane

  normal[dim] = -1;
  while (!same_side(point,normal)) point[dim]++;

  // lower loop bound = current loc (add 1 more for safety)

  *hi = point[dim] + 1;
}

/* ----------------------------------------------------------------------
   test if all 8 corner points of my sub-box are on "correct" side of a plane
   plane is defined by point[3] it goes thru and a normal[3]
   normal also defines the correct side 
------------------------------------------------------------------------- */

int CreateAtoms::same_side(int *point, int *normal)
{
  // p1 = plane center point in box coords
  // p2 = point on correct side of plane, in box coords

  double p1x = point[0];
  double p1y = point[1];
  double p1z = point[2];
  domain->lattice2box(&p1x,&p1y,&p1z);

  double p2x = point[0] + normal[0];
  double p2y = point[1] + normal[1];
  double p2z = point[2] + normal[2];
  domain->lattice2box(&p2x,&p2y,&p2z);

  // for each of 8 sub-box corner points, dot these 2 vectors:
  //   v1 = from plane center point to point on correct side of plane
  //   v2 = from plane center point to box corner point
  // negative result = portion of box is on wrong side of plane, return 0

  double v1[3],v2[3];

  points2vec(p1x,p1y,p1z,p2x,p2y,p2z,v1);

  points2vec(p1x,p1y,p1z,subxlo,subylo,subzlo,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxhi,subylo,subzlo,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxlo,subyhi,subzlo,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxhi,subyhi,subzlo,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxlo,subylo,subzhi,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxhi,subylo,subzhi,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxlo,subyhi,subzhi,v2);
  if (dot(v1,v2) < 0.0) return 0;
  points2vec(p1x,p1y,p1z,subxhi,subyhi,subzhi,v2);
  if (dot(v1,v2) < 0.0) return 0;

  // all 8 points were on correct side

  return 1;
}

/* ---------------------------------------------------------------------- */
				 
double CreateAtoms::dot(double *vec1, double *vec2)
{
  double sum = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
  return sum;
}

/* ---------------------------------------------------------------------- */

void CreateAtoms::points2vec(double p1x, double p1y, double p1z,
			     double p2x, double p2y, double p2z, double *v)
{
  v[0] = p2x - p1x;
  v[1] = p2y - p1y;
  v[2] = p2z - p1z;
}
