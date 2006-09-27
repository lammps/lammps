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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "math.h"
#include "domain.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "region.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#define RegionInclude
#include "style.h"
#undef RegionInclude

#define BIG   1.0e20
#define SMALL 1.0e-4
#define DELTA 1
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ----------------------------------------------------------------------
   default is periodic 
------------------------------------------------------------------------- */

Domain::Domain()
{
  box_exist = 0;

  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;

  boxxlo = boxylo = boxzlo = -0.5;
  boxxhi = boxyhi = boxzhi = 0.5;

  char *str = "none";
  int n = strlen(str) + 1;
  lattice_style = new char[n];
  strcpy(lattice_style,str);

  origin_x = origin_y = origin_z = 0.0;
  orient_x[0] = 1; orient_x[1] = 0; orient_x[2] = 0;
  orient_y[0] = 0; orient_y[1] = 1; orient_y[2] = 0;
  orient_z[0] = 0; orient_z[1] = 0; orient_z[2] = 1;

  nregion = maxregion = 0;
  regions = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete [] lattice_style;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  // set box_change if box dimensions ever change
  // due to shrink-wrapping, fix nph, npt, volume/rescale, uniaxial

  box_change = 0;
  if (nonperiodic == 2) box_change = 1;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nph") == 0) box_change = 1;
    if (strcmp(modify->fix[i]->style,"npt") == 0) box_change = 1;
    if (strcmp(modify->fix[i]->style,"volume/rescale") == 0) box_change = 1;
    if (strcmp(modify->fix[i]->style,"uniaxial") == 0) box_change = 1;
  }
}

/* ----------------------------------------------------------------------
   setup initial global box, boxxlo-boxzhi are already set
   adjust for any shrink-wrapped boundaries
   store min values if boundary type = 3 = "m"
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  if (boundary[0][0] == 2) boxxlo -= SMALL;
  else if (boundary[0][0] == 3) minxlo = boxxlo;
  if (boundary[0][1] == 2) boxxhi += SMALL;
  else if (boundary[0][1] == 3) minxhi = boxxhi;

  if (boundary[1][0] == 2) boxylo -= SMALL;
  else if (boundary[1][0] == 3) minylo = boxylo;
  if (boundary[1][1] == 2) boxyhi += SMALL;
  else if (boundary[1][1] == 3) minyhi = boxyhi;

  if (boundary[2][0] == 2) boxzlo -= SMALL;
  else if (boundary[2][0] == 3) minzlo = boxzlo;
  if (boundary[2][1] == 2) boxzhi += SMALL;
  else if (boundary[2][1] == 3) minzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   setup global box parameters
   set prd, prd_half, prd[], boxlo/hi[], periodicity[]
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  xprd_half = 0.5*xprd;
  yprd_half = 0.5*yprd;
  zprd_half = 0.5*zprd;

  prd[0]   = xprd;    prd[1]   = yprd;    prd[2]   = zprd;
  boxlo[0] = boxxlo;  boxlo[1] = boxylo;  boxlo[2] = boxzlo;
  boxhi[0] = boxxhi;  boxhi[1] = boxyhi;  boxhi[2] = boxzhi;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;
}

/* ----------------------------------------------------------------------
   set local subbox from global boxxlo-boxzhi and proc grid
   set subxlo-subzhi, sublo/hi[] 
   for uppermost proc, insure subhi = boxhi (in case round-off occurs)
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  subxlo = boxxlo + comm->myloc[0] * xprd / comm->procgrid[0];
  if (comm->myloc[0] < comm->procgrid[0]-1)
    subxhi = boxxlo + (comm->myloc[0]+1) * xprd / comm->procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + comm->myloc[1] * yprd / comm->procgrid[1];
  if (comm->myloc[1] < comm->procgrid[1]-1)
    subyhi = boxylo + (comm->myloc[1]+1) * yprd / comm->procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo + comm->myloc[2] * zprd / comm->procgrid[2];
  if (comm->myloc[2] < comm->procgrid[2]-1)
    subzhi = boxzlo + (comm->myloc[2]+1) * zprd / comm->procgrid[2];
  else subzhi = boxzhi;

  sublo[0] = subxlo;  sublo[1] = subylo;  sublo[2] = subzlo;
  subhi[0] = subxhi;  subhi[1] = subyhi;  subhi[2] = subzhi;
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxxlo thru boxzhi
   call set_global_box and set_local_box 
------------------------------------------------------------------------- */

void Domain::reset_box()
{
  if (nonperiodic == 2) {

    // compute extent of atoms on this proc

    double extent[3][2],all[3][2];

    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;
    
    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      extent[0][0] = MIN(extent[0][0],x[i][0]);
      extent[0][1] = MAX(extent[0][1],x[i][0]);
      extent[1][0] = MIN(extent[1][0],x[i][1]);
      extent[1][1] = MAX(extent[1][1],x[i][1]);
      extent[2][0] = MIN(extent[2][0],x[i][2]);
      extent[2][1] = MAX(extent[2][1],x[i][2]);
    }

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX
    // set box by extent in shrink-wrapped dims

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // if any of 6 dims is shrink-wrapped, set it to extent of atoms +/- SMALL
    // enforce min extent for boundary type = 3 = "m"
    
    if (xperiodic == 0) {
      if (boundary[0][0] == 2) boxxlo = -all[0][0] - SMALL;
      else if (boundary[0][0] == 3) boxxlo = MIN(-all[0][0]-SMALL,minxlo);
      if (boundary[0][1] == 2) boxxhi = all[0][1] + SMALL;
      else if (boundary[0][1] == 3) boxxhi = MAX(all[0][1]+SMALL,minxhi);
    }
    if (yperiodic == 0) {
      if (boundary[1][0] == 2) boxylo = -all[1][0] - SMALL;
      else if (boundary[1][0] == 3) boxylo = MIN(-all[1][0]-SMALL,minylo);
      if (boundary[1][1] == 2) boxyhi = all[1][1] + SMALL;
      else if (boundary[1][1] == 3) boxyhi = MAX(all[1][1]+SMALL,minyhi);
    }
    if (zperiodic == 0) {
      if (boundary[2][0] == 2) boxzlo = -all[2][0] - SMALL;
      else if (boundary[2][0] == 3) boxzlo = MIN(-all[2][0]-SMALL,minzlo);
      if (boundary[2][1] == 2) boxzhi = all[2][1] + SMALL;
      else if (boundary[2][1] == 3) boxzhi = MAX(all[2][1]+SMALL,minzhi);
    }
  }

  set_global_box();
  set_local_box();
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::pbc()
{
  int i,idim,otherdims;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int *image = atom->image;

  if (xperiodic) {
    for (i = 0; i < nlocal; i++) {
      if (x[i][0] < boxxlo) {
	x[i][0] += xprd;
	idim = image[i] & 1023;
        otherdims = image[i] ^ idim;
	idim--;
	idim &= 1023;
	image[i] = otherdims | idim;
      }
      if (x[i][0] >= boxxhi) {
	x[i][0] -= xprd;
	x[i][0] = MAX(x[i][0],boxxlo);
	idim = image[i] & 1023;
	otherdims = image[i] ^ idim;
	idim++;
	idim &= 1023;
	image[i] = otherdims | idim;
      }
    }
  }

  if (yperiodic) {
    for (i = 0; i < nlocal; i++) {
      if (x[i][1] < boxylo) {
	x[i][1] += yprd;
	idim = (image[i] >> 10) & 1023;
        otherdims = image[i] ^ (idim << 10);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 10);
      }
      if (x[i][1] >= boxyhi) {
	x[i][1] -= yprd;
	x[i][1] = MAX(x[i][1],boxylo);
	idim = (image[i] >> 10) & 1023;
        otherdims = image[i] ^ (idim << 10);
	idim++;
	idim &= 1023;
	image[i] = otherdims | (idim << 10);
      }
    }
  }

  if (zperiodic) {
    for (i = 0; i < nlocal; i++) {
      if (x[i][2] < boxzlo) {
	x[i][2] += zprd;
	idim = image[i] >> 20;
        otherdims = image[i] ^ (idim << 20);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 20);
      }
      if (x[i][2] >= boxzhi) {
	x[i][2] -= zprd;
	x[i][2] = MAX(x[i][2],boxzlo);
	idim = image[i] >> 20;
        otherdims = image[i] ^ (idim << 20);
	idim++;
	idim &= 1023;
	image[i] = otherdims | (idim << 20);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test 
------------------------------------------------------------------------- */

void Domain::minimum_image(double *dx, double *dy, double *dz)
{
  if (xperiodic) {
    if (fabs(*dx) > xprd_half) {
      if (*dx < 0.0) *dx += xprd;
      else *dx -= xprd;
    }
  }
  if (yperiodic) {
    if (fabs(*dy) > yprd_half) {
      if (*dy < 0.0) *dy += yprd;
      else *dy -= yprd;
    }
  }
  if (zperiodic) {
    if (fabs(*dz) > zprd_half) {
      if (*dz < 0.0) *dz += zprd;
      else *dz -= zprd;
    }
  }
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test 
------------------------------------------------------------------------- */

void Domain::minimum_image(double *x)
{
  if (xperiodic) {
    if (fabs(x[0]) > xprd_half) {
      if (x[0] < 0.0) x[0] += xprd;
      else x[0] -= xprd;
    }
  }
  if (yperiodic) {
    if (fabs(x[1]) > yprd_half) {
      if (x[1] < 0.0) x[1] += yprd;
      else x[1] -= yprd;
    }
  }
  if (zperiodic) {
    if (fabs(x[2]) > zprd_half) {
      if (x[2] < 0.0) x[2] += zprd;
      else x[2] -= zprd;
    }
  }
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   adjust image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap(double &x, double &y, double &z, int &image)
{
  if (xperiodic) {
    while (x < boxxlo) {
      x += xprd;
      int idim = image & 1023;
      int otherdims = image ^ idim;
      idim--;
      idim &= 1023;
      image = otherdims | idim;
    }
    while (x >= boxxhi) {
      x -= xprd;
      int idim = image & 1023;
      int otherdims = image ^ idim;
      idim++;
      idim &= 1023;
      image = otherdims | idim;
    }
    x = MAX(x,boxxlo);
  }

  if (yperiodic) {
    while (y < boxylo) {
      y += yprd;
      int idim = (image >> 10) & 1023;
      int otherdims = image ^ (idim << 10);
      idim--;
      idim &= 1023;
      image = otherdims | (idim << 10);
    }
    while (y >= boxyhi) {
      y -= yprd;
      int idim = (image >> 10) & 1023;
      int otherdims = image ^ (idim << 10);
      idim++;
      idim &= 1023;
      image = otherdims | (idim << 10);
    }
    y = MAX(y,boxylo);
  }

  if (zperiodic) {
    while (z < boxzlo) {
      z += zprd;
      int idim = image >> 20;
      int otherdims = image ^ (idim << 20);
      idim--;
      idim &= 1023;
      image = otherdims | (idim << 20);
    }
    while (z >= boxzhi) {
      z -= zprd;
      int idim = image >> 20;
      int otherdims = image ^ (idim << 20);
      idim--;
      idim &= 1023;
      image = otherdims | (idim << 20);
    }
    z = MAX(z,boxzlo);
  }
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   don't reset image flag
------------------------------------------------------------------------- */

void Domain::unmap(double &x, double &y, double &z, int image)
{
  int xbox = (image & 1023) - 512;
  int ybox = (image >> 10 & 1023) - 512;
  int zbox = (image >> 20) - 512;

  x = x + xbox*xprd;
  y = y + ybox*yprd;
  z = z + zbox*zprd;
}

/* ----------------------------------------------------------------------
   set lattice constant 
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal lattice command");

  delete [] lattice_style;
  int n = strlen(arg[0]) + 1;
  lattice_style = new char[n];
  strcpy(lattice_style,arg[0]);

  if (strcmp(arg[0],"none") == 0) return;

  if (narg != 2) error->all("Illegal lattice command");

  if (strcmp(arg[0],"sc") && strcmp(arg[0],"bcc") && strcmp(arg[0],"fcc") &&
      strcmp(arg[0],"sq") && strcmp(arg[0],"sq2") && strcmp(arg[0],"hex") &&
      strcmp(arg[0],"diamond"))
    error->all("Illegal lattice command");

  // check that lattice matches dimension

  int dim = force->dimension;

  if (dim == 2) {
    if (strcmp(arg[0],"sq") && strcmp(arg[0],"sq2") && strcmp(arg[0],"hex"))
      error->all("Lattice style incompatible with dimension");
  }

  if (dim == 3) {
    if (strcmp(arg[0],"sc") && strcmp(arg[0],"bcc") &&
	strcmp(arg[0],"fcc") && strcmp(arg[0],"diamond"))
      error->all("Lattice style incompatible with dimension");
  }

  // set lattice constants depending on # of atoms per unit cell
  // hex is only case where xlattice = ylattice = zlattice is not true

  double value = atof(arg[1]);

  if (!strcmp(arg[0],"sc") || !strcmp(arg[0],"sq")) {
    if (strcmp(update->unit_style,"lj") == 0)
      xlattice = ylattice = zlattice = pow(1.0/value,1.0/dim);
    else xlattice = ylattice = zlattice = value;
  }

  if (!strcmp(arg[0],"bcc") || !strcmp(arg[0],"sq2")) {
    if (strcmp(update->unit_style,"lj") == 0)
      xlattice = ylattice = zlattice = pow(2.0/value,1.0/dim);
    else xlattice = ylattice = zlattice = value;
  }

  if (strcmp(arg[0],"fcc") == 0) {
    if (strcmp(update->unit_style,"lj") == 0)
      xlattice = ylattice = zlattice = pow(4.0/value,1.0/dim);
    else xlattice = ylattice = zlattice = value;
  }

  if (strcmp(arg[0],"hex") == 0) {
    if (strcmp(update->unit_style,"lj") == 0)
      xlattice = zlattice = pow(2.0/sqrt(3.0)/value,1.0/dim);
    else xlattice = zlattice = value;
    ylattice = xlattice * sqrt(3.0);
  }

  if (strcmp(arg[0],"diamond") == 0) {
    if (strcmp(update->unit_style,"lj") == 0)
      xlattice = ylattice = zlattice = pow(8.0/value,1.0/dim);
    else xlattice = ylattice = zlattice = value;
  }
}

/* ----------------------------------------------------------------------
   check if orientation vectors are mutually orthogonal 
------------------------------------------------------------------------- */

int Domain::orthogonality()
{
  if (orient_x[0]*orient_y[0] + orient_x[1]*orient_y[1] + 
      orient_x[2]*orient_y[2]) return 0;

  if (orient_y[0]*orient_z[0] + orient_y[1]*orient_z[1] + 
      orient_y[2]*orient_z[2]) return 0;

  if (orient_x[0]*orient_z[0] + orient_x[1]*orient_z[1] + 
      orient_x[2]*orient_z[2]) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   check righthandedness of orientation vectors
   x cross y must be in same direction as z 
------------------------------------------------------------------------- */

int Domain::right_handed()
{
  int xy0 = orient_x[1]*orient_y[2] - orient_x[2]*orient_y[1];
  int xy1 = orient_x[2]*orient_y[0] - orient_x[0]*orient_y[2];
  int xy2 = orient_x[0]*orient_y[1] - orient_x[1]*orient_y[0];
  if (xy0*orient_z[0] + xy1*orient_z[1] + xy2*orient_z[2] <= 0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   convert lattice coords into box coords
   x,y,z = point in lattice coords
   orient_xyz = lattice vectors that point in positive x,y,z box directions
   origin_xyz = origin of lattice in lattice units
   return xnew,ynew,znew = point in box coords
   method:
     compute projection of vector from (0,0,0) to lattice onto each
       orient vector via dot product scaled by length of
       orient vector in lattice coords
     this projection (offset by origin and scaled by lattice constant)
       gives x,y,z box coords 
------------------------------------------------------------------------- */

void Domain::lattice2box(double *x, double *y, double *z)
{
  double length;
  double dotprod;

  length = orient_x[0]*orient_x[0] + orient_x[1]*orient_x[1] +
    orient_x[2]*orient_x[2];
  length = sqrt(length);
  dotprod = *x * orient_x[0] + *y * orient_x[1] + *z * orient_x[2];
  double xnew = (dotprod/length + origin_x) * xlattice;

  length = orient_y[0]*orient_y[0] + orient_y[1]*orient_y[1] +
    orient_y[2]*orient_y[2];
  length = sqrt(length);
  dotprod = *x * orient_y[0] + *y * orient_y[1] + *z * orient_y[2];
  double ynew = (dotprod/length + origin_y) * ylattice;

  length = orient_z[0]*orient_z[0] + orient_z[1]*orient_z[1] +
     orient_z[2]*orient_z[2];
  length = sqrt(length);
  dotprod = *x * orient_z[0] + *y * orient_z[1] + *z * orient_z[2];
  double znew = (dotprod/length + origin_z) * zlattice;

  *x = xnew;
  *y = ynew;
  *z = znew;
}

/* ----------------------------------------------------------------------
   create a new region 
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal region command");

  // error checks

  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(arg[0],regions[iregion]->id) == 0)
      error->all("Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTA;
    regions = (Region **) 
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if (strcmp(arg[1],"none") == 0) error->all("Invalid region style");

#define RegionClass
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) regions[nregion] = new Class(narg,arg);
#include "style.h"
#undef RegionClass

  else error->all("Invalid region style");

  nregion++;
}

/* ----------------------------------------------------------------------
   boundary settings from the input script 
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal boundary command");

  char c;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      if (c == 'p') boundary[idim][iside] = 0;
      else if (c == 'f') boundary[idim][iside] = 1;
      else if (c == 's') boundary[idim][iside] = 2;
      else if (c == 'm') boundary[idim][iside] = 3;
      else error->all("Illegal boundary command");
    }

  for (int idim = 0; idim < 3; idim++)
    if ((boundary[idim][0] == 0 && boundary[idim][1]) ||
	(boundary[idim][0] && boundary[idim][1] == 0))
      error->all("Both sides of boundary must be periodic");

  if (boundary[0][0] == 0) xperiodic = 1;
  else xperiodic = 0;
  if (boundary[1][0] == 0) yperiodic = 1;
  else yperiodic = 0;
  if (boundary[2][0] == 0) zperiodic = 1;
  else zperiodic = 0;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
	boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
	boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }
}
