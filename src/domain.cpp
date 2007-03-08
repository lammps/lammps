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
#include "stdio.h"
#include "math.h"
#include "domain.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "region.h"
#include "lattice.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#define RegionInclude
#include "style.h"
#undef RegionInclude

using namespace LAMMPS_NS;

#define BIG   1.0e20
#define SMALL 1.0e-4
#define DELTA 1
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ----------------------------------------------------------------------
   default is periodic 
------------------------------------------------------------------------- */

Domain::Domain(LAMMPS *lmp) : Pointers(lmp)
{
  box_exist = 0;

  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;

  triclinic = 0;
  boxxlo = boxylo = boxzlo = -0.5;
  boxxhi = boxyhi = boxzhi = 0.5;
  xy = xz = yz = 0.0;

  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;

  lattice = NULL;
  nregion = maxregion = 0;
  regions = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete lattice;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  // set box_change if box dimensions/shape ever changes
  // due to shrink-wrapping, fixes that change volume (npt, vol/rescale, etc)

  box_change = 0;
  if (nonperiodic == 2) box_change = 1;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->box_change) box_change = 1;
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  // error checks for orthogonal and triclinic domains

  if (boxxlo >= boxxhi || boxylo >= boxyhi || boxzlo >= boxzhi)
    error->one("Box bounds are invalid");

  if (triclinic) {
    if (force->dimension == 2 && (xz != 0.0 || yz != 0.0))
      error->all("Cannot skew triclinic box in z for 2d simulation");
    if (xy != 0.0 && (!xperiodic || !yperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");
    if (xz != 0.0 && (!xperiodic || !zperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");
    if (yz != 0.0 && (!yperiodic || !zperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");

    if (fabs(xy/(boxyhi-boxylo)) > 0.5)
      error->all("Triclinic box skew is too large");
    if (fabs(xz/(boxzhi-boxzlo)) > 0.5)
      error->all("Triclinic box skew is too large");
    if (fabs(yz/(boxzhi-boxzlo)) > 0.5)
      error->all("Triclinic box skew is too large");
  }

  // adjust box lo/hi for shrink-wrapped dims

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
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxxhi - boxxlo;
  prd[1] = yprd = boxyhi - boxylo;
  prd[2] = zprd = boxzhi - boxzlo;

  xprd_half = 0.5*xprd;
  yprd_half = 0.5*yprd;
  zprd_half = 0.5*zprd;

  boxlo[0] = boxxlo;  boxlo[1] = boxylo;  boxlo[2] = boxzlo;
  boxhi[0] = boxxhi;  boxhi[1] = boxyhi;  boxhi[2] = boxzhi;

  if (triclinic) {
    h[0] = xprd;
    h[1] = yprd;
    h[2] = zprd;
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
  
    h_inv[0] = 1.0/h[0];
    h_inv[1] = 1.0/h[1];
    h_inv[2] = 1.0/h[2];
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);
    
    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];

    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }
}

/* ----------------------------------------------------------------------
   set lamda box params, only need be done one time
   assumes global box is defined and proc assignment has been made by comm
   for uppermost proc, insure subhi = 1.0 (in case round-off occurs)
------------------------------------------------------------------------- */

void Domain::set_lamda_box()
{
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;

  sublo_lamda[0] = 1.0*myloc[0] / procgrid[0];
  sublo_lamda[1] = 1.0*myloc[1] / procgrid[1];
  sublo_lamda[2] = 1.0*myloc[2] / procgrid[2];

  if (myloc[0] < procgrid[0]-1)
    subhi_lamda[0] = 1.0*(myloc[0]+1) / procgrid[0];
  else subhi_lamda[0] = 1.0;

  if (myloc[1] < procgrid[1]-1)
    subhi_lamda[1] = 1.0*(myloc[1]+1) / procgrid[1];
  else subhi_lamda[1] = 1.0;

  if (myloc[2] < procgrid[2]-1)
    subhi_lamda[2] = 1.0*(myloc[2]+1) / procgrid[2];
  else subhi_lamda[2] = 1.0;
}

/* ----------------------------------------------------------------------
   set local subbox params
   assumes global box is defined and proc assignment has been made
   for uppermost proc, insure subhi = boxhi (in case round-off occurs)
   for triclinic, lo/hi_bound is a bbox around all 8 tilted pts of sub-box
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;

  if (domain->triclinic == 0) {
    sublo[0] = boxlo[0] + myloc[0] * xprd / procgrid[0];
    if (myloc[0] < procgrid[0]-1)
      subhi[0] = boxlo[0] + (myloc[0]+1) * xprd / procgrid[0];
    else subhi[0] = boxhi[0];
    
    sublo[1] = boxlo[1] + myloc[1] * yprd / procgrid[1];
    if (myloc[1] < procgrid[1]-1)
      subhi[1] = boxlo[1] + (myloc[1]+1) * yprd / procgrid[1];
    else subhi[1] = boxhi[1];

    sublo[2] = boxlo[2] + myloc[2] * zprd / procgrid[2];
    if (myloc[2] < procgrid[2]-1)
      subhi[2] = boxlo[2] + (myloc[2]+1) * zprd / procgrid[2];
    else subhi[2] = boxhi[2];

  } else {
    sublo[0] = h[0]*sublo_lamda[0] + h[5]*sublo_lamda[1] + 
      h[4]*sublo_lamda[2] + boxlo[0];
    sublo[1] = h[1]*sublo_lamda[1] + h[3]*sublo_lamda[2] + boxlo[1];
    sublo[2] = h[2]*sublo_lamda[2] + boxlo[2];

    subhi[0] = sublo[0] + xprd/procgrid[0];
    subhi[1] = sublo[1] + yprd/procgrid[1];
    subhi[2] = sublo[2] + zprd/procgrid[2];

    sublo_bound[0] = MIN(sublo[0],sublo[0] + xy/procgrid[1]);
    sublo_bound[0] = MIN(sublo_bound[0],sublo_bound[0] + xz/procgrid[2]);
    sublo_bound[1] = MIN(sublo[1],sublo[1] + yz/procgrid[2]);
    sublo_bound[2] = sublo[2];

    subhi_bound[0] = MAX(subhi[0],subhi[0] + xy/procgrid[1]);
    subhi_bound[0] = MAX(subhi_bound[0],subhi_bound[0] + xz/procgrid[2]);
    subhi_bound[1] = MAX(subhi[1],subhi[1] + yz/procgrid[2]);
    subhi_bound[2] = subhi[2];
  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   shrink-wrapping only occurs in non-periodic, non-triclinic dims
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

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // in shrink-wrapped dims, set box by atom extent
    // if set, observe min box size settings
    
    if (xperiodic == 0) {
      if (boundary[0][0] == 2) boxxlo = -all[0][0] - SMALL;
      else if (boundary[0][0] == 3) boxxlo = MIN(-all[0][0]-SMALL,minxlo);
      if (boundary[0][1] == 2) boxxhi = all[0][1] + SMALL;
      else if (boundary[0][1] == 3) boxxhi = MAX(all[0][1]+SMALL,minxhi);
      if (boxxlo > boxxhi) error->all("Illegal simulation box");
    }
    if (yperiodic == 0) {
      if (boundary[1][0] == 2) boxylo = -all[1][0] - SMALL;
      else if (boundary[1][0] == 3) boxylo = MIN(-all[1][0]-SMALL,minylo);
      if (boundary[1][1] == 2) boxyhi = all[1][1] + SMALL;
      else if (boundary[1][1] == 3) boxyhi = MAX(all[1][1]+SMALL,minyhi);
      if (boxylo > boxyhi) error->all("Illegal simulation box");
    }
    if (zperiodic == 0) {
      if (boundary[2][0] == 2) boxzlo = -all[2][0] - SMALL;
      else if (boundary[2][0] == 3) boxzlo = MIN(-all[2][0]-SMALL,minzlo);
      if (boundary[2][1] == 2) boxzhi = all[2][1] + SMALL;
      else if (boundary[2][1] == 3) boxzhi = MAX(all[2][1]+SMALL,minzhi);
      if (boxzlo > boxzhi) error->all("Illegal simulation box");
    }
  }

  set_global_box();
  set_local_box();
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::pbc()
{
  int i,idim,otherdims;
  double *lo,*hi,*period;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int *image = atom->image;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
	x[i][0] += period[0];
	idim = image[i] & 1023;
        otherdims = image[i] ^ idim;
	idim--;
	idim &= 1023;
	image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
	x[i][0] -= period[0];
	x[i][0] = MAX(x[i][0],lo[0]);
	idim = image[i] & 1023;
	otherdims = image[i] ^ idim;
	idim++;
	idim &= 1023;
	image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i][1] < lo[1]) {
	x[i][1] += period[1];
	idim = (image[i] >> 10) & 1023;
        otherdims = image[i] ^ (idim << 10);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 10);
      }
      if (x[i][1] >= hi[1]) {
	x[i][1] -= period[1];
	x[i][1] = MAX(x[i][1],lo[1]);
	idim = (image[i] >> 10) & 1023;
        otherdims = image[i] ^ (idim << 10);
	idim++;
	idim &= 1023;
	image[i] = otherdims | (idim << 10);
      }
    }

    if (zperiodic) {
      if (x[i][2] < lo[2]) {
	x[i][2] += period[2];
	idim = image[i] >> 20;
        otherdims = image[i] ^ (idim << 20);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 20);
      }
      if (x[i][2] >= hi[2]) {
	x[i][2] -= period[2];
	x[i][2] = MAX(x[i][2],lo[2]);
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
   for triclinic, also add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::minimum_image(double &dx, double &dy, double &dz)
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
	if (dy < 0.0) dy += yprd;
	else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
	if (dz < 0.0) dz += zprd;
	else dz -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
	if (dz < 0.0) {
	  dz += zprd;
	  dy += yz;
	  dx += xz;
	} else {
	  dz -= zprd;
	  dy -= yz;
	  dx -= xz;
	}
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
	if (dy < 0.0) {
	  dy += yprd;
	  dx += xy;
	} else {
	  dy -= yprd;
	  dx -= xy;
	}
      }
    }
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::minimum_image(double *delta)
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(delta[0]) > xprd_half) {
	if (delta[0] < 0.0) delta[0] += xprd;
	else delta[0] -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(delta[1]) > yprd_half) {
	if (delta[1] < 0.0) delta[1] += yprd;
	else delta[1] -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(delta[2]) > zprd_half) {
	if (delta[2] < 0.0) delta[2] += zprd;
	else delta[2] -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(delta[2]) > zprd_half) {
	if (delta[2] < 0.0) {
	  delta[2] += zprd;
	  delta[1] += yz;
	  delta[0] += xz;
	} else {
	  delta[2] -= zprd;
	  delta[1] -= yz;
	  delta[0] -= xz;
	}
      }
    }
    if (yperiodic) {
      if (fabs(delta[1]) > yprd_half) {
	if (delta[1] < 0.0) {
	  delta[1] += yprd;
	  delta[0] += xy;
	} else {
	  delta[1] -= yprd;
	  delta[0] -= xy;
	}
      }
    }
    if (xperiodic) {
      if (fabs(delta[0]) > xprd_half) {
	if (delta[0] < 0.0) delta[0] += xprd;
	else delta[0] -= xprd;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   adjust image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, convert atom to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap(double *x, int &image)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      int idim = image & 1023;
      int otherdims = image ^ idim;
      idim--;
      idim &= 1023;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      int idim = image & 1023;
      int otherdims = image ^ idim;
      idim++;
      idim &= 1023;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      int idim = (image >> 10) & 1023;
      int otherdims = image ^ (idim << 10);
      idim--;
      idim &= 1023;
      image = otherdims | (idim << 10);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      int idim = (image >> 10) & 1023;
      int otherdims = image ^ (idim << 10);
      idim++;
      idim &= 1023;
      image = otherdims | (idim << 10);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      int idim = image >> 20;
      int otherdims = image ^ (idim << 20);
      idim--;
      idim &= 1023;
      image = otherdims | (idim << 20);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      int idim = image >> 20;
      int otherdims = image ^ (idim << 20);
      idim++;
      idim &= 1023;
      image = otherdims | (idim << 20);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(double *x, int image)
{
  int xbox = (image & 1023) - 512;
  int ybox = (image >> 10 & 1023) - 512;
  int zbox = (image >> 20) - 512;

  if (triclinic == 0) {
    x[0] += xbox*xprd;
    x[1] += ybox*yprd;
    x[2] += zbox*zprd;
  } else {
    x[0] += h[0]*xbox + h[5]*ybox + h[4]*zbox;
    x[1] += h[1]*ybox + h[3]*zbox;
    x[2] += h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   create a lattice
   delete it if style = none
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  delete lattice;
  lattice = new Lattice(lmp,narg,arg);
  if (lattice->style == 0) {
    delete lattice;
    lattice = NULL;
  }
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
  else if (strcmp(arg[1],#key) == 0) \
    regions[nregion] = new Class(lmp,narg,arg);
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

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
	boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
	boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }
}

/* ----------------------------------------------------------------------
   print box info, orthogonal or triclinic
------------------------------------------------------------------------- */

void Domain::print_box(char *str)
{
  if (comm->me == 0) {
    if (screen) {
      if (domain->triclinic == 0)
	fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
		str,boxxlo,boxylo,boxzlo,boxxhi,boxyhi,boxzhi);
      else {
	char *format = 
	  "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
	fprintf(screen,format,
		str,boxxlo,boxylo,boxzlo,boxxhi,boxyhi,boxzhi,xy,xz,yz);
      }
    }
    if (logfile) {
      if (triclinic == 0)
	fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
		str,boxxlo,boxylo,boxzlo,boxxhi,boxyhi,boxzhi);
      else {
	char *format = 
	  "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
	fprintf(logfile,format,
		str,boxxlo,boxylo,boxzlo,boxxhi,boxyhi,boxzhi,xy,xz,yz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void Domain::lamda2x(int n)
{
  double **x = atom->x;

  for (int i = 0; i < n; i++) { 
    x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
    x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
    x[i][2] = h[2]*x[i][2] + boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void Domain::x2lamda(int n)
{
  double delta[3];
  double **x = atom->x;

  for (int i = 0; i < n; i++) { 
    delta[0] = x[i][0] - boxlo[0];
    delta[1] = x[i][1] - boxlo[1];
    delta[2] = x[i][2] - boxlo[2];

    x[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
    x[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
    x[i][2] = h_inv[2]*delta[2];
  }
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for one atom
   x = H lamda + x0;
   lamda and x can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::lamda2x(double *lamda, double *x)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}
