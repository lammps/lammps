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

/* ----------------------------------------------------------------------
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
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
#include "fix_deform.h"
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

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ----------------------------------------------------------------------
   default is periodic 
------------------------------------------------------------------------- */

Domain::Domain(LAMMPS *lmp) : Pointers(lmp)
{
  box_exist = 0;

  dimension = 3;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;

  triclinic = 0;
  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  xy = xz = yz = 0.0;

  h[3] = h[4] = h[5] = 0.0;
  h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
  h_rate[0] = h_rate[1] = h_rate[2] = 
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
  
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

  // check for fix deform

  deform_flag = deform_remap = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      deform_flag = 1;
      if (((FixDeform *) modify->fix[i])->remapflag == V_REMAP) {
	deform_remap = 1;
	deform_groupbit = modify->fix[i]->groupbit;
      } else deform_remap = 0;
    }
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  // error checks for orthogonal and triclinic domains

  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one("Box bounds are invalid");

  if (triclinic) {
    if (domain->dimension == 2 && (xz != 0.0 || yz != 0.0))
      error->all("Cannot skew triclinic box in z for 2d simulation");
    if (xy != 0.0 && (!xperiodic || !yperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");
    if (xz != 0.0 && (!xperiodic || !zperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");
    if (yz != 0.0 && (!yperiodic || !zperiodic))
      error->all("Triclinic box must be periodic in skewed dimensions");

    if (fabs(xy/(boxhi[0]-boxlo[0])) > 0.5)
      error->all("Triclinic box skew is too large");
    if (fabs(xz/(boxhi[0]-boxlo[0])) > 0.5)
      error->all("Triclinic box skew is too large");
    if (fabs(yz/(boxhi[1]-boxlo[1])) > 0.5)
      error->all("Triclinic box skew is too large");
  }

  // adjust box lo/hi for shrink-wrapped dims

  if (boundary[0][0] == 2) boxlo[0] -= SMALL;
  else if (boundary[0][0] == 3) minxlo = boxlo[0];
  if (boundary[0][1] == 2) boxhi[0] += SMALL;
  else if (boundary[0][1] == 3) minxhi = boxhi[0];

  if (boundary[1][0] == 2) boxlo[1] -= SMALL;
  else if (boundary[1][0] == 3) minylo = boxlo[1];
  if (boundary[1][1] == 2) boxhi[1] += SMALL;
  else if (boundary[1][1] == 3) minyhi = boxhi[1];

  if (boundary[2][0] == 2) boxlo[2] -= SMALL;
  else if (boundary[2][0] == 3) minzlo = boxlo[2];
  if (boundary[2][1] == 2) boxhi[2] += SMALL;
  else if (boundary[2][1] == 3) minzhi = boxhi[2];
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  xprd_half = 0.5*xprd;
  yprd_half = 0.5*yprd;
  zprd_half = 0.5*zprd;

  if (triclinic) {
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
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
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;

  if (triclinic == 0) {
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
  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   if shrink-wrapped and triclinic, perform shrink-wrap in box coords
------------------------------------------------------------------------- */

void Domain::reset_box()
{
  if (nonperiodic == 2) {

    // convert back to box coords for shrink-wrap operation

    if (triclinic) lamda2x(atom->nlocal);

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
    // if minimum set, enforce min box size settings

    if (xperiodic == 0) {
      if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - SMALL;
      else if (boundary[0][0] == 3) boxlo[0] = MIN(-all[0][0]-SMALL,minxlo);
      if (boundary[0][1] == 2) boxhi[0] = all[0][1] + SMALL;
      else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+SMALL,minxhi);
      if (boxlo[0] > boxhi[0]) error->all("Illegal simulation box");
    }
    if (yperiodic == 0) {
      if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - SMALL;
      else if (boundary[1][0] == 3) boxlo[1] = MIN(-all[1][0]-SMALL,minylo);
      if (boundary[1][1] == 2) boxhi[1] = all[1][1] + SMALL;
      else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+SMALL,minyhi);
      if (boxlo[1] > boxhi[1]) error->all("Illegal simulation box");
    }
    if (zperiodic == 0) {
      if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - SMALL;
      else if (boundary[2][0] == 3) boxlo[2] = MIN(-all[2][0]-SMALL,minzlo);
      if (boundary[2][1] == 2) boxhi[2] = all[2][1] + SMALL;
      else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+SMALL,minzhi);
      if (boxlo[2] > boxhi[2]) error->all("Illegal simulation box");
    }
  }

  set_global_box();
  set_local_box();

  // if shrink-wrapped, convert to lamda coords for new box
  // must re-invoke pbc() b/c x2lamda result can be outside 0,1 due to roundoff

  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
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
  double **v = atom->v;
  int *mask = atom->mask;
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
	if (deform_remap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
	idim = image[i] & 1023;
        otherdims = image[i] ^ idim;
	idim--;
	idim &= 1023;
	image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
	x[i][0] -= period[0];
	x[i][0] = MAX(x[i][0],lo[0]);
	if (deform_remap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
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
	if (deform_remap && mask[i] & deform_groupbit) {
	  v[i][0] += h_rate[5];
	  v[i][1] += h_rate[1];
	}
	idim = (image[i] >> 10) & 1023;
        otherdims = image[i] ^ (idim << 10);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 10);
      }
      if (x[i][1] >= hi[1]) {
	x[i][1] -= period[1];
	x[i][1] = MAX(x[i][1],lo[1]);
	if (deform_remap && mask[i] & deform_groupbit) {
	  v[i][0] -= h_rate[5];
	  v[i][1] -= h_rate[1];
	}
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
	if (deform_remap && mask[i] & deform_groupbit) {
	  v[i][0] += h_rate[4];
	  v[i][1] += h_rate[3];
	  v[i][2] += h_rate[2];
	}
	idim = image[i] >> 20;
        otherdims = image[i] ^ (idim << 20);
	idim--;
	idim &= 1023;
	image[i] = otherdims | (idim << 20);
      }
      if (x[i][2] >= hi[2]) {
	x[i][2] -= period[2];
	x[i][2] = MAX(x[i][2],lo[2]);
	if (deform_remap && mask[i] & deform_groupbit) {
	  v[i][0] -= h_rate[4];
	  v[i][1] -= h_rate[3];
	  v[i][2] -= h_rate[2];
	}
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
   for triclinic, point is converted to lamda coords (0-1) before doing remap
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
   remap the point into the periodic box no matter how far away
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
------------------------------------------------------------------------- */

void Domain::remap(double *x)
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
    while (coord[0] < lo[0]) coord[0] += period[0];
    while (coord[0] >= hi[0]) coord[0] -= period[0];
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) coord[1] += period[1];
    while (coord[1] >= hi[1]) coord[1] -= period[1];
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) coord[2] += period[2];
    while (coord[2] >= hi[2]) coord[2] -= period[2];
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   x overwritten with result, don't reset image flag
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
   unmap the point via image flags
   result returned in y, don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(double *x, int image, double *y)
{
  int xbox = (image & 1023) - 512;
  int ybox = (image >> 10 & 1023) - 512;
  int zbox = (image >> 20) - 512;

  if (triclinic == 0) {
    y[0] = x[0] + xbox*xprd;
    y[1] = x[1] + ybox*yprd;
    y[2] = x[2] + zbox*zprd;
  } else {
    y[0] = x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    y[1] = x[1] + h[1]*ybox + h[3]*zbox;
    y[2] = x[2] + h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   create a lattice
   delete it if style = none
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
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

  if (find_region(arg[0]) >= 0) error->all("Reuse of region ID");

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
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int Domain::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
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

void Domain::print_box(const char *str)
{
  if (comm->me == 0) {
    if (screen) {
      if (triclinic == 0)
	fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
		str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
	char *format = (char *)
	  "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
	fprintf(screen,format,
		str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
		xy,xz,yz);
      }
    }
    if (logfile) {
      if (triclinic == 0)
	fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
		str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
	char *format = (char *)
	  "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
	fprintf(logfile,format,
		str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
		xy,xz,yz);
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

/* ----------------------------------------------------------------------
   convert 8 lamda corner pts of lo/hi box to box coords
   return bboxlo/hi = bounding box around 8 corner pts in box coords
------------------------------------------------------------------------- */

void Domain::bbox(double *lo, double *hi, double *bboxlo, double *bboxhi)
{
  double x[3];

  bboxlo[0] = bboxlo[1] = bboxlo[2] = BIG;
  bboxhi[0] = bboxhi[1] = bboxhi[2] = -BIG;

  x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
}
