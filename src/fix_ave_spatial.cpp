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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_ave_spatial.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{LOWER,CENTER,UPPER,COORD};
enum{DENSITY_MASS,DENSITY_NUM,VX,VY,VZ,FX,FY,FZ,COMPUTE};
enum{SAMPLE,ALL};
enum{BOX,LATTICE,REDUCED};

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all("Illegal fix ave/spatial command");

  nevery = atoi(arg[3]);
  nfreq = atoi(arg[4]);

  if (strcmp(arg[5],"x") == 0) dim = 0;
  else if (strcmp(arg[5],"y") == 0) dim = 1;
  else if (strcmp(arg[5],"z") == 0) dim = 2;
  else error->all("Illegal fix ave/spatial command");

  if (strcmp(arg[6],"lower") == 0) originflag = LOWER;
  if (strcmp(arg[6],"center") == 0) originflag = CENTER;
  if (strcmp(arg[6],"upper") == 0) originflag = UPPER;
  else originflag = COORD;
  if (originflag == COORD) origin = atof(arg[6]);

  delta = atof(arg[7]);

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[8],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix ave/spatial file %s",arg[8]);
      error->one(str);
    }
  }

  if (strcmp(arg[9],"density") == 0) {
    if (strcmp(arg[10],"mass") == 0) which = DENSITY_MASS;
    else if (strcmp(arg[10],"number") == 0) which = DENSITY_NUM;
    else error->all("Illegal fix ave/spatial command");

  } else if (strcmp(arg[9],"atom") == 0) {
    if (strcmp(arg[10],"vx") == 0) which = VX;
    else if (strcmp(arg[10],"vy") == 0) which = VY;
    else if (strcmp(arg[10],"vz") == 0) which = VZ;
    else if (strcmp(arg[10],"fx") == 0) which = FX;
    else if (strcmp(arg[10],"fy") == 0) which = FY;
    else if (strcmp(arg[10],"fz") == 0) which = FZ;
    else error->all("Illegal fix ave/spatial command");

  } else if (strcmp(arg[9],"compute") == 0) {
    which = COMPUTE;
    int n = strlen(arg[10]) + 1;
    id_compute = new char[n];
    strcpy(id_compute,arg[10]);

  } else error->all("Illegal fix ave/spatial command");

  // parse optional args

  normflag = ALL;
  scaleflag = BOX;

  int iarg = 11;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"all") == 0) normflag = ALL;
      else if (strcmp(arg[iarg+1],"sample") == 0) normflag = SAMPLE;
      else error->all("Illegal fix ave/spatial command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"reduced") == 0) scaleflag = REDUCED;
      else error->all("Illegal fix ave/spatial command");
      iarg += 2;
    } else error->all("Illegal fix ave/spatial command");
  }

  // if density, no normalization by atom count should be done
  // thus ALL and SAMPLE should give same answer, but code does normalize
  // thus only ALL is computed correctly, so force norm to be ALL

  if (which == DENSITY_MASS || which == DENSITY_NUM) normflag = ALL;

  // setup scaling

  int triclinic = domain->triclinic;
  if (triclinic == 1 && scaleflag != REDUCED)
    error->all("Fix ave/spatial for triclinic boxes requires units reduced");

  if (scaleflag == LATTICE && domain->lattice == NULL)
    error->all("Use of fix ave/spatial with undefined lattice");

  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors

  double scale;
  if (dim == 0) scale = xscale;
  if (dim == 1) scale = yscale;
  if (dim == 2) scale = zscale;
  delta *= scale;
  if (originflag == COORD) origin *= scale;

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/spatial command");
  if (nfreq < nevery || nfreq % nevery)
    error->all("Illegal fix ave/spatial command");

  if (delta <= 0.0) error->all("Illegal fix ave/spatial command");
  invdelta = 1.0/delta;

  nvalues = 1;
  if (which == COMPUTE) {
    int icompute = modify->find_compute(id_compute);
    if (icompute < 0)
      error->all("Compute ID for fix ave/spatial does not exist");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all("Fix ave/spatial compute does not calculate per-atom info");
    nvalues = compute_size_peratom = modify->compute[icompute]->size_peratom;
    if (nvalues == 0) nvalues = 1;
  }

  if (me == 0) {
    fprintf(fp,"Spatial-averaged data for fix %s, group %s, and %s %s\n",
	    id,group->names[igroup],arg[9],arg[10]);
    fprintf(fp,"TimeStep Number-of-layers (one per snapshot)\n");
    fprintf(fp,"Layer Coord Atoms Value(s) (one per layer)\n");
  }
  
  nsum = nlayers = maxlayer = 0;
  coord = NULL;
  count_one = count_many = count_total = NULL;
  values_one = values_many = values_total = NULL;
}

/* ---------------------------------------------------------------------- */

FixAveSpatial::~FixAveSpatial()
{
  if (which == COMPUTE) delete [] id_compute;
  if (me == 0) fclose(fp);

  memory->sfree(coord);
  memory->sfree(count_one);
  memory->sfree(count_many);
  memory->sfree(count_total);
  memory->destroy_2d_double_array(values_one);
  memory->destroy_2d_double_array(values_many);
  memory->destroy_2d_double_array(values_total);
}

/* ---------------------------------------------------------------------- */

int FixAveSpatial::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::init()
{
  // set ptrs to current compute and precompute

  if (which == COMPUTE) {
    int icompute = modify->find_compute(id_compute);
    if (icompute < 0) 
      error->all("Compute ID for fix ave/spatial does not exist");
    compute = modify->compute[icompute];
    
    if (compute->id_pre) {
      icompute = modify->find_compute(compute->id_pre);
      if (icompute < 0)
	error->all("Precompute ID for fix ave/spatial does not exist");
      precompute = modify->compute[icompute];
    } else precompute = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::end_of_step()
{
  int i,j,m,ilayer;
  double lo,hi;

  // if computing the first sample, setup layers
  // compute current origin = boundary for some layer
  // lo = layer boundary immediately below boxlo
  // hi = layer boundary immediately above boxhi
  // allocate and initialize arrays based on new layer count

  if (nsum == 0) {
    double *boxlo,*boxhi,*prd;
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }

    if (originflag == LOWER) origin = boxlo[dim];
    else if (originflag == UPPER) origin = boxhi[dim];
    else if (originflag == CENTER) origin = 0.5 * (boxlo[dim] + boxhi[dim]);

    if (origin < boxlo[dim]) {
      m = static_cast<int> ((boxlo[dim] - origin) * invdelta);
      lo = origin + m*delta;
    } else {
      m = static_cast<int> ((origin - boxlo[dim]) * invdelta);
      lo = origin - m*delta;
      if (lo > boxlo[dim]) lo -= delta;
    }
    if (origin < boxhi[dim]) {
      m = static_cast<int> ((boxhi[dim] - origin) * invdelta);
      hi = origin + m*delta;
      if (hi < boxhi[dim]) hi += delta;
    } else {
      m = static_cast<int> ((origin - boxhi[dim]) * invdelta);
      hi = origin - m*delta;
    }

    offset = lo;
    nlayers = static_cast<int> ((hi-lo) * invdelta + 0.5);
    double volume = domain->xprd * domain->yprd * domain->zprd;
    layer_volume = delta * volume/prd[dim];

    if (nlayers > maxlayer) {
      maxlayer = nlayers;
      coord = (double *) memory->srealloc(coord,nlayers*sizeof(double),
					  "ave/spatial:coord");
      count_one = (double *) 
	memory->srealloc(count_one,nlayers*sizeof(double),
			 "ave/spatial:count_one");
      count_many = (double *) 
	memory->srealloc(count_many,nlayers*sizeof(double),
			 "ave/spatial:count_many");
      count_total = (double *) 
	memory->srealloc(count_total,nlayers*sizeof(double),
			 "ave/spatial:count_total");
      values_one = memory->grow_2d_double_array(values_one,nlayers,nvalues,
						"ave/spatial:values_one");
      values_many = memory->grow_2d_double_array(values_many,nlayers,nvalues,
						 "ave/spatial:values_many");
      values_total = memory->grow_2d_double_array(values_total,nlayers,nvalues,
						"ave/spatial:values_total");
    }

    for (m = 0; m < nlayers; m++) {
      coord[m] = offset + (m+0.5)*delta;
      count_many[m] = count_total[m] = 0.0;
      for (i = 0; i < nvalues; i++) values_many[m][i] = 0.0;
    }
  }
  
  // zero out arrays for one sample

  nsum++;
  for (m = 0; m < nlayers; m++) {
    count_one[m] = 0.0;
    for (i = 0; i < nvalues; i++) values_one[m][i] = 0.0;
  }

  // perform the computation for one sample
  // sum within each layer, only include atoms in fix group
  // insure array index is within bounds (since atoms can be outside box)
  // if scaleflag = REDUCED, convert box coords to lamda coords
  // DENSITY_MASS adds mass to values
  // DENSITY_NUM adds 1 to values
  // ATOM adds atom vector to values
  // COMPUTE adds its vector to values

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  if (which == DENSITY_MASS) {
    int *type = atom->type;
    double *mass = atom->mass;
    double *rmass = atom->rmass;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	if (mass) values_one[ilayer][0] += mass[type[i]];
	else values_one[ilayer][0] += rmass[i];
      }
    }

  } else if (which == DENSITY_NUM) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	values_one[ilayer][0] += 1.0;
      }
    }

  } else if (which != COMPUTE) {
    double *vector;
    int nstride = 3;
    if (which == VX) vector = &atom->v[0][0];
    else if (which == VY) vector = &atom->v[0][1];
    else if (which == VZ) vector = &atom->v[0][2];
    else if (which == FX) vector = &atom->f[0][0];
    else if (which == FY) vector = &atom->f[0][1];
    else if (which == FZ) vector = &atom->f[0][2];

    m = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	values_one[ilayer][0] += vector[m];
      }
      m += nstride;
    }

  } else {
    if (precompute) precompute->compute_peratom();
    compute->compute_peratom();
    double *scalar = compute->scalar_atom;
    double **vector = compute->vector_atom;

    m = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	if (compute_size_peratom == 0) values_one[ilayer][0] += scalar[i];
	else
	  for (j = 0; j < nvalues; j++)
	    values_one[ilayer][j] += vector[i][j];
      }
    }
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  // average a single sample

  if (normflag == ALL) {
    for (m = 0; m < nlayers; m++) {
      count_many[m] += count_one[m];
      for (j = 0; j < nvalues; j++)
	values_many[m][j] += values_one[m][j];
    }
  } else {
    MPI_Allreduce(count_one,count_many,nlayers,MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nlayers; m++) {
      if (count_many[m] > 0.0)
	for (j = 0; j < nvalues; j++)
	  values_many[m][j] += values_one[m][j]/count_many[m];
      count_total[m] += count_many[m];
    }
  }

  // output the results
  // time average across samples
  // if density, also normalize by volume

  if (update->ntimestep % nfreq == 0) {
    if (normflag == ALL) {
      MPI_Allreduce(count_many,count_total,nlayers,MPI_DOUBLE,MPI_SUM,world);
      MPI_Allreduce(&values_many[0][0],&values_total[0][0],nlayers*nvalues,
		    MPI_DOUBLE,MPI_SUM,world);
      for (m = 0; m < nlayers; m++) {
	if (count_total[m] > 0.0)
	  for (j = 0; j < nvalues; j++)
	    values_total[m][j] /= count_total[m];
	count_total[m] /= nsum;
      }
    } else {
      MPI_Allreduce(&values_many[0][0],&values_total[0][0],nlayers*nvalues,
		    MPI_DOUBLE,MPI_SUM,world);
      for (m = 0; m < nlayers; m++) {
	for (j = 0; j < nvalues; j++)
	  values_total[m][j] /= nsum;
	count_total[m] /= nsum;
      }
    }

    if (which == DENSITY_MASS || which == DENSITY_NUM) {
      for (m = 0; m < nlayers; m++)
	values_total[m][0] *= count_total[m] / layer_volume;
    }

    if (me == 0) {
      fprintf(fp,"%d %d\n",update->ntimestep,nlayers);
      for (m = 0; m < nlayers; m++) {
	fprintf(fp,"  %d %g %g",m+1,coord[m],count_total[m]);
	for (i = 0; i < nvalues; i++) fprintf(fp," %g",values_total[m][i]);
	fprintf(fp,"\n");
      }
      fflush(fp);
    }

    nsum = 0;
 }
}
