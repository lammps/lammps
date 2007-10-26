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
enum{DENSITY_MASS,DENSITY_NUM,COMPUTE,FIX};
enum{SAMPLE,ALL};
enum{BOX,LATTICE,REDUCED};
enum{ONE,RUNNING,WINDOW};

#define BIG 1000000000

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 11) error->all("Illegal fix ave/spatial command");

  MPI_Comm_rank(world,&me);

  no_change_box = 1;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  nfreq = atoi(arg[5]);

  if (strcmp(arg[6],"x") == 0) dim = 0;
  else if (strcmp(arg[6],"y") == 0) dim = 1;
  else if (strcmp(arg[6],"z") == 0) dim = 2;
  else error->all("Illegal fix ave/spatial command");

  if (strcmp(arg[7],"lower") == 0) originflag = LOWER;
  if (strcmp(arg[7],"center") == 0) originflag = CENTER;
  if (strcmp(arg[7],"upper") == 0) originflag = UPPER;
  else originflag = COORD;
  if (originflag == COORD) origin = atof(arg[6]);

  delta = atof(arg[8]);

  if (strcmp(arg[9],"density") == 0) {
    if (strcmp(arg[10],"mass") == 0) which = DENSITY_MASS;
    else if (strcmp(arg[10],"number") == 0) which = DENSITY_NUM;
    else error->all("Illegal fix ave/spatial command");
  } else if (strcmp(arg[9],"compute") == 0) {
    which = COMPUTE;
    int n = strlen(arg[10]) + 1;
    id_compute = new char[n];
    strcpy(id_compute,arg[10]);
  } else if (strcmp(arg[9],"fix") == 0) {
    which = FIX;
    int n = strlen(arg[10]) + 1;
    id_fix = new char[n];
    strcpy(id_fix,arg[10]);
  } else error->all("Illegal fix ave/spatial command");

  // parse optional args

  normflag = ALL;
  scaleflag = BOX;
  fp = NULL;
  ave = ONE;

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
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      if (me == 0) {
	fp = fopen(arg[iarg+1],"w");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix ave/spatial file %s",arg[iarg+1]);
	  error->one(str);
	}
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all("Illegal fix ave/spatial command");
      if (ave == WINDOW) {
	if (iarg+3 > narg) error->all("Illegal fix ave/spatial command");
	nwindow = atoi(arg[iarg+2]);
	if (nwindow <= 0) error->all("Illegal fix ave/spatial command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
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
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/spatial command");

  if (delta <= 0.0) error->all("Illegal fix ave/spatial command");
  invdelta = 1.0/delta;

  // nvalues = # of quantites per line of output file
  // for COMPUTE, setup list of computes to call, including pre-computes

  nvalues = 1;
  compute = NULL;

  if (which == COMPUTE) {
    int icompute = modify->find_compute(id_compute);
    if (icompute < 0)
      error->all("Compute ID for fix ave/spatial does not exist");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all("Fix ave/spatial compute does not calculate per-atom info");
    nvalues = size_peratom = modify->compute[icompute]->size_peratom;
    if (nvalues == 0) nvalues = 1;
    ncompute = 1 + modify->compute[icompute]->npre;
    compute = new Compute*[ncompute];
  }

  if (which == FIX) {
    int ifix = modify->find_fix(id_fix);
    if (ifix < 0)
      error->all("Fix ID for fix ave/spatial does not exist");
    if (modify->fix[ifix]->peratom_flag == 0)
      error->all("Fix ave/spatial fix does not calculate per-atom info");
    nvalues = size_peratom = modify->fix[ifix]->size_peratom;
    if (nvalues == 0) nvalues = 1;
  }

  // print header into file

  if (fp && me == 0) {
    fprintf(fp,"Spatial-averaged data for fix %s, group %s, and %s %s\n",
	    id,group->names[igroup],arg[10],arg[11]);
    fprintf(fp,"TimeStep Number-of-layers (one per snapshot)\n");
    fprintf(fp,"Layer Coord Atoms Value(s) (one per layer)\n");
  }
  
  // enable this fix to produce a global vector
  // set size_vector to BIG since compute_vector() will check bounds

  vector_flag = 1;
  size_vector = BIG;
  scalar_vector_freq = nfreq;
  extensive = 0;

  // initializations

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;
  nlayers = maxlayer = 0;
  coord = NULL;
  count_one = count_many = count_sum = count_total = NULL;
  count_list = NULL;
  values_one = values_many = values_sum = values_total = NULL;
  values_list = NULL;

  // nvalid = next step on which end_of_step does something
  // can be this timestep if multiple of nfreq and nrepeat = 1
  // else backup from next multiple of nfreq

  nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;

  if (nvalid < update->ntimestep)
    error->all("Fix ave/spatial cannot be started on this timestep");
}

/* ---------------------------------------------------------------------- */

FixAveSpatial::~FixAveSpatial()
{
  if (which == COMPUTE) delete [] id_compute;
  if (which == FIX) delete [] id_fix;
  if (fp && me == 0) fclose(fp);

  delete [] compute;

  memory->sfree(coord);
  memory->sfree(count_one);
  memory->sfree(count_many);
  memory->sfree(count_sum);
  memory->sfree(count_total);
  memory->destroy_2d_double_array(count_list);
  memory->destroy_2d_double_array(values_one);
  memory->destroy_2d_double_array(values_many);
  memory->destroy_2d_double_array(values_sum);
  memory->destroy_2d_double_array(values_total);
  memory->destroy_3d_double_array(values_list);
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
  // # of layers cannot vary for ave = RUNNING or WINDOW

  if (ave == RUNNING || ave == WINDOW) {
    if (scaleflag != REDUCED && domain->box_change)
      error->all("Fix ave/spatial settings invalid with changing box");
  }

  // set ptrs to compute and its pre-computes called each end-of-step
  // put pre-computes in list before compute

  if (which == COMPUTE) {
    int icompute = modify->find_compute(id_compute);
    if (icompute < 0)
      error->all("Compute ID for fix ave/spatial does not exist");
    
    ncompute = 0;
    if (modify->compute[icompute]->npre)
      for (int i = 0; i < modify->compute[icompute]->npre; i++) {
	int ic = modify->find_compute(modify->compute[icompute]->id_pre[i]);
	if (ic < 0)
	  error->all("Precompute ID for fix ave/spatial does not exist");
	compute[ncompute++] = modify->compute[ic];
      }
    
    compute[ncompute++] = modify->compute[icompute];
  }

  // set ptr to fix ID
  // check that fix frequency is acceptable

  if (which == FIX) {
    int ifix = modify->find_fix(id_fix);
    if (ifix < 0) 
      error->all("Fix ID for fix ave/spatial does not exist");
    fix = modify->fix[ifix];
    if (nevery % fix->peratom_freq)
      error->all("Fix ave/spatial and fix not computed at compatible times");
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSpatial::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::end_of_step()
{
  int i,j,m,ilayer;
  double lo,hi;

  // skip if not step which requires doing something

  if (update->ntimestep != nvalid) return;

  // if computing the first sample, setup layers
  // compute current origin = boundary for some layer
  // lo = layer boundary immediately below boxlo
  // hi = layer boundary immediately above boxhi
  // allocate and initialize arrays based on new layer count

  if (irepeat == 0) {
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
      count_sum = (double *) 
	memory->srealloc(count_sum,nlayers*sizeof(double),
			 "ave/spatial:count_sum");
      count_total = (double *) 
	memory->srealloc(count_total,nlayers*sizeof(double),
			 "ave/spatial:count_total");

      values_one = memory->grow_2d_double_array(values_one,nlayers,nvalues,
						"ave/spatial:values_one");
      values_many = memory->grow_2d_double_array(values_many,nlayers,nvalues,
						 "ave/spatial:values_many");
      values_sum = memory->grow_2d_double_array(values_sum,nlayers,nvalues,
						"ave/spatial:values_sum");
      values_total = memory->grow_2d_double_array(values_total,nlayers,nvalues,
						  "ave/spatial:values_total");

      // initialize count and values total to zero since they accumulate

      for (m = 0; m < nlayers; m++) {
	for (i = 0; i < nvalues; i++) values_total[m][i] = 0.0;
	count_total[m] = 0.0;
      }

      // only allocate count and values list for ave = WINDOW
      // only happens once since nlayers never changes for these ave settings
      
      if (ave == WINDOW) {
	count_list =
	  memory->create_2d_double_array(nwindow,nlayers,
					 "ave/spatial:count_list");
	values_list =
	  memory->create_3d_double_array(nwindow,nlayers,nvalues,
					 "ave/spatial:values_list");
      }
    }

    for (m = 0; m < nlayers; m++) {
      coord[m] = offset + (m+0.5)*delta;
      count_many[m] = count_sum[m] = 0.0;
      for (i = 0; i < nvalues; i++) values_many[m][i] = 0.0;
    }
  }
  
  // zero out arrays for one sample

  for (m = 0; m < nlayers; m++) {
    count_one[m] = 0.0;
    for (i = 0; i < nvalues; i++) values_one[m][i] = 0.0;
  }

  // perform the computation for one sample
  // sum within each layer, only include atoms in fix group
  // insure array index is within bounds (since atoms can be outside box)
  // if scaleflag = REDUCED, box coords -> lamda coords before computing layer

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // DENSITY_MASS adds mass to values

  if (which == DENSITY_MASS) {
    int *type = atom->type;
    double *mass = atom->mass;
    double *rmass = atom->rmass;

    if (scaleflag == REDUCED) domain->x2lamda(nlocal);

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

    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  // DENSITY_NUM adds 1 to values

  } else if (which == DENSITY_NUM) {

    if (scaleflag == REDUCED) domain->x2lamda(nlocal);

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	values_one[ilayer][0] += 1.0;
      }
    }

    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  // COMPUTE adds its scalar or vector quantity to values

  } else if (which == COMPUTE) {
    for (i = 0; i < ncompute; i++) compute[i]->compute_peratom();
    double *scalar = compute[ncompute-1]->scalar_atom;
    double **vector = compute[ncompute-1]->vector_atom;

    if (scaleflag == REDUCED) domain->x2lamda(nlocal);

    m = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	if (size_peratom == 0) values_one[ilayer][0] += scalar[i];
	else
	  for (j = 0; j < nvalues; j++)
	    values_one[ilayer][j] += vector[i][j];
      }
    }

    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  // FIX adds its scalar or vector quantity to values

  } else if (which == FIX) {
    double *scalar = fix->scalar_atom;
    double **vector = fix->vector_atom;

    if (scaleflag == REDUCED) domain->x2lamda(nlocal);

    m = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ilayer = static_cast<int> ((x[i][dim] - offset) * invdelta);
	if (ilayer < 0) ilayer = 0;
	if (ilayer >= nlayers) ilayer = nlayers-1;
	count_one[ilayer] += 1.0;
	if (size_peratom == 0) values_one[ilayer][0] += scalar[i];
	else
	  for (j = 0; j < nvalues; j++)
	    values_one[ilayer][j] += vector[i][j];
      }
    }

    if (scaleflag == REDUCED) domain->lamda2x(nlocal);
  }

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
      count_sum[m] += count_many[m];
    }
  }

  // done if irepeat < nrepeat

  irepeat++;
  nvalid += nevery;
  if (irepeat < nrepeat) return;

  // time average across samples
  // if density, also normalize by volume

  double repeat = nrepeat;

  if (normflag == ALL) {
    MPI_Allreduce(count_many,count_sum,nlayers,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nlayers*nvalues,
		  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nlayers; m++) {
      if (count_sum[m] > 0.0)
	for (j = 0; j < nvalues; j++)
	  values_sum[m][j] /= count_sum[m];
      count_sum[m] /= repeat;
    }
  } else {
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nlayers*nvalues,
		  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nlayers; m++) {
      for (j = 0; j < nvalues; j++)
	values_sum[m][j] /= repeat;
      count_sum[m] /= repeat;
    }
  }

  if (which == DENSITY_MASS || which == DENSITY_NUM) {
    for (m = 0; m < nlayers; m++)
      values_sum[m][0] *= count_sum[m] / layer_volume;
  }

  // reset irepeat and nvalid

  irepeat = 0;
  nvalid = update->ntimestep+nfreq - (nrepeat-1)*nevery;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, comine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (m = 0; m < nlayers; m++) {
      for (i = 0; i < nvalues; i++) 
	values_total[m][i] = values_sum[m][i];
      count_total[m] = count_sum[m];
    }
    norm = 1;

  } else if (ave == RUNNING) {
    for (m = 0; m < nlayers; m++) {
      for (i = 0; i < nvalues; i++)
	values_total[m][i] += values_sum[m][i];
      count_total[m] += count_sum[m];
    }
    norm++;

  } else if (ave == WINDOW) {
    for (m = 0; m < nlayers; m++) {
      for (i = 0; i < nvalues; i++) {
	values_total[m][i] += values_sum[m][i];
	if (window_limit) values_total[m][i] -= values_list[iwindow][m][i];
	values_list[iwindow][m][i] = values_sum[m][i];
      }
      count_total[m] += count_sum[m];
      if (window_limit) count_total[m] -= count_list[iwindow][m];
      count_list[iwindow][m] = count_sum[m];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // output result to file
  
  if (fp && me == 0) {
    fprintf(fp,"%d %d\n",update->ntimestep,nlayers);
    for (m = 0; m < nlayers; m++) {
      fprintf(fp,"  %d %g %g",m+1,coord[m],count_total[m]/norm);
      for (i = 0; i < nvalues; i++) fprintf(fp," %g",values_total[m][i]/norm);
      fprintf(fp,"\n");
    }
    fflush(fp);
  }
}

/* ----------------------------------------------------------------------
   return Nth vector value
   since values_sum is 2d array, map N into ilayer and ivalue
   if ilayer >= nlayers, just return 0, since nlayers can vary with time
------------------------------------------------------------------------- */

double FixAveSpatial::compute_vector(int n)
{
  int ivalue = n % nvalues;
  int ilayer = n / nvalues;
  if (ilayer < nlayers && norm) return values_total[ilayer][ivalue]/norm;
  return 0.0;
}
