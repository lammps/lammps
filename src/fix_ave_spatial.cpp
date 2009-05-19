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
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{LOWER,CENTER,UPPER,COORD};
enum{X,V,F,DENSITY_NUMBER,DENSITY_MASS,COMPUTE,FIX,VARIABLE};
enum{SAMPLE,ALL};
enum{BOX,LATTICE,REDUCED};
enum{ONE,RUNNING,WINDOW};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

#define BIG 1000000000

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all("Illegal fix ave/spatial command");

  MPI_Comm_rank(world,&me);

  no_change_box = 1;
  time_depend = 1;

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
  if (originflag == COORD) origin = atof(arg[7]);

  delta = atof(arg[8]);

  // parse values until one isn't recognized

  which = new int[narg-9];
  argindex = new int[narg-9];
  ids = new char*[narg-9];
  value2index = new int[narg-9];
  nvalues = 0;

  int iarg = 9;
  while (iarg < narg) {
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"x") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"y") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"z") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"density/number") == 0) {
      which[nvalues] = DENSITY_NUMBER;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"density/mass") == 0) {
      which[nvalues] = DENSITY_MASS;
      argindex[nvalues++] = 0;

    } else if ((strncmp(arg[iarg],"c_",2) == 0) || 
	       (strncmp(arg[iarg],"f_",2) == 0) || 
	       (strncmp(arg[iarg],"v_",2) == 0)) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Illegal fix ave/spatial command");
	argindex[nvalues] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // optional args

  normflag = ALL;
  scaleflag = LATTICE;
  fp = NULL;
  ave = ONE;

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

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/spatial command");
  if (nfreq < nevery || nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all("Illegal fix ave/spatial command");

  if (delta <= 0.0) error->all("Illegal fix ave/spatial command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/spatial does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
	error->all("Fix ave/spatial compute does not calculate per-atom values");
      if (argindex[i] == 0 && modify->compute[icompute]->size_peratom != 0)
	error->all("Fix ave/spatial compute does not calculate a per-atom scalar");
      if (argindex[i] && modify->compute[icompute]->size_peratom == 0)
	error->all("Fix ave/spatial compute does not calculate a per-atom vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_peratom)
	error->all("Fix ave/spatial compute vector is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all("Fix ID for fix ave/spatial does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
	error->all("Fix ave/spatial fix does not calculate per-atom values");
      if (argindex[i] && modify->fix[ifix]->size_peratom != 0)
	error->all("Fix ave/spatial fix does not calculate a per-atom scalar");
      if (argindex[i] && modify->fix[ifix]->size_peratom == 0)
	error->all("Fix ave/spatial fix does not calculate a per-atom vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom)
	error->all("Fix ave/spatial fix vector is accessed out-of-range");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all("Variable name for fix ave/spatial does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
	error->all("Fix ave/spatial variable is not atom-style variable");
    }
  }

  // print header into file

  if (fp && me == 0) {
    fprintf(fp,"# Spatial-averaged data for fix %s and group %s\n",id,arg[1]);
    fprintf(fp,"# TimeStep Number-of-layers\n");
    fprintf(fp,"# Layer Coordinate Natoms");
    for (int i = 0; i < nvalues; i++)
      if (which[i] == COMPUTE) fprintf(fp," c_%s",ids[i]);
      else if (which[i] == FIX) fprintf(fp," f_%s",ids[i]);
      else if (which[i] == VARIABLE) fprintf(fp," v_%s",ids[i]);
      else fprintf(fp," %s",arg[9+i]);
    fprintf(fp,"\n");
  }
  
  // this fix produces a global vector
  // set size_vector to BIG since compute_vector() checks bounds on-the-fly

  vector_flag = 1;
  size_vector = BIG;
  scalar_vector_freq = nfreq;
  extvector = 0;

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

  invdelta = 1.0/delta;

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

  maxatomvar = 0;
  varatom = NULL;

  maxatomlayer = 0;
  layer = NULL;

  // nvalid = next step on which end_of_step does something
  // can be this timestep if multiple of nfreq and nrepeat = 1
  // else backup from next multiple of nfreq

  nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;

  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveSpatial::~FixAveSpatial()
{
  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;

  if (fp && me == 0) fclose(fp);

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

  memory->sfree(varatom);
  memory->sfree(layer);
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

  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
	error->all("Compute ID for fix ave/spatial does not exist");
      value2index[m] = icompute;
      
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) 
	error->all("Fix ID for fix ave/spatial does not exist");
      value2index[m] = ifix;

      if (nevery % modify->fix[ifix]->peratom_freq)
	error->all("Fix for fix ave/spatial not computed at compatible time");

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) 
	error->all("Variable name for fix ave/spatial does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSpatial::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::end_of_step()
{
  int i,j,m,n,ilayer;
  double lo,hi;

  // skip if not step which requires doing something

  int ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

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

  // assign each atom to a layer
  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;


  if (nlocal > maxatomlayer) {
    maxatomlayer = atom->nmax;
    memory->sfree(layer);
    layer = (int *) 
      memory->smalloc(maxatomlayer*sizeof(int),"ave/spatial:layer");
  }

  double *boxlo,*boxhi,*prd;
  double xremap;
  int periodicity = domain->periodicity[dim];

  if (periodicity) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);
  
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xremap = x[i][dim];
      if (periodicity) {
	if (xremap < boxlo[dim]) xremap += prd[dim];
	if (xremap >= boxhi[dim]) xremap -= prd[dim];
      }
      ilayer = static_cast<int> ((xremap - offset) * invdelta);
      if (ilayer < 0) ilayer = 0;
      if (ilayer >= nlayers) ilayer = nlayers-1;
      layer[i] = ilayer;
      count_one[ilayer] += 1.0;
    }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  // perform the computation for one sample
  // accumulate results of attributes,computes,fixes,variables to local copy
  // sum within each layer, only include atoms in fix group
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    // X,V,F adds coords,velocities,forces to values

    if (which[m] == X || which[m] == V || which[m] == F) {
      double **attribute;
      if (which[m] == X) attribute = x;
      else if (which[m] == V) attribute = atom->v;
      else attribute = atom->f;

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  values_one[layer[i]][m] += attribute[i][j];

    // DENSITY_NUMBER adds 1 to values

    } else if (which[m] == DENSITY_NUMBER) {

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  values_one[layer[i]][m] += 1.0;

    // DENSITY_MASS adds mass to values

    } else if (which[m] == DENSITY_MASS) {
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  if (rmass) values_one[layer[i]][m] += rmass[i];
	  else values_one[layer[i]][m] += mass[type[i]];

    // COMPUTE adds its scalar or vector component to values
    // invoke compute if not previously invoked

    } else if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	compute->compute_peratom();
	compute->invoked_flag |= INVOKED_PERATOM;
      }
      double *scalar = compute->scalar_atom;
      double **vector = compute->vector_atom;
      int jm1 = j - 1;

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  if (j == 0) values_one[layer[i]][m] += scalar[i];
	  else values_one[layer[i]][m] += vector[i][jm1];
      
    // FIX adds its scalar or vector component to values
    // access fix fields, guaranteed to be ready

    } else if (which[m] == FIX) {
      double *scalar = modify->fix[n]->scalar_atom;
      double **vector = modify->fix[n]->vector_atom;
      int jm1 = j - 1;

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  if (j == 0) values_one[layer[i]][m] += scalar[i];
	  else values_one[layer[i]][m] += vector[i][jm1];
	}

    // VARIABLE adds its per-atom quantities to values
    // evaluate atom-style variable

    } else if (which[m] == VARIABLE) {
      if (nlocal > maxatomvar) {
	maxatomvar = atom->nmax;
	memory->sfree(varatom);
	varatom = (double *) 
	  memory->smalloc(maxatomvar*sizeof(double),"ave/spatial:varatom");
      }

      input->variable->compute_atom(n,igroup,varatom,1,0);

      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  values_one[layer[i]][m] += varatom[i];
    }
  }

  // process a single sample
  // if normflag = ALL, accumulate values,count separately to many
  // if normflag = SAMPLE, one = value/count, accumulate one to many
  // exception is SAMPLE density: no normalization by atom count

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
	for (j = 0; j < nvalues; j++) {
	  if (which[j] == DENSITY_NUMBER || which[j] == DENSITY_MASS)
	    values_many[m][j] += values_one[m][j];
	  else values_many[m][j] += values_one[m][j]/count_many[m];
	}
      count_sum[m] += count_many[m];
    }
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // time average across samples
  // if normflag = ALL, final is total value / total count
  // if normflag = SAMPLE, final is sum of ave / repeat
  // exception is ALL density: normalized by repeat, not total count

  double repeat = nrepeat;

  if (normflag == ALL) {
    MPI_Allreduce(count_many,count_sum,nlayers,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nlayers*nvalues,
		  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nlayers; m++) {
      if (count_sum[m] > 0.0)
	for (j = 0; j < nvalues; j++)
	  if (which[j] == DENSITY_NUMBER || which[j] == DENSITY_MASS)
	    values_sum[m][j] /= repeat;
	  else values_sum[m][j] /= count_sum[m];
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

  // density is additionally normalized by layer volume

  for (j = 0; j < nvalues; j++)
    if (which[j] == DENSITY_NUMBER || which[j] == DENSITY_MASS)
      for (m = 0; m < nlayers; m++)
	values_sum[m][j] /= layer_volume;

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
    fprintf(fp,"%d %d\n",ntimestep,nlayers);
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
   if ilayer exceeds current layers, return 0.0 instead of generate an error
------------------------------------------------------------------------- */

double FixAveSpatial::compute_vector(int n)
{
  int ivalue = n % nvalues;
  int ilayer = n / nvalues;
  if (ilayer >= nlayers) return 0.0;
  if (values_total == NULL) return 0.0;
  return values_total[ilayer][ivalue]/norm;
}

/* ----------------------------------------------------------------------
   memory usage of varatom and layer
------------------------------------------------------------------------- */

double FixAveSpatial::memory_usage()
{
  double bytes = maxatomvar * sizeof(double);
  bytes += maxatomlayer * sizeof(int);
  return bytes;
}
