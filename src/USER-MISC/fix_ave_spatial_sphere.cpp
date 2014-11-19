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
   Adapted from fix ave/spatial by Niall Jackson <niall.jackson@gmail.com>
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "fix_ave_spatial_sphere.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

//the different types of variables that we can track
enum{V,F,DENSITY_NUMBER,DENSITY_MASS,COMPUTE,FIX,VARIABLE,CONSTANT};
//normalisation types
enum{SAMPLE,ALL};
//unit scaling options
enum{BOX,LATTICE,REDUCED};
//averaging types
enum{ONE,RUNNING,WINDOW};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_PERATOM 8
#define BIG 1000000000

/* ---------------------------------------------------------------------- */

FixAveSpatialSphere::FixAveSpatialSphere(LAMMPS* lmp, int narg, char** arg) :
  Fix(lmp, narg, arg) {
  
  if (narg < 12) error->all(FLERR, "Illegal fix ave/spatial/sphere command");
  
  MPI_Comm_rank(world, &me);
  
  nevery= force->inumeric(FLERR, arg[3]);
  nrepeat= force->inumeric(FLERR, arg[4]);
  nfreq= force->inumeric(FLERR, arg[5]);
  
  global_freq= nfreq;
  
  for (int iarg = 6; iarg < 9; ++iarg) {
    int dim= iarg-6;
    origin_ids[dim]= NULL;
    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"v_",2) == 0) {
      if(arg[iarg][0] == 'c') {
        origin_type[dim]= COMPUTE;
      } else {
        origin_type[dim]= VARIABLE;
      }

      //strip the v_/c_ off the front of the variable name
      int n= strlen(arg[iarg]);
      char *suffix= new char[n];
      strcpy(suffix, &arg[iarg][2]);
      //does the variable name contain a [?
      //(in other words, is this an array element?)
      char *ptr= strchr(suffix, '[');
      if(ptr) {
        //if the last character ISN'T ], then the array syntax
        //has been used wrongly
        if(suffix[strlen(suffix)-1] != ']') {
          error->all(FLERR, "Illegal fix ave/spatial/sphere command");
        }
        origin_index[dim]= atoi(ptr+1);
        *ptr= '\0';
      } else {
        origin_index[dim]= 0;
      }
      //store the name of the variable
      n= strlen(suffix)+1;
      origin_ids[dim]= new char[n];
      strcpy(origin_ids[dim], suffix);
      //tidy up the array we allocated earlier
      delete[] suffix;

      if(origin_type[dim] == COMPUTE) {
        int icompute= modify->find_compute(origin_ids[dim]);
        if(icompute < 0) 
          error->all(FLERR,
                     "Compute ID for fix ave/spatial/sphere does not exist");

        if(origin_index[dim] == 0 && 
           modify->compute[icompute]->scalar_flag != 0)
          error->all(FLERR,"Compute for fix ave/spatial/sphere does not calculate a scalar");
        else if(origin_index[dim] && modify->compute[icompute]->vector_flag == 0)
          error->all(FLERR, "Compute for fix ave/spatial/sphere does not calculate a vector");
        else if(origin_index[dim] && origin_index[dim] > modify->compute[icompute]->size_vector)
          error->all(FLERR, "Fix ave/spatial/sphere compute vector accessed out of range");
      } else if(origin_type[dim] == VARIABLE) {
        int ivariable = input->variable->find(origin_ids[dim]);
        if (ivariable < 0)
          error->all(FLERR,"Variable name for fix ave/spatial/sphere does not exist");
        else if (input->variable->equalstyle(ivariable) == 0)
          error->all(FLERR,"Fix ave/spatial/sphere variable is not equal-style variable");
      }
    } else {
      origin_type[dim]= CONSTANT;
      origin[dim]= force->numeric(FLERR, arg[iarg]);
    }
  }
  
  //the extrema of the radius
  r_min= force->numeric(FLERR, arg[9]);
  r_max= force->numeric(FLERR, arg[10]);
 
  //and finally, the number of spherical bins
  nbins= force->inumeric(FLERR, arg[11]);
  
  if (r_min >= r_max) {
      error->all(FLERR, "r_min must be less than r_max in fix ave/spatial/sphere");
  }

  if (nbins <= 0) {
      error->all(FLERR, "nbins must be positive in fix ave/spatial/sphere");
  }
  
  int iarg= 12;
  which= new int[narg-12];
  argindex= new int[narg-12];
  ids= new char*[narg-12];
  value2index= new int[narg-12];
  nvalues= 0;
  
  while(iarg < narg) {
    ids[nvalues]= NULL;
    
    if(strcmp(arg[iarg], "vx") == 0) {
        which[nvalues] = V;
        argindex[nvalues++]= 0;
    } else if(strcmp(arg[iarg], "vy") == 0) {
        which[nvalues] = V;
        argindex[nvalues++]= 1;
    } else if(strcmp(arg[iarg], "vz") == 0) {
        which[nvalues] = V;
        argindex[nvalues++]= 2;
    } else if(strcmp(arg[iarg], "fx") == 0) {
        which[nvalues] = F;
        argindex[nvalues++]= 0;
    } else if(strcmp(arg[iarg], "fy") == 0) {
        which[nvalues] = F;
        argindex[nvalues++]= 1;
    } else if(strcmp(arg[iarg], "fz") == 0) {
        which[nvalues] = F;
        argindex[nvalues++]= 2;
    } else if(strcmp(arg[iarg], "density/number") == 0) {
        which[nvalues] = DENSITY_NUMBER;
        argindex[nvalues++]= 0;
    } else if(strcmp(arg[iarg], "density/mass") == 0) {
        which[nvalues] = DENSITY_MASS;
        argindex[nvalues++]= 0;
    } else if (strncmp(arg[iarg],"c_",2) == 0 ||
              strncmp(arg[iarg],"f_",2) == 0 ||
              strncmp(arg[iarg],"v_",2) == 0) {
      if(arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if(arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if(arg[iarg][0] == 'v') which[nvalues] = VARIABLE;
      
      //strip the v_/c_/f_ off the front of the variable name
      int n= strlen(arg[iarg]);
      char *suffix= new char[n];
      strcpy(suffix, &arg[iarg][2]);
    
      //does the variable name contain a [?
      //(in other words, is this an array element?)
      char *ptr= strchr(suffix, '[');
      if(ptr) {
        //if the last character ISN'T ], then the array syntax
        //has been used wrongly
        if(suffix[strlen(suffix)-1] != ']') {
          error->all(FLERR, "Illegal fix ave/spatial/sphere command");
        }
        argindex[nvalues]= atoi(ptr+1);
        *ptr= '\0';
      } else {
        argindex[nvalues]= 0;
      }
      
      //store the name of the variable
      n= strlen(suffix)+1;
      ids[nvalues]= new char[n];
      strcpy(ids[nvalues], suffix);
      nvalues++;
      //tidy up the array we allocated earlier
      delete[] suffix;
    } else break;
    
    iarg++;
  }
  
  //process the optional arguments
  normflag= ALL; //normalise the average right at the end
  scaleflag= LATTICE; // lattice units are the default
  regionflag= 0; //are we restricted to a particular region?
  idregion= NULL; //name of the region to sue, if so
  fp= NULL; //output file
  ave= ONE; //averaging mode
  nwindow= 0; //number of averaging windows
  overwrite= 0; //continuously overwrite output?
  char *title1= NULL;
  char *title2= NULL;
  char *title3= NULL;
  
  while (iarg < narg) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      if (strcmp(arg[iarg+1],"all") == 0) normflag = ALL;
      else if (strcmp(arg[iarg+1],"sample") == 0) normflag = SAMPLE;
      else error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"reduced") == 0) scaleflag = REDUCED;
      else error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix ave/spatial does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix ave/spatial file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial/sphere command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/spatial/sphere command");
  }
  
  //setup and error checking
  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/spatial/sphere command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix ave/spatial/sphere command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/spatial/sphere command");
  
  //setup the variables that we need to average
  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/spatial/sphere does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,"Fix ave/spatial/sphere compute does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/spatial/sphere compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/spatial/sphere compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,
                   "Fix ave/spatial/sphere compute vector is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for  fix ave/spatial/sphere does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix ave/spatial/sphere fix does not calculate per-atom values");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,
                   "Fix ave/spatial/spherefix does not calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Fix ave/spatial/sphere fix does not calculate a per-atom array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,"Fix ave/spatial/sphere fix vector is accessed out-of-range");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/spatial/sphere does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/spatial/sphere variable is not atom-style variable");
    }
  }

  // setup the scaling
  int triclinic= domain->triclinic;
  if(scaleflag != REDUCED) {
    if (triclinic) {
      error->all(FLERR, "Fix ave/spatial/sphere for triclinic boxes requires units reduced");
    } else {
      // if not in reduced units, the box must not become triclinic
      no_change_box= 1;
    }
  }
  
  if(scaleflag == LATTICE) {
    scale= domain->lattice->xlattice;
    if(domain->lattice->xlattice != domain->lattice->ylattice
      && domain->lattice->xlattice != domain->lattice->zlattice)
      error->all(FLERR, "Fix ave/spatial/sphere with lattice units requires that the lattice "
                        "spacing be equal in all dimensions");
  } else if(scaleflag == REDUCED) {
    scale= 1.0;
  } else {
    scale= 1.0;
  }
  
  //apply the scaling factors
  origin[0]*= scale;
  origin[1]*= scale;
  origin[2]*= scale;
  r_min*= scale;
  r_max*= scale;
  r_minsq= r_min*r_min;
  r_maxsq= r_max*r_max;
  deltar= (r_max-r_min)/nbins;
  inv_deltar= 1.0/deltar;
  
  //print file comment lines
  if (fp && me == 0) {
    if (title1) {
      fprintf(fp,"%s\n",title1);
    } else {
      fprintf(fp,"# Spatial-averaged data for fix %s and group %s\n", id,arg[1]);
      fprintf(fp,"# Spherical bins centred at (");
      for (int i= 0; i < 3; ++i) {
        if(origin_type[i] == CONSTANT) fprintf(fp,"%g",origin[i]);
        else fprintf(fp,"%s",arg[i+6]);
        
        if(i != 2) fprintf(fp, ", ");
        else fprintf(fp, ")");
      }
      if (scaleflag == REDUCED) {
        fprintf(fp, " [reduced units]");
      }
      fprintf(fp, "\n");
    }
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# Timestep Number-of-bins\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else {
      fprintf(fp,"# Bin r Ncount");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",arg[12+i]);
      fprintf(fp,"\n");
    }
    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;
  
  //this fix produces a global array
  array_flag= 1;
  size_array_rows= BIG;
  size_array_cols= 1 + 1 + nvalues;
  extarray= 0; //intensive value (don't divide by N when we output this)
  
  //initialization of remaining variables
  irepeat= 0;
  iwindow= window_limit= 0;
  norm= 0;
  maxvar= 0;
  varatom= NULL;
  maxatom= 0;
  bin= NULL;
  binvol= NULL;
  maxbin= 0;
  count_one= count_many= count_sum= count_total= NULL;
  coord= NULL;
  count_list= NULL;
  values_one= values_many= values_sum= values_total= NULL;
  values_list= NULL;
  
  //nvalid = next step on which end_of_step does something
  nvalid= nextvalid();
  //make every compute update itself on step nvalid
  //since we don't yet know which ones we need
  //can correct this when we call end_of_step
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveSpatialSphere::~FixAveSpatialSphere()
{
  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;
  delete [] idregion;
  for (int i= 0; i < 3; ++i) delete [] origin_ids[i];

  if (fp && me == 0) fclose(fp);

  memory->destroy(varatom);
  memory->destroy(bin);
  memory->destroy(binvol);

  memory->destroy(count_one);
  memory->destroy(count_many);
  memory->destroy(count_sum);
  memory->destroy(count_total);
  memory->destroy(coord);
  memory->destroy(count_list);
  memory->destroy(values_one);
  memory->destroy(values_many);
  memory->destroy(values_sum);
  memory->destroy(values_total);
  memory->destroy(values_list);
}

/* ----------------------------------------------------------------------
   tell LAMMPS at which stages in the timestep we want to be called
   we only care about running at the end of a step
------------------------------------------------------------------------- */

int FixAveSpatialSphere::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatialSphere::init()
{
  //set up the region and check it's valid
  if(regionflag) {
    int iregion= domain->find_region(idregion);
    if(iregion == -1) 
      error->all(FLERR,"Region ID for fix ave/spatial/sphere does not exist");
    region = domain->regions[iregion];
  }
  
  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  for (int m = 0; m < 3; ++m) {
    if (origin_type[m] == COMPUTE) {
      int icompute = modify->find_compute(origin_ids[m]);
      if (icompute < 0)
        error->all(FLERR, "Compute ID for fix ave/spatial/sphere does not exist");
      origin_val2idx[m]= icompute;
    } else if (origin_type[m] == VARIABLE) {
      int ivariable = input->variable->find(origin_ids[m]);
      if (ivariable < 0)
        error->all(FLERR, "Variable name for fix ave/spatial/sphere does not exist");
      origin_val2idx[m]= ivariable;
    }
  }

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/spatial/sphere does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/spatial/sphere does not exist");
      value2index[m] = ifix;

      if (nevery % modify->fix[ifix]->peratom_freq)
        error->all(FLERR,
                   "Fix for fix ave/spatial/sphere not computed at compatible time");

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/spatial/sphere does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
  
  
  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed
  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   setup initial bins
   only does averaging if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSpatialSphere::setup(int vflag)
{
  setup_bins();
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSpatialSphere::end_of_step()
{
  //skip this step if we don't need to do anything
  bigint ntimestep= update->ntimestep;
  if(ntimestep != nvalid) return;

  //update region if necessary
  if (regionflag) region->prematch();
  
  //zero out arrays that accumulate over many samples
  if(irepeat == 0) {
    for(int m= 0; m < nbins; m++) {
      count_many[m]= count_sum[m]= 0.0;
      for(int i= 0; i < nvalues; i++) values_many[m][i]= 0.0;
    }
  }

  //if in reduced units, we need to update the bin volumes
  if (scaleflag == REDUCED && domain->box_change) set_bin_volumes();
  
  //zero out arrays for one sample
  for(int m= 0; m < nbins; m++) {
    count_one[m]= 0.0;
    for (int i= 0; i < nvalues; i++) values_one[m][i]= 0.0;
  }
  
  //bin each atom
  double **x= atom->x;
  int *mask= atom->mask;
  int nlocal= atom->nlocal;
  
  if (nlocal > maxatom) {
    maxatom= atom->nmax;
    memory->destroy(bin);
    memory->create(bin,maxatom,"ave/spatial/sphere:bin");
  }
  
  // perform the computation for one sample
  // accumulate results of attributes,computes,fixes,variables to local copy
  // sum within each bin, only include atoms in fix group
  // compute/fix/variable may invoke computes so wrap with clear/add
  modify->clearstep_compute();

  //allocates the bins array
  //each atom gets either the number of the bin into which
  //it fits, or -1
  bin_atoms();
  
  for(int m= 0; m < nvalues; m++) {
    int n= value2index[m];
    int j= argindex[m];
    
    if(which[m] == V || which[m] == F) {
      double **attribute;
      if(which[m] == V) attribute= atom->v;
      else attribute= atom->f;

      for(int i= 0; i < nlocal; i++) {
        if(mask[i] & groupbit && bin[i] > -1) {
          values_one[bin[i]][m] += attribute[i][j];
        }
      }
    } else if(which[m] == DENSITY_NUMBER) {
      for(int i= 0; i < nlocal; i++) {
        if(mask[i] & groupbit && bin[i] > -1) {
          values_one[bin[i]][m] += 1.0;
        }
      }
      // density is additionally normalized by bin volume
      for (int j = 0; j < nbins; j++)
            values_one[j][m] /= binvol[j];
    } else if(which[m] == DENSITY_MASS) {
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;
      for(int i= 0; i < nlocal; i++) {
        if(mask[i] & groupbit && bin[i] > -1) {
          if (rmass) values_one[bin[i]][m] += rmass[i];
          else values_one[bin[i]][m] += mass[type[i]];
        }
      }
      // density is additionally normalized by bin volume
      for (int j = 0; j < nbins; j++)
            values_one[j][m] /= binvol[j];
    } else if(which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }
      double *vector = compute->vector_atom;
      double **array = compute->array_atom;
      int jm1 = j - 1;

      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && bin[i] > -1) {
          if (j == 0) values_one[bin[i]][m] += vector[i];
          else values_one[bin[i]][m] += array[i][jm1];
        }

    } else if(which[m] == FIX) {
      double *vector = modify->fix[n]->vector_atom;
      double **array = modify->fix[n]->array_atom;
      int jm1 = j - 1;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && bin[i] > -1) {
          if (j == 0) values_one[bin[i]][m] += vector[i];
          else values_one[bin[i]][m] += array[i][jm1];
        }
    } else if(which[m] == VARIABLE) {
      if (nlocal > maxvar) {
        maxvar = atom->nmax;
        memory->destroy(varatom);
        memory->create(varatom,maxvar,"ave/spatial:varatom");
      }

      input->variable->compute_atom(n,igroup,varatom,1,0);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit && bin[i] > -1)
          values_one[bin[i]][m] += varatom[i];
    }
  }
  
  // process a single sample
  // if normflag = ALL, accumulate values,count separately to many
  // if normflag = SAMPLE, one = value/count, accumulate one to many
  // exception is SAMPLE density: no normalization by atom count

  if (normflag == ALL) {
    for (int m = 0; m < nbins; m++) {
      count_many[m] += count_one[m];
      for (int j = 0; j < nvalues; j++)
        values_many[m][j] += values_one[m][j];
    }
  } else {
    MPI_Allreduce(count_one,count_many,nbins,MPI_DOUBLE,MPI_SUM,world);
    for (int m = 0; m < nbins; m++) {
      if (count_many[m] > 0.0)
        for (int j = 0; j < nvalues; j++) {
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
  // exception is densities: normalized by repeat, not total count
  
  double repeat = nrepeat;
  double mv2d = force->mv2d;
  if (normflag == ALL) {
    MPI_Allreduce(count_many,count_sum,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nbins*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
    for (int m = 0; m < nbins; m++) {
      if (count_sum[m] > 0.0)
        for (int j = 0; j < nvalues; j++)
          if (which[j] == DENSITY_NUMBER) values_sum[m][j] /= repeat;
          else if (which[j] == DENSITY_MASS) values_sum[m][j] *= mv2d/repeat;
          else values_sum[m][j] /= count_sum[m];
      count_sum[m] /= repeat;
    }
  } else {
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nbins*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
    for (int m = 0; m < nbins; m++) {
      for (int j = 0; j < nvalues; j++)
        values_sum[m][j] /= repeat;
      count_sum[m] /= repeat;
    }
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, comine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++)
        values_total[m][i] = values_sum[m][i];
      count_total[m] = count_sum[m];
    }
    norm = 1;

  } else if (ave == RUNNING) {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++)
        values_total[m][i] += values_sum[m][i];
      count_total[m] += count_sum[m];
    }
    norm++;

  } else if (ave == WINDOW) {
    for (int m = 0; m < nbins; m++) {
      for (int i = 0; i < nvalues; i++) {
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
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d\n",ntimestep,nbins);
    for (int m = 0; m < nbins; m++) {
      fprintf(fp,"  %d %g %g",m+1,coord[m],
              count_total[m]/norm);
      for (int i = 0; i < nvalues; i++)
        fprintf(fp," %g",values_total[m][i]/norm);
      fprintf(fp,"\n");
    }

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      ftruncate(fileno(fp),fileend);
    }
  }
}

/* ----------------------------------------------------------------------
   setup the bins
   called at setup()
------------------------------------------------------------------------- */

void FixAveSpatialSphere::setup_bins()
{
  if(nbins > maxbin) {
    maxbin= nbins;
    memory->grow(count_one,nbins,"ave/spatial/sphere:count_one");
    memory->grow(count_many,nbins,"ave/spatial/sphere:count_many");
    memory->grow(count_sum,nbins,"ave/spatial/sphere:count_sum");
    memory->grow(count_total,nbins,"ave/spatial/sphere:count_total");

    memory->grow(coord,nbins,"ave/spatial/sphere:coord");
    memory->grow(binvol,nbins, "ave/spatial/sphere:binvol");
    memory->grow(values_one,nbins,nvalues,"ave/spatial/sphere:values_one");
    memory->grow(values_many,nbins,nvalues,"ave/spatial/sphere:values_many");
    memory->grow(values_sum,nbins,nvalues,"ave/spatial/sphere:values_sum");
    memory->grow(values_total,nbins,nvalues,"ave/spatial/sphere:values_total");
    
    //the count and value lists are only relevant for windowed averaging
    if(ave == WINDOW) {
      memory->create(count_list,nwindow,nbins,"ave/spatial/sphere:count_list");
      memory->create(values_list,nwindow,nbins,nvalues,
                     "ave/spatial/sphere:values_list");
    }
    
    //reinitialize the regrown count and values
    for(int m= 0; m < nbins; m++) {
      for(int i= 0; i < nvalues; i++) values_total[m][i]= 0.0;
      count_total[m]= 0.0;
    }
  }
  
  //set the bin coordinates
  for(int i= 0; i < nbins; i++) {
    coord[i]= r_min + (i+0.5)*deltar;
    double temp_r= r_min + (i+1)*deltar;
  }

  set_bin_volumes();
}

/* ----------------------------------------------------------------------
   set the bin volumes
   called at setup() and when averaging occurs if box size changes
------------------------------------------------------------------------- */

void FixAveSpatialSphere::set_bin_volumes()
{
  double last_vol= THIRD*MY_4PI*r_minsq*r_min;
  for(int i= 0; i < nbins; i++) {
    double temp_r= r_min + (i+1)*deltar;
    double temp_vol= THIRD*MY_4PI*temp_r*temp_r*temp_r;
    binvol[i]= temp_vol-last_vol;
    last_vol= temp_vol;
  }

  //if in reduced co-ordinates, need to adjust the bin volumes appropriately
  if (scaleflag == REDUCED) {
    double xprd= domain->xprd;
    double yprd= domain->yprd;
    double zprd= domain->zprd;
    for(int i= 0; i < nbins; i++) {
      binvol[i]*= xprd*yprd*zprd;
    }
  }
}

/* ----------------------------------------------------------------------
   assign atoms to bins
------------------------------------------------------------------------- */

void FixAveSpatialSphere::bin_atoms() 
{
  //update the origin coordinates if necessary
  //unscale the existing coordinates
  if(scaleflag == REDUCED) domain->lamda2x(origin, origin);

  for (int i = 0; i < 3; ++i) {
    if (origin_type[i] == COMPUTE) {
      Compute *compute = modify->compute[origin_val2idx[i]];
      if (origin_index[i] == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        origin[i] = compute->scalar;
       } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        origin[i] = compute->vector[origin_index[i]-1];
      }
    } else if (origin_type[i] == VARIABLE) {
      origin[i] = input->variable->compute_equal(origin_val2idx[i]);
    }
  }
  //make sure that the origin point is INSIDE the box
  domain->remap(origin);
  //rescale if necessary
  if(scaleflag == REDUCED) domain->x2lamda(origin, origin);

  double **x= atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double *boxlo,*boxhi,*prd,*prd_half;
  //deal with the periodic boundary conditions
  if(domain->periodicity[0] || domain->periodicity[1] || domain->periodicity[2]) {
    if(scaleflag == REDUCED) {
      boxlo= domain->boxlo_lamda;
      boxhi= domain->boxhi_lamda;
      prd= domain->prd_lamda;
      prd_half= domain->prd_half_lamda;
    } else {
      boxlo= domain->boxlo;
      boxhi= domain->boxhi;
      prd= domain->prd;
      prd_half= domain->prd_half;
    }
  }
 
  double rsq, tmp;
  
  if(!regionflag) {
    if(scaleflag == REDUCED) {
      //convert the atomic coordinates into reduced coordinates
      domain->x2lamda(nlocal);
    }
    for(int i= 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {      
        rsq= 0.0;
        for(int dim= 0; dim < 3; dim++) {
          tmp= x[i][dim] - origin[dim];
          if(domain->periodicity[dim]) {
            if(tmp < -prd_half[dim]) tmp+= prd[dim];
            else if(tmp >= prd_half[dim]) tmp-= prd[dim];
          }
          rsq+= tmp*tmp;
        }
        if(rsq >= r_minsq && rsq <= r_maxsq) {
          double r= sqrt(rsq);
          int ibin= static_cast<int>((r-r_min)*inv_deltar);
          ibin = MAX(ibin,0);
          ibin = MIN(ibin,nbins-1);
          bin[i]= ibin;
          count_one[ibin]+= 1.0;
        } else {
          bin[i]= -1;
        }
      }
    }
    if(scaleflag == REDUCED) {
      //convert the atomic coordinates back into "real" coordinates
      domain->lamda2x(nlocal);
    }
  } else {
    for(int i= 0; i < nlocal; i++) {
      if(mask[i] & groupbit && region->match(x[i][0], x[i][1], x[i][2])) {
        rsq= 0.0;
        for(int dim= 0; dim < 3; dim++) {
          tmp= x[i][dim] - origin[dim];
          if(domain->periodicity[dim]) {
            if(tmp < -0.5*prd[dim]) tmp+= prd[dim];
            else if(tmp >= 0.5*prd[dim]) tmp-= prd[dim];
          }
          rsq+= tmp*tmp;
        }
        if(rsq >= r_minsq && rsq <= r_maxsq) {
          double r= sqrt(rsq);
          int ibin= static_cast<int>((r-r_min)*inv_deltar);
          ibin = MAX(ibin,0);
          ibin = MIN(ibin,nbins-1);
          bin[i]= ibin;
          count_one[ibin]+= 1.0;
        } else {
          bin[i]= -1;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return I,J array value
   if I exceeds current bins, return 0.0 instead of generating an error
   column 1 = bin coords, next column = count, remaining columns = Nvalues
------------------------------------------------------------------------- */

double FixAveSpatialSphere::compute_array(int i, int j)
{
  if (values_total == NULL) return 0.0;
  if (i >= nbins || j > nvalues) return 0.0;
  if (!norm) return 0.0;
  if (j == 0) return coord[i];
  if (j == 1) return count_total[i]/norm;
  return values_total[i][j]/norm;
}


/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveSpatialSphere::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ----------------------------------------------------------------------
   memory usage of varatom and bins
------------------------------------------------------------------------- */

double FixAveSpatialSphere::memory_usage()
{
  double bytes = maxvar * sizeof(double);         // varatom
  bytes += maxatom * sizeof(int);                 // bin
  bytes += 4*nbins * sizeof(double);              // count one,many,sum,total
  bytes += nbins * sizeof(double);           // coord
  bytes += nvalues*nbins * sizeof(double);        // values one,many,sum,total
  bytes += nwindow*nbins * sizeof(double);          // count_list
  bytes += nwindow*nbins*nvalues * sizeof(double);  // values_list
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatialSphere::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) 
    error->all(FLERR,"Fix ave/spatial/sphere missed timestep");
}
