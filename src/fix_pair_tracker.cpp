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

#include <string.h>
#include "fix_pair_tracker.h"
#include "pair_tracker.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "group.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 100

/* ---------------------------------------------------------------------- */

FixPairTracker::FixPairTracker(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0),
  array(NULL), vector(NULL), pack_choice(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix pair/tracker command");
  local_flag = 1;
  nvalues = narg - 4; 
  tmin = -1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix pair/tracker command");  
    
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pack_choice = new FnPtrPack[nvalues];
  
  int iarg = 4;
  int i = 0;
  while (iarg < narg) {
    
    if (strcmp(arg[iarg],"id1") == 0) {
      pack_choice[i] = &FixPairTracker::pack_id1;
    } else if (strcmp(arg[iarg],"id2") == 0) {
      pack_choice[i] = &FixPairTracker::pack_id2;
           
    } else if (strcmp(arg[iarg],"time/created") == 0) {
      pack_choice[i] = &FixPairTracker::pack_time_created;
    } else if (strcmp(arg[iarg],"time/broken") == 0) {
      pack_choice[i] = &FixPairTracker::pack_time_broken;
    } else if (strcmp(arg[iarg],"time/total") == 0) {
      pack_choice[i] = &FixPairTracker::pack_time_total;
      
    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &FixPairTracker::pack_x;    
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &FixPairTracker::pack_y;    
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &FixPairTracker::pack_z;    
      
    } else if (strcmp(arg[iarg],"rmin") == 0) {
      pack_choice[i] = &FixPairTracker::pack_rmin;    
    } else if (strcmp(arg[iarg],"rave") == 0) {
      pack_choice[i] = &FixPairTracker::pack_rave;     

    } else if (strcmp(arg[iarg],"tmin") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Invalid keyword in fix pair/tracker command");
      tmin = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      i -= 1;
      nvalues -= 2;
      iarg ++;
   
    } else error->all(FLERR, "Invalid keyword in fix pair/tracker command");
    
    iarg ++;
    i ++;
  }

  nmax = 0;
  ncount = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

FixPairTracker::~FixPairTracker()
{  
  delete [] pack_choice;

  memory->destroy(vector);  
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixPairTracker::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::init()
{
  // Set size of array/vector
  ncount = 0;
  
  if (ncount > nmax) {
      reallocate(ncount);}
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::lost_contact(int i, int j, double n, double rs, double rm)
{    
  if (update->ntimestep-n > tmin) {
    if (ncount == nmax) reallocate(ncount);
    
    index_i = i;
    index_j = j;

    rmin = rm;
    rsum = rs;
    ntimestep = n;
    
    // fill vector or array with local values
    if (nvalues == 1) {
      (this->*pack_choice[0])(0);
    } else {
      for (int k = 0; k < nvalues; k++) {
        (this->*pack_choice[k])(k); 
      }
    }  

    ncount += 1;  
  }
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::post_force(int /*vflag*/) 
{    
  if (update->ntimestep % nevery == 0) {
    size_local_rows = ncount;
    ncount = 0;  
  }
}


/* ---------------------------------------------------------------------- */

void FixPairTracker::reallocate(int n)
{
  // grow vector or array 
  while (nmax <= n) nmax += DELTA;
  
  if (nvalues == 1) {
    memory->grow(vector,nmax,"fix_broken_bonds:vector");
    vector_local = vector;
  } else {
    memory->grow(array,nmax,nvalues,"fix_broken_bonds:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double FixPairTracker::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  bytes += nmax*2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword fix property/local can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_created(int n) 
{
  if (nvalues == 1)
    vector[ncount] = ntimestep;
  else
    array[ncount][n] = ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_broken(int n) 
{
  if (nvalues == 1)
    vector[ncount] = update->ntimestep;
  else
    array[ncount][n] = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_time_total(int n) 
{
  if (nvalues == 1)
    vector[ncount] = update->ntimestep-ntimestep;
  else
    array[ncount][n] = update->ntimestep-ntimestep;
}


/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_id1(int n) 
{
  tagint *tag = atom->tag;
  
  if (nvalues == 1)
    vector[ncount] = tag[index_i];
  else
    array[ncount][n] = tag[index_i];
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_id2(int n)
{
  tagint *tag = atom->tag;

  if (nvalues == 1)
    vector[ncount] = tag[index_j];
  else
    array[ncount][n] = tag[index_j];
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_x(int n)
{
  double lx_new = domain->xprd;
  double **x = atom->x; 
    
  if (nvalues == 1)
    vector[ncount] = (x[index_i][0] + x[index_j][0])/2;
  else
    array[ncount][n] = (x[index_i][0] + x[index_j][0])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_y(int n)
{
  double lx_new = domain->yprd;
  double **x = atom->x; 
  
  if (nvalues == 1)
    vector[ncount] = (x[index_i][1] + x[index_j][1])/2;
  else
    array[ncount][n] = (x[index_i][1] + x[index_j][1])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_z(int n)
{
  double lx_new = domain->zprd;
  double **x = atom->x; 
    
  if (nvalues == 1)
    vector[ncount] = (x[index_i][2] + x[index_j][2])/2;
  else
    array[ncount][n] = (x[index_i][2] + x[index_j][2])/2;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_rmin(int n)
{
  if (nvalues == 1)
    vector[ncount] = rmin;
  else
    array[ncount][n] = rmin;
}

/* ---------------------------------------------------------------------- */

void FixPairTracker::pack_rave(int n)
{
  if (nvalues == 1)
    vector[ncount] = rsum/(update->ntimestep-ntimestep);
  else
    array[ncount][n] = rsum/(update->ntimestep-ntimestep);
}
