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
#include "fix_pair_tracking.h"
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

FixPairTracking::FixPairTracking(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0),
  array(NULL), vector(NULL), pack_choice(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix pair/tracker command");
  store_flag = 0;
  local_flag = 1;
  nvalues = narg - 4; 
  
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix pair/tracker command");  
    
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pack_choice = new FnPtrPack[nvalues];
  
  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg - 4;
    
    if (strcmp(arg[iarg],"id1") == 0) {
      pack_choice[i] = &FixPairTracking::pack_id1;
    } else if (strcmp(arg[iarg],"id2") == 0) {
      pack_choice[i] = &FixPairTracking::pack_id2;
           
    } else if (strcmp(arg[iarg],"time/created") == 0) {
      pack_choice[i] = &FixPairTracking::pack_time_created;
    } else if (strcmp(arg[iarg],"time/broken") == 0) {
      pack_choice[i] = &FixPairTracking::pack_time_broken;
    } else if (strcmp(arg[iarg],"time/total") == 0) {
      pack_choice[i] = &FixPairTracking::pack_time_total;
      
    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &FixPairTracking::pack_x;    
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &FixPairTracking::pack_y;    
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &FixPairTracking::pack_z;    
      
    } else if (strcmp(arg[iarg],"xstore") == 0) {
      pack_choice[i] = &FixPairTracking::pack_xstore;    
      store_flag = 1;
    } else if (strcmp(arg[iarg],"ystore") == 0) {
      pack_choice[i] = &FixPairTracking::pack_ystore;    
      store_flag = 1;
    } else if (strcmp(arg[iarg],"zstore") == 0) {
      pack_choice[i] = &FixPairTracking::pack_zstore;      
      store_flag = 1;      
   
    } else if (strcmp(arg[iarg],"rmin") == 0) {
      pack_choice[i] = &FixPairTracking::pack_rmin;    
    } else if (strcmp(arg[iarg],"rave") == 0) {
      pack_choice[i] = &FixPairTracking::pack_rave;     
   
    } else error->all(FLERR, "Invalid keyword in fix pair/tracker command");
  }

  nmax = 0;
  ncount = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

FixPairTracking::~FixPairTracking()
{
  if (modify->nfix & store_flag == 1) {
    modify->delete_fix(id_fix);
    delete [] id_fix;
  }
  
  delete [] pack_choice;

  memory->destroy(vector);  
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixPairTracking::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::post_constructor()
{
  //If use stored x,y,z values, use store property (can transfer to ghost atoms) to store positions
  
  if(store_flag == 1){
      
    lx = domain->xprd;
    ly = domain->yprd;
    lz = domain->zprd;        

    int nn = strlen(id) + strlen("_FIX_PROP_ATOM") + 1;
    id_property_fix = new char[nn];
    strcpy(id_property_fix,id);
    strcat(id_property_fix,"_FIX_PROP_ATOM");
    
    int ifix = modify->find_fix(id_property_fix);
    if (ifix < 0) {
    
      int n_x = strlen(id) + 4;
      
      char * lab1 = new char[n_x];
      strcpy(lab1, "d_x");
      strcat(lab1, id);
      char * lab2 = new char[n_x];
      strcpy(lab2, "d_y");
      strcat(lab2, id);
      char * lab3 = new char[n_x];
      strcpy(lab3, "d_z");
      strcat(lab3, id);
        
      char **newarg = new char*[8];
      newarg[0] = id_property_fix;
      newarg[1] = group->names[igroup];
      newarg[2] = (char *) "property/atom"; 
      newarg[3] = (char *) lab1;         
      newarg[4] = (char *) lab2;         
      newarg[5] = (char *) lab3; 
      newarg[6] = (char *) "ghost";
      newarg[7] = (char *) "yes";

      modify->add_fix(8,newarg); 
      //Needs ghost atoms to calculate CoM

      int type_flag;
      int col_flag;

      strcpy(lab1, "x");
      strcat(lab1, id);
      strcpy(lab2, "y");
      strcat(lab2, id);
      strcpy(lab3, "z");
      strcat(lab3, id);

      index_x = atom->find_custom(lab1, type_flag, col_flag);
      index_y = atom->find_custom(lab2, type_flag, col_flag);
      index_z = atom->find_custom(lab3, type_flag, col_flag);
      delete [] newarg;    
      delete [] lab1;
      delete [] lab2;
      delete [] lab3;
    } 
    
    ifix = modify->find_fix(id_fix);
    if (ifix < 0) error->all(FLERR,"Could not find fix ID for fix broken/bond");
    if (modify->fix[ifix]->restart_reset) {
        modify->fix[ifix]->restart_reset = 0;
    } else {

      double *xi = atom->dvector[index_x];
      double *yi = atom->dvector[index_y];
      double *zi = atom->dvector[index_z];
      
      double **xs = atom->x;
      int nlocal = atom->nlocal;
      
      for (int i = 0; i < nlocal; i++) {
        xi[i] = xs[i][0];
        yi[i] = xs[i][1];
        zi[i] = xs[i][2];   
      }
    }    
  }
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::init()
{
  // Set size of array/vector
  ncount = 0;
  
  if (ncount > nmax) {
      reallocate(ncount);}
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::lost_contact(int i, int j, double n, double rs, double rm)
{    
  if (ncount == nmax) reallocate(ncount);
  
  index_i = i;
  index_j = j;
  rmin = rm;
  rave = ra;
  ntimestep = n;
  
  // fill vector or array with local values
  if (nvalues == 1) {
    (this->*pack_choice[0])(0);
  } else {
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n); 
  }  
  
  ncount += 1;  
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::post_force(int /*vflag*/) 
{    
  if (update->ntimestep % nevery == 0) {
    size_local_rows = ncount;
    ncount = 0;  
  }
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::reallocate(int n)
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

double FixPairTracking::memory_usage()
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

void FixPairTracking::pack_time_created(int n) 
{
  if (nvalues == 1)
    vector[ncount] = ntimestep;
  else
    array[ncount][n] = ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_time_broken(int n) 
{
  if (nvalues == 1)
    vector[ncount] = update->ntimestep;
  else
    array[ncount][n] = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_time_total(int n) 
{
  if (nvalues == 1)
    vector[ncount] = update->ntimestep-ntimestep;
  else
    array[ncount][n] = update->ntimestep-ntimestep;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_id1(int n) 
{
  tagint *tag = atom->tag;
  
  if (nvalues == 1)
    vector[ncount] = tag[index_i];
  else
    array[ncount][n] = tag[index_i];
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_id2(int n)
{
  tagint *tag = atom->tag;

  if (nvalues == 1)
    vector[ncount] = tag[index_j];
  else
    array[ncount][n] = tag[index_j];
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_x(int n)
{
  double lx_new = domain->xprd;
  double **x = atom->x; 
    
  if (nvalues == 1)
    vector[ncount] = (x[index_i][0] + x[index_j][0])/2;
  else
    array[ncount][n] = (x[index_i][0] + x[index_j][0])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_y(int n)
{
  double lx_new = domain->yprd;
  double **x = atom->x; 
  
  if (nvalues == 1)
    vector[ncount] = (x[index_i][1] + x[index_j][1])/2;
  else
    array[ncount][n] = (x[index_i][1] + x[index_j][1])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_z(int n)
{
  double lx_new = domain->zprd;
  double **x = atom->x; 
    
  if (nvalues == 1)
    vector[ncount] = (x[index_i][2] + x[index_j][2])/2;
  else
    array[ncount][n] = (x[index_i][2] + x[index_j][2])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_xstore(int n)
{
  double *x = atom->dvector[index_x];
        
  if (nvalues == 1)
    vector[ncount] = (x[index_i] + x[index_j])/2;
  else
    array[ncount][n] = (x[index_i] + x[index_j])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_ystore(int n)
{
  double *y = atom->dvector[index_y];
        
  if (nvalues == 1)
    vector[ncount] = (y[index_i] + y[index_j])/2;
  else
    array[ncount][n] = (y[index_i] + y[index_j])/2;
}


/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_zstore(int n)
{
  double *z = atom->dvector[index_z];
        
  if (nvalues == 1)
    vector[ncount] = (z[index_i] + z[index_j])/2;
  else
    array[ncount][n] = (z[index_i] + z[index_j])/2;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_rmin(int n)
{
  if (nvalues == 1)
    vector[ncount] = rmin;
  else
    array[ncount][n] = rmin;
}

/* ---------------------------------------------------------------------- */

void FixPairTracking::pack_rave(int n)
{
  if (nvalues == 1)
    vector[ncount] = rsum/(update->ntimestep-ntimestep);
  else
    array[ncount][n] = rsum/(update->ntimestep-ntimestep);
}
