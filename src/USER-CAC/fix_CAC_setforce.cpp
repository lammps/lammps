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

#include <cstring>
#include <cstdlib>
#include "fix_CAC_setforce.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixCAC_Set_Force::FixCAC_Set_Force(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), idregion(NULL), sforce(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix cac/setforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix cac/setforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix cac/setforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix cac/setforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,3,"cac_setforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixCAC_Set_Force::~FixCAC_Set_Force()
{
  if (copymode) return;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixCAC_Set_Force::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_Force::init()
{
	// check variables
  if (!atom->CAC_flag) error->all(FLERR,"fix cac/setforce requires a CAC atom style");
	if (xstr) {
		xvar = input->variable->find(xstr);
		if (xvar < 0)
			error->all(FLERR, "Variable name for fix cac/setforce does not exist");
		if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
		else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
		else error->all(FLERR, "Variable for fix cac/setforce is invalid style");
	}
	if (ystr) {
		yvar = input->variable->find(ystr);
		if (yvar < 0)
			error->all(FLERR, "Variable name for fix cac/setforce does not exist");
		if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
		else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
		else error->all(FLERR, "Variable for fix cac/setforce is invalid style");
	}
	if (zstr) {
		zvar = input->variable->find(zstr);
		if (zvar < 0)
			error->all(FLERR, "Variable name for fix cac/setforce does not exist");
		if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
		else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
		else error->all(FLERR, "Variable for fix cac/setforce is invalid style");
	}

	// set index and check validity of region

	if (iregion >= 0) {
		iregion = domain->find_region(idregion);
		if (iregion == -1)
			error->all(FLERR, "Region ID for fix cac/setforce does not exist");
	}

	if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
		varflag = ATOM;
	else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
		varflag = EQUAL;
	else varflag = CONSTANT;


	// cannot use non-zero forces for a minimization since no energy is integrated
	// use fix addforce instead

	int flag = 0;
	if (update->whichflag == 2) {
		if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
		if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
		if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
		if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
		if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
		if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
	}
	if (flag)
		error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_Force::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "Cannot use respa with cac/setforce");
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_Force::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_Force::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double ****nodal_forces = atom->nodal_forces;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;	
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"cac_setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;
  double ****nodal_positions = atom->nodal_positions;
  //double ****initial_nodal_positions = atom->initial_nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
 

  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;

  if (varflag == CONSTANT) {
	  for (int i = 0; i < nlocal; i++) {
		  if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
			  for (int j = 0; j < nodes_per_element; j++) {
				  for (int l = 0; l < poly_count[i]; l++) {
					  if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
					  foriginal[0] += nodal_forces[i][j][l][0];
					  foriginal[1] += nodal_forces[i][j][l][1];
					  foriginal[2] += nodal_forces[i][j][l][2];
					  if (xstyle) nodal_forces[i][j][l][0] = xvalue;
					  if (ystyle) nodal_forces[i][j][l][1] = yvalue;
					  if (zstyle) nodal_forces[i][j][l][2] = zvalue;
				  }
			  }
		  }
	  }
  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
			  for (int j = 0; j < nodes_per_element; j++) {
				 for (int l = 0; l < poly_count[i]; l++) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += nodal_forces[i][j][l][0];
        foriginal[1] += nodal_forces[i][j][l][1];
        foriginal[2] += nodal_forces[i][j][l][2];
        if (xstyle == ATOM) nodal_forces[i][j][l][0] = sforce[i][0];
        else if (xstyle) nodal_forces[i][j][l][0] = xvalue;
        if (ystyle == ATOM) nodal_forces[i][j][l][1] = sforce[i][1];
        else if (ystyle) nodal_forces[i][j][l][1] = yvalue;
        if (zstyle == ATOM) nodal_forces[i][j][l][2] = sforce[i][2];
        else if (zstyle) nodal_forces[i][j][l][2] = zvalue;
         }
        }
      }
    }
  }
  
}


/* ---------------------------------------------------------------------- */

void FixCAC_Set_Force::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixCAC_Set_Force::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixCAC_Set_Force::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*3 * sizeof(double);
  return bytes;
}
