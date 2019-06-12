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
#include <stdlib.h>
#include "fix_CAC_setdof.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include <cmath>
#define PI 3.14159265

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixCAC_Set_DOF::FixCAC_Set_DOF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix setforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = nlevels_respa = 0;
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
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix setforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix setforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix setforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,3,"setforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixCAC_Set_DOF::~FixCAC_Set_DOF()
{
  if (copymode) return;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixCAC_Set_DOF::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::init()
{
	// check variables

	if (xstr) {
		xvar = input->variable->find(xstr);
		if (xvar < 0)
			error->all(FLERR, "Variable name for fix setforce does not exist");
		if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
		else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
		else error->all(FLERR, "Variable for fix setforce is invalid style");
	}
	if (ystr) {
		yvar = input->variable->find(ystr);
		if (yvar < 0)
			error->all(FLERR, "Variable name for fix setforce does not exist");
		if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
		else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
		else error->all(FLERR, "Variable for fix setforce is invalid style");
	}
	if (zstr) {
		zvar = input->variable->find(zstr);
		if (zvar < 0)
			error->all(FLERR, "Variable name for fix setforce does not exist");
		if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
		else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
		else error->all(FLERR, "Variable for fix setforce is invalid style");
	}

	// set index and check validity of region

	if (iregion >= 0) {
		iregion = domain->find_region(idregion);
		if (iregion == -1)
			error->all(FLERR, "Region ID for fix setforce does not exist");
	}

	if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
		varflag = ATOM;
	else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
		varflag = EQUAL;
	else varflag = CONSTANT;

	if (strstr(update->integrate_style, "respa")) {
		nlevels_respa = ((Respa *)update->integrate)->nlevels;
		if (respa_level >= 0) ilevel_respa = MIN(respa_level, nlevels_respa - 1);
		else ilevel_respa = nlevels_respa - 1;
	}

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

	double ****nodal_positions = atom->nodal_positions;
	double ****nodal_velocities = atom->nodal_velocities;
	int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;	
	double ****initial_nodal_positions = atom->initial_nodal_positions;
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	int **node_types = atom->node_types;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	double *boxlo = domain->boxlo;
	double *boxhi = domain->boxhi;
	double origin[3];
	double distance;
	double delr[3];
	double center_radius_sq;
	double wave_vector;
	double Amplitude = 0.003;
	double cell = 10 * 5.430950;
	double omega;

	omega = 43 * wave_vector;
	double unit_vector[2];
	center_radius_sq = 56.25 * 4*cell * cell;
	wave_vector = 2 * PI / (6 * cell);
	origin[0] = (boxhi[0] + boxlo[0]) / 2;
	origin[1] = (boxhi[1] + boxlo[1]) / 2;
	origin[2] = (boxhi[2] + boxlo[2]) / 2;
	
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
      nodes_per_element = nodes_count_list[element_type[i]];
			for (int j = 0; j < nodes_per_element; j++) {
				for (int l = 0; l < poly_count[i]; l++) {

					delr[0] = (nodal_positions[i][j][l][0] - origin[0]);
					delr[1] = (nodal_positions[i][j][l][1] - origin[1]);
					delr[2] = (nodal_positions[i][j][l][2] - origin[2]);
					distance = delr[0] * delr[0] + delr[1] * delr[1];
					if (delr[0] == 0 && delr[1] == 0) {
						unit_vector[0] = sqrt(2) / 2;
						unit_vector[1] = sqrt(2) / 2;
					}
					else {
						unit_vector[0] = delr[0] / sqrt(distance);
						unit_vector[1] = delr[1] / sqrt(distance);
					}
					if (distance <= center_radius_sq) {
						//nodal_positions[i][j][l][0] = initial_nodal_positions[i][j][l][0] + Amplitude * cos(sqrt(distance)*wave_vector )*unit_vector[0];
						//nodal_positions[i][j][l][1] = initial_nodal_positions[i][j][l][1] + Amplitude * cos(sqrt(distance)*wave_vector )*unit_vector[1];
						//nodal_velocities[i][j][l][0] = Amplitude * omega * sin(sqrt(distance)*wave_vector )*unit_vector[0];
						//nodal_velocities[i][j][l][1] = Amplitude * omega * sin(sqrt(distance)*wave_vector )*unit_vector[1];
						//n1z += Amplitude*cos(sqrt(distance)*wave_vector)*delr[2]/ distance;
					//if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

						nodal_velocities[i][j][l][2] += 0;
					}
				}
			}
		}
	}
	
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::post_force(int vflag)
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
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;
  double ****nodal_positions = atom->nodal_positions;
  double ****initial_nodal_positions = atom->initial_nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
 

  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double origin[3];
  double distance;
  double distancex;
  double delr[3];
  double delrx[3];
  double center_radius_sq;
  double wave_vector;
  double Amplitude = 0.01;
  double cell = 10 * 5.430950;
  double unit_vector[2];
  double omega;
  double current_time = update->dt*update->ntimestep;
  center_radius_sq = 56.25 * 4*cell * cell;
  wave_vector = 2 * PI / (6 * cell);
  omega = 43 * wave_vector;
  origin[0] = (boxhi[0] + boxlo[0]) / 2;
  origin[1] = (boxhi[1] + boxlo[1]) / 2;
  origin[2] = (boxhi[2] + boxlo[2]) / 2;
  /*
  for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
		  x[i][0] = 0;
		  x[i][1] = 0;
		  x[i][2] = 0;
		  v[i][0] = 0;
		  v[i][1] = 0;
		  v[i][2] = 0;
		  for (int j = 0; j < nodes_per_element; j++) {
			  for (int l = 0; l < poly_count[i]; l++) {

				  delr[0] = (nodal_positions[i][j][l][0] - origin[0]);
				  delr[1] = (nodal_positions[i][j][l][1] - origin[1]);
				  delr[2] = (nodal_positions[i][j][l][2] - origin[2]);
				  delrx[0] = (x[i][0] - origin[0]);
				  delrx[1] = (x[i][1] - origin[1]);
				  delrx[2] = (x[i][2] - origin[2]);

				  distance = delr[0] * delr[0] + delr[1] * delr[1];
				  distancex = delrx[0] * delrx[0] + delrx[1] * delrx[1];
				  if (delr[0] == 0 && delr[1] == 0) {
					  unit_vector[0] = -sqrt(2) / 2;
					  unit_vector[1] = sqrt(2) / 2;
				  }
				  else {
					  unit_vector[0] = -delr[1] / sqrt(distance);
					  unit_vector[1] = delr[0] / sqrt(distance);
				  }
				 

					  //Amplitude = 0.01*exp(-2*distance / center_radius_sq);
					  nodal_positions[i][j][l][0] = initial_nodal_positions[i][j][l][0]+Amplitude * cos(sqrt(distance)*wave_vector- omega*current_time)*unit_vector[0];
					  nodal_positions[i][j][l][1] = initial_nodal_positions[i][j][l][1]+Amplitude * cos(sqrt(distance)*wave_vector - omega*current_time)*unit_vector[1];
					  nodal_velocities[i][j][l][0] = Amplitude *  omega * sin(sqrt(distance)*wave_vector - omega*current_time)*unit_vector[0];
					  nodal_velocities[i][j][l][1] = Amplitude *  omega * sin(sqrt(distance)*wave_vector - omega*current_time)*unit_vector[1];
					  //n1z += Amplitude*cos(sqrt(distance)*wave_vector)*delr[2]/ distance;
					  //if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

					  //nodal_velocities[i][j][l][2] += 0;
					  x[i][0] += nodal_positions[i][j][l][0];
					  x[i][1] += nodal_positions[i][j][l][1];
					  x[i][2] += nodal_positions[i][j][l][2];
					  v[i][0] += nodal_velocities[i][j][l][0];
					  v[i][1] += nodal_velocities[i][j][l][1];
					  v[i][2] += nodal_velocities[i][j][l][2];
			  }
		  }
		  x[i][0] = x[i][0] / nodes_per_element / poly_count[i];
		  x[i][1] = x[i][1] / nodes_per_element / poly_count[i];
		  x[i][2] = x[i][2] / nodes_per_element / poly_count[i];
		  v[i][0] = v[i][0] / nodes_per_element / poly_count[i];
		  v[i][1] = v[i][1] / nodes_per_element / poly_count[i];
		  v[i][2] = v[i][2] / nodes_per_element / poly_count[i];
	  }
  }
  */
  
  
  if (varflag == CONSTANT) {
	  for (int i = 0; i < nlocal; i++) {
		  nodes_per_element = nodes_count_list[element_type[i]];
		  if (mask[i] & groupbit) {
			  for (int j = 0; j < nodes_per_element; j++) {
				  for (int l = 0; l < poly_count[i]; l++) {
					  if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
					  foriginal[0] += f[i][0];
					  foriginal[1] += f[i][1];
					  foriginal[2] += f[i][2];
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

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (xstyle == ATOM) f[i][0] = sforce[i][0];
        else if (xstyle) f[i][0] = xvalue;
        if (ystyle == ATOM) f[i][1] = sforce[i][1];
        else if (ystyle) f[i][1] = yvalue;
        if (zstyle == ATOM) f[i][2] = sforce[i][2];
        else if (zstyle) f[i][2] = zvalue;
      }
  }
  
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::post_force_respa(int vflag, int ilevel, int iloop)
{
  // set force to desired value on requested level, 0.0 on other levels

  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = NULL;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }

    double **x = atom->x;
    double **f = atom->f;
	
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if (xstyle) f[i][0] = 0.0;
        if (ystyle) f[i][1] = 0.0;
        if (zstyle) f[i][2] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixCAC_Set_DOF::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixCAC_Set_DOF::compute_vector(int n)
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

double FixCAC_Set_DOF::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*3 * sizeof(double);
  return bytes;
}
