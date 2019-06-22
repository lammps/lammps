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
#include <cmath>
#include "fix_cac_temp_rescale.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempRescale_CAC::FixTempRescale_CAC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  tstr(NULL), id_temp(NULL), tflag(0)
{
  if (narg < 8) error->all(FLERR,"Illegal fix cac/temp/rescale command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix cac/temp/rescale command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;
  dynamic_group_allow = 1;

  tstr = NULL;
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[4][2]);
    tstyle = EQUAL;
  } else {
    t_start = force->numeric(FLERR,arg[4]);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = force->numeric(FLERR,arg[5]);
  t_window = force->numeric(FLERR,arg[6]);
  fraction = force->numeric(FLERR,arg[7]);

  // create a new compute temp
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 16;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_cac_nodal_temp");

  char **newarg = new char*[6];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "cac/nodal/temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempRescale_CAC::~FixTempRescale_CAC()
{
  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempRescale_CAC::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale_CAC::init()
{
  // check variable
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix cac/temp/rescale does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR,"Variable for fix cac/temp/rescale is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix cac/temp/rescale does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale_CAC::end_of_step()
{
  double t_current = temperature->compute_scalar();

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // protect against division by zero

  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix cac/temp/rescale cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix cac/temp/rescale variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    double ****nodal_velocities = atom->nodal_velocities;
	int nodes_per_element;
	int *nodes_count_list = atom->nodes_per_element_list;	
	int *element_type = atom->element_type;
	int *poly_count = atom->poly_count;
	int **node_types = atom->node_types;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    energy += (t_current-t_target) * efactor;

    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
			for (int j = 0; j < nodes_per_element; j++) {
				for (int l = 0; l < poly_count[i]; l++) {
					//temperature->remove_bias(i,v[i]);
					nodal_velocities[i][j][l][0] *= factor;
					nodal_velocities[i][j][l][1] *= factor;
					nodal_velocities[i][j][l][2] *= factor;
					// temperature->restore_bias(i,v[i]);
				}
			}
        }
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];  
			for (int j = 0; j < nodes_per_element; j++) {
				for (int l = 0; l < poly_count[i]; l++) {
					//temperature->remove_bias(i,v[i]);
					nodal_velocities[i][j][l][0] *= factor;
					nodal_velocities[i][j][l][1] *= factor;
					nodal_velocities[i][j][l][2] *= factor;
					// temperature->restore_bias(i,v[i]);
				}
			}
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixTempRescale_CAC::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"cac/nodal_temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale_CAC::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempRescale_CAC::compute_scalar()
{
  return energy;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixTempRescale_CAC::extract(const char *str, int &dim)
{
  if (strcmp(str,"t_target") == 0) {
    dim = 0;
    return &t_target;
  }
  return NULL;
}
