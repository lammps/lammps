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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_press_berendsen.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_deform.h"
#include "compute.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

enum{NOBIAS,BIAS};
enum{XYZ,XY,YZ,XZ,ANISO};

/* ---------------------------------------------------------------------- */

FixPressBerendsen::FixPressBerendsen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix press/berendsen command");

  box_change = 1;

  // Berendsen barostat should be applied every step

  nevery = 1;

  if (strcmp(arg[3],"xyz") == 0) {
    if (narg < 7) error->all("Illegal fix press/berendsen command");

    press_couple = XYZ;
    p_start[0] = p_start[1] = p_start[2] = atof(arg[4]);
    p_stop[0] = p_stop[1] = p_stop[2] = atof(arg[5]);
    p_period[0] = p_period[1] = p_period[2] = atof(arg[6]);
    p_flag[0] = p_flag[1] = p_flag[2] = 1;
    if (domain->dimension == 2) {
      p_start[2] = p_stop[2] = p_period[2] = 0.0;
      p_flag[2] = 0;
    }

  } else {
    if (strcmp(arg[3],"xy") == 0) press_couple = XY;
    else if (strcmp(arg[3],"yz") == 0) press_couple = YZ;
    else if (strcmp(arg[3],"xz") == 0) press_couple = XZ;
    else if (strcmp(arg[3],"aniso") == 0) press_couple = ANISO;
    else error->all("Illegal fix press/berendsen command");

    if (narg < 8) error->all("Illegal fix press/berendsen command");

    if (domain->dimension == 2 && 
	(press_couple == XY || press_couple == YZ || press_couple == XZ))
      error->all("Invalid fix press/berendsen command for a 2d simulation");

    if (strcmp(arg[4],"NULL") == 0) {
      p_start[0] = p_stop[0] = p_period[0] = 0.0;
      p_flag[0] = 0;
    } else {
      p_start[0] = atof(arg[4]);
      p_stop[0] = atof(arg[5]);
      p_flag[0] = 1;
    }
    if (strcmp(arg[6],"NULL") == 0) {
      p_start[1] = p_stop[1] = p_period[1] = 0.0;
      p_flag[1] = 0;
    } else {
      p_start[1] = atof(arg[6]);
      p_stop[1] = atof(arg[7]);
      p_flag[1] = 1;
    }
    if (strcmp(arg[8],"NULL") == 0) {
      p_start[2] = p_stop[2] = p_period[2] = 0.0;
      p_flag[2] = 0;
    } else {
      if (domain->dimension == 2)
	error->all("Invalid fix press/berendsen command for a 2d simulation");
      p_start[2] = atof(arg[8]);
      p_stop[2] = atof(arg[9]);
      p_flag[2] = 1;
    }

    double period = atof(arg[10]);
    if (p_flag[0]) p_period[0] = period;
    if (p_flag[1]) p_period[1] = period;
    if (p_flag[2]) p_period[2] = period;
  }

  // process extra keywords

  bulkmodulus = 10.0;
  allremap = 1;

  int iarg;
  if (press_couple == XYZ) iarg = 7;
  else iarg = 11;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"modulus") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix press/berendsen command");
      bulkmodulus = atof(arg[iarg+1]);
      if (bulkmodulus <= 0.0) 
	error->all("Illegal fix press/berendsen command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix press/berendsen command");
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1],"partial") == 0) allremap = 0;
      else error->all("Illegal fix press/berendsen command");
      iarg += 2;
    } else error->all("Illegal fix press/berendsen command");
  }

  // error checks

  if (press_couple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all("Invalid fix press/berendsen command pressure settings");
  if (press_couple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all("Invalid fix press/berendsen command pressure settings");
  if (press_couple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all("Invalid fix press/berendsen command pressure settings");

  if (press_couple == XY && 
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1]))
    error->all("Invalid fix press/berendsen command pressure settings");
  if (press_couple == YZ && 
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2]))
    error->all("Invalid fix press/berendsen command pressure settings");
  if (press_couple == XZ && 
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2]))
    error->all("Invalid fix press/berendsen command pressure settings");

  if (p_flag[0] && domain->xperiodic == 0)
    error->all("Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all("Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all("Cannot use fix press/berendsen on a non-periodic dimension");

  if (p_flag[0] && p_period[0] <= 0.0)
    error->all("Fix press/berendsen period must be > 0.0");
  if (p_flag[1] && p_period[1] <= 0.0)
    error->all("Fix press/berendsen period must be > 0.0");
  if (p_flag[2] && p_period[2] <= 0.0)
    error->all("Fix press/berendsen period must be > 0.0");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  modify->add_compute(4,newarg);
  delete [] newarg;
  pflag = 1;

  nrigid = 0;
  rfix = NULL;
}

/* ---------------------------------------------------------------------- */

FixPressBerendsen::~FixPressBerendsen()
{
  delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixPressBerendsen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::init()
{
  if (domain->triclinic) 
    error->all("Cannot use fix press/berendsen with triclinic box");

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      int *dimflag = ((FixDeform *) modify->fix[i])->dimflag;
      if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) || 
	  (p_flag[2] && dimflag[2]))
	error->all("Cannot use fix press/berendsen and "
		   "fix deform on same dimension");
    }

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) 
    error->all("Temperature ID for fix press/berendsen does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  icompute = modify->find_compute(id_press);
  if (icompute < 0)
    error->all("Pressure ID for fix press/berendsen does not exist");
  pressure = modify->compute[icompute];

  // Kspace setting

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any rigid fixes exist so rigid bodies move when box is remapped
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts 
------------------------------------------------------------------------- */

void FixPressBerendsen::setup(int vflag)
{
  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::end_of_step()
{
  // compute new T,P

  if (press_couple == XYZ) {
    double tmp = temperature->compute_scalar();
    tmp = pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      dilation[i] = 
	pow(1.0 - update->dt/p_period[i] * 
	    (p_target[i]-p_current[i])/bulkmodulus,1.0/3.0);
    }
  }

  // remap simulation box and atoms
  // redo KSpace coeffs since volume has changed

  remap();
  if (kspace_flag) force->kspace->setup();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::couple()
{
  double *tensor = pressure->vector;

  if (press_couple == XYZ)
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  else if (press_couple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (press_couple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (press_couple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } if (press_couple == ANISO) {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixPressBerendsen::remap()
{
  int i;
  double oldlo,oldhi,ctr;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap) domain->x2lamda(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	domain->x2lamda(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo-ctr)*dilation[i] + ctr;
      domain->boxhi[i] = (oldhi-ctr)*dilation[i] + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap) domain->lamda2x(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ---------------------------------------------------------------------- */

int FixPressBerendsen::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all("Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all("Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning("Temperature for NPT is not for group all");

    // reset id_temp of pressure to new temperature ID
    
    icompute = modify->find_compute(id_press);
    if (icompute < 0) 
      error->all("Pressure ID for fix press/berendsen does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all("Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all("Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}
