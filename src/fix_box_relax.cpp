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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_box_relax.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

enum{XYZ,XY,YZ,XZ,ANISO};

/* ---------------------------------------------------------------------- */

FixBoxRelax::FixBoxRelax(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix box/relax command");

  box_change = 1;

  if (strcmp(arg[3],"xyz") == 0) {
    if (narg < 5) error->all("Illegal fix box/relax command");
    press_couple = XYZ;
    p_target[0] = p_target[1] = p_target[2] = atof(arg[4]);
    p_flag[0] = p_flag[1] = p_flag[2] = 1;
    if (domain->dimension == 2) p_flag[2] = 0;

  } else {
    if (strcmp(arg[3],"xy") == 0) press_couple = XY;
    else if (strcmp(arg[3],"yz") == 0) press_couple = YZ;
    else if (strcmp(arg[3],"xz") == 0) press_couple = XZ;
    else if (strcmp(arg[3],"aniso") == 0) press_couple = ANISO;
    else error->all("Illegal fix box/relax command");

    if (narg < 7) error->all("Illegal fix box/relax command");

    if (domain->dimension == 2 && 
	(press_couple == XY || press_couple == YZ || press_couple == XZ))
      error->all("Invalid fix box/relax command for a 2d simulation");

    if (strcmp(arg[4],"NULL") == 0) {
      p_target[0] = 0.0;
      p_flag[0] = 0;
    } else {
      p_target[0] = atof(arg[4]);
      p_flag[0] = 1;
    }
    if (strcmp(arg[5],"NULL") == 0) {
      p_target[1] = 0.0;
      p_flag[1] = 0;
    } else {
      p_target[1] = atof(arg[5]);
      p_flag[1] = 1;
    }
    if (strcmp(arg[6],"NULL") == 0) {
      p_target[2] = 0.0;
      p_flag[2] = 0;
    } else {
      if (domain->dimension == 2)
	error->all("Invalid fix box/relax command for a 2d simulation");
      p_target[2] = atof(arg[6]);
      p_flag[2] = 1;
    }
  }

  pflagsum = p_flag[0] + p_flag[1] + p_flag[2];

  // process extra keywords

  allremap = 1;
  vmax = 0.0001;

  int iarg;
  if (press_couple == XYZ) iarg = 5;
  else iarg = 7;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix box/relax command");
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1],"partial") == 0) allremap = 0;
      else error->all("Illegal fix box/relax command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vmax") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix box/relax command");
      vmax = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal fix box/relax command");
  }

  // error checks

  if (press_couple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all("Illegal fix box/relax command");
  if (press_couple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all("Illegal fix box/relax command");
  if (press_couple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all("Illegal fix box/relax command");

  if (press_couple == XY && p_target[0] != p_target[1])
    error->all("Illegal fix box/relax command");
  if (press_couple == YZ && p_target[1] != p_target[2])
    error->all("Illegal fix box/relax command");
  if (press_couple == XZ && p_target[0] != p_target[2])
    error->all("Illegal fix box/relax command");

  if (p_flag[0] && domain->xperiodic == 0)
    error->all("Cannot use fix box/relax on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all("Cannot use fix box/relax on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all("Cannot use fix box/relax on a non-periodic dimension");

  if (vmax <= 0.0) error->all("Illegal fix box/relax command");

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

  // create a new compute pressure style (virial only)
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[5];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  newarg[4] = (char *) "virial";
  modify->add_compute(5,newarg);
  delete [] newarg;
  pflag = 1;

  dimension = domain->dimension;
  nrigid = 0;
  rfix = 0;

  current_lifo = 0;
}

/* ---------------------------------------------------------------------- */

FixBoxRelax::~FixBoxRelax()
{
  delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixBoxRelax::setmask()
{
  int mask = 0;
  mask |= MIN_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBoxRelax::init()
{
  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) 
    error->all("Temperature ID for fix box/relax does not exist");
  temperature = modify->compute[icompute];

  icompute = modify->find_compute(id_press);
  if (icompute < 0) error->all("Pressure ID for fix box/relax does not exist");
  pressure = modify->compute[icompute];

  // initial box dimensions

  xprdinit = domain->xprd;
  yprdinit = domain->yprd;
  zprdinit = domain->zprd;
  if (dimension == 3) volinit = domain->xprd*domain->yprd*domain->zprd;
  else volinit = domain->xprd*domain->yprd;
  pv2e = 1.0 / force->nktv2p;

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
   compute energy and force due to extra degrees of freedom
------------------------------------------------------------------------- */

double FixBoxRelax::min_energy(double *fextra)
{
  double eng,scale,scalex,scaley,scalez;

  double t_current = temperature->compute_scalar();
  if (press_couple == XYZ) {
    double tmp = pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  // trigger virial computation on every iteration of minimizer

  pressure->addstep(update->ntimestep+1);

  // compute energy, forces for each extra degree of freedom
  // returned eng = PV must be in units of energy
  // returned fextra must likewise be in units of energy

  if (press_couple == XYZ) {
    scale = domain->xprd/xprdinit;
    if (dimension == 3) {
      eng = pv2e * p_target[0] * (scale*scale*scale-1.0)*volinit;
      fextra[0] = pv2e * (p_current[0] - p_target[0])*3.0*scale*scale*volinit;
    } else {
      eng = pv2e * p_target[0] * (scale*scale-1.0)*volinit;
      fextra[0] = pv2e * (p_current[0] - p_target[0])*2.0*scale*volinit;
    }

  } else {
    fextra[0] = fextra[1] = fextra[2] = 0.0;
    scalex = scaley = scalez = 1.0;
    if (p_flag[0]) scalex = domain->xprd/xprdinit;
    if (p_flag[1]) scaley = domain->yprd/yprdinit;
    if (p_flag[2]) scalez = domain->zprd/zprdinit;
    eng = pv2e * (p_flag[0]*p_target[0] + p_flag[1]*p_target[1] + 
		  p_flag[2]*p_target[2])/pflagsum * 
      (scalex*scaley*scalez-1.0)*volinit;
    if (p_flag[0])
      fextra[0] = pv2e * (p_current[0] - p_target[0])*scaley*scalez*volinit;
    if (p_flag[1])
      fextra[1] = pv2e * (p_current[1] - p_target[1])*scalex*scalez*volinit;
    if (p_flag[2])
      fextra[2] = pv2e * (p_current[2] - p_target[2])*scalex*scaley*volinit;
  }

  return eng;
}

/* ----------------------------------------------------------------------
   store extra dof values for minimization linesearch starting point
   boxlo0,boxhi0 = box dimensions
   s0 = ratio of current boxsize to initial boxsize
   box values are pushed onto a LIFO stack so nested calls can be made
   values are popped by calling min_step(0.0)
------------------------------------------------------------------------- */

void FixBoxRelax::min_store()
{
  for (int i = 0; i < 3; i++) {
    boxlo0[current_lifo][i] = domain->boxlo[i];
    boxhi0[current_lifo][i] = domain->boxhi[i];
  }
  s0[0] = (boxhi0[current_lifo][0]-boxlo0[current_lifo][0])/xprdinit;
  s0[1] = (boxhi0[current_lifo][1]-boxlo0[current_lifo][1])/yprdinit;
  s0[2] = (boxhi0[current_lifo][2]-boxlo0[current_lifo][2])/zprdinit;
}

/* ----------------------------------------------------------------------
   clear the LIFO stack for min_store
------------------------------------------------------------------------- */

void FixBoxRelax::min_clearstore()
{
  current_lifo = 0;
}

/* ----------------------------------------------------------------------
   push the LIFO stack for min_store
------------------------------------------------------------------------- */

void FixBoxRelax::min_pushstore()
{
  if (current_lifo >= MAX_LIFO_DEPTH) {
    error->all("Attempt to push beyond stack limit <FixBoxRelax>");
    return;
  }

  current_lifo++;
}


/* ----------------------------------------------------------------------
   pop the LIFO stack for min_store
------------------------------------------------------------------------- */

void FixBoxRelax::min_popstore()
{
  if (current_lifo <= 0) {
    error->all("Attempt to pop empty stack <FixBoxRelax>");
    return;
  }

  current_lifo--;
}

/* ----------------------------------------------------------------------
   change the box dimensions by fraction ds = alpha*hextra
------------------------------------------------------------------------- */

void FixBoxRelax::min_step(double alpha, double *hextra)
{
  if (press_couple == XYZ) {
    ds[0] = ds[1] = ds[2] = alpha*hextra[0];
  } else {
    if (p_flag[0]) ds[0] = alpha*hextra[0];
    if (p_flag[1]) ds[1] = alpha*hextra[1];
    if (p_flag[2]) ds[2] = alpha*hextra[2];
  }

  remap();
}

/* ----------------------------------------------------------------------
   max allowed step size along hextra
------------------------------------------------------------------------- */

double FixBoxRelax::max_alpha(double *hextra)
{
  double alpha = 0.0;
  if (press_couple == XYZ) alpha = vmax/fabs(hextra[0]);
  else {
    alpha = vmax/fabs(hextra[0]);
    alpha = MIN(alpha,vmax/fabs(hextra[1]));
    alpha = MIN(alpha,vmax/fabs(hextra[2]));
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   return number of degrees of freedom added by this fix
------------------------------------------------------------------------- */

int FixBoxRelax::min_dof()
{
  if (press_couple == XYZ) return 1;
  return 3;
}

/* ----------------------------------------------------------------------
   dilate the box and owned/ghost atoms around center of box
------------------------------------------------------------------------- */

void FixBoxRelax::remap()
{
  int i,n;
  double ctr;
  
  // ctr = geometric center of box in a dimension
  // rescale simulation box from linesearch starting point
  // scale atom coords for all atoms or only for fix group atoms

  double **x = atom->x;
  int *mask = atom->mask;
  n = atom->nlocal + atom->nghost;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap) domain->x2lamda(n);
  else {
    for (i = 0; i < n; i++)
      if (mask[i] & groupbit)
	domain->x2lamda(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++)
    if (p_flag[i]) {
      double currentBoxLo0 = boxlo0[current_lifo][i];
      double currentBoxHi0 = boxhi0[current_lifo][i];
      ctr = 0.5 * (currentBoxLo0 + currentBoxHi0);
      domain->boxlo[i] = currentBoxLo0 + (currentBoxLo0-ctr)*ds[i]/s0[i];
      domain->boxhi[i] = currentBoxHi0 + (currentBoxHi0-ctr)*ds[i]/s0[i];
    }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap) domain->lamda2x(n);
  else {
    for (i = 0; i < n; i++)
      if (mask[i] & groupbit)
	domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ---------------------------------------------------------------------- */

void FixBoxRelax::couple()
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
  } else if (press_couple == ANISO) {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ---------------------------------------------------------------------- */

int FixBoxRelax::modify_param(int narg, char **arg)
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
      error->warning("Temperature for fix modify is not for group all");

    // reset id_temp of pressure to new temperature ID
    
    icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all("Pressure ID for fix modify does not exist");
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
