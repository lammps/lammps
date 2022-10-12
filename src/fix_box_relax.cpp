// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "fix_box_relax.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_extra.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

#define MAX_LIFO_DEPTH 2     // 3 box0 arrays in *.h dimensioned to this

/* ---------------------------------------------------------------------- */

FixBoxRelax::FixBoxRelax(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_temp(nullptr), id_press(nullptr), tflag(0), pflag(0)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix box/relax", error);

  scalar_flag = 1;
  extscalar = 1;
  global_freq = 1;
  no_change_box = 1;

  // default values

  pcouple = NONE;
  allremap = 1;
  vmax = 0.0001;
  deviatoric_flag = 0;
  nreset_h0 = 0;

  p_target[0] = p_target[1] = p_target[2] =
    p_target[3] = p_target[4] = p_target[5] = 0.0;
  p_flag[0] = p_flag[1] = p_flag[2] =
    p_flag[3] = p_flag[4] = p_flag[5] = 0;

  // turn on tilt factor scaling, whenever applicable

  dimension = domain->dimension;

  scaleyz = scalexz = scalexy = 0;
  if (domain->yperiodic && domain->xy != 0.0) scalexy = 1;
  if (domain->zperiodic && dimension == 3) {
    if (domain->yz != 0.0) scaleyz = 1;
    if (domain->xz != 0.0) scalexz = 1;
  }

  // set fixed-point to default = center of cell

  fixedpoint[0] = 0.5*(domain->boxlo[0]+domain->boxhi[0]);
  fixedpoint[1] = 0.5*(domain->boxlo[1]+domain->boxhi[1]);
  fixedpoint[2] = 0.5*(domain->boxlo[2]+domain->boxhi[2]);

  // process keywords

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"iso") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax iso", error);
      pcouple = XYZ;
      p_target[0] = p_target[1] = p_target[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_target[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"aniso") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax aniso", error);
      pcouple = NONE;
      p_target[0] = p_target[1] = p_target[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_target[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax tri", error);
      pcouple = NONE;
      scalexy = scalexz = scaleyz = 0;
      p_target[0] = p_target[1] = p_target[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      p_target[3] = p_target[4] = p_target[5] = 0.0;
      p_flag[3] = p_flag[4] = p_flag[5] = 1;
      if (dimension == 2) {
        p_target[2] = p_target[3] = p_target[4] = 0.0;
        p_flag[2] = p_flag[3] = p_flag[4] = 0;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax x", error);
      p_target[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[0] = 1;
      deviatoric_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax y", error);
      p_target[1] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[1] = 1;
      deviatoric_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax z", error);
      p_target[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[2] = 1;
      deviatoric_flag = 1;
      iarg += 2;
      if (dimension == 2) error->all(FLERR,"Fix box/relax z not allowed for a 2d simulation");

    } else if (strcmp(arg[iarg],"yz") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax yz", error);
      p_target[3] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[3] = 1;
      deviatoric_flag = 1;
      scaleyz = 0;
      iarg += 2;
      if (dimension == 2) error->all(FLERR,"Fix box/relax yz not allowed for a 2d simulation");
    } else if (strcmp(arg[iarg],"xz") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax xz", error);
      p_target[4] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[4] = 1;
      deviatoric_flag = 1;
      scalexz = 0;
      iarg += 2;
      if (dimension == 2) error->all(FLERR,"Fix box/relax xz not allowed for a 2d simulation");
    } else if (strcmp(arg[iarg],"xy") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax xy", error);
      p_target[5] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_flag[5] = 1;
      deviatoric_flag = 1;
      scalexy = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax couple", error);
      if (strcmp(arg[iarg+1],"xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1],"xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1],"yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1],"xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1],"none") == 0) pcouple = NONE;
      else error->all(FLERR,"Illegal fix box/relax couple command: {} unknown", arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax dilate", error);
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1],"partial") == 0) allremap = 0;
      else error->all(FLERR,"Illegal fix box/relax dilate command: {} unkown", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vmax") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax vmax", error);
      vmax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (vmax <= 0.0) error->all(FLERR,"Fix box/relax vmax value {} must be > 0", vmax);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nreset") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax nreset", error);
      nreset_h0 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nreset_h0 < 0) error->all(FLERR,"Fix box/relax nreset value {} must be >= 0", nreset_h0);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scalexy") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax scalexy", error);
      scalexy = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scalexz") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax scalexz", error);
      scalexz = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scaleyz") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix box/relax scaleyz", error);
      scaleyz = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"fixedpoint") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix box/relax fixedpoint", error);
      fixedpoint[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      fixedpoint[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      fixedpoint[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else error->all(FLERR,"Unknown fix box/relax keyword {}", arg[iarg]);
  }

  if (p_flag[0]) box_change |= BOX_CHANGE_X;
  if (p_flag[1]) box_change |= BOX_CHANGE_Y;
  if (p_flag[2]) box_change |= BOX_CHANGE_Z;
  if (p_flag[3]) box_change |= BOX_CHANGE_YZ;
  if (p_flag[4]) box_change |= BOX_CHANGE_XZ;
  if (p_flag[5]) box_change |= BOX_CHANGE_XY;

  if (allremap == 0) restart_pbc = 1;

  // error checks

  if (dimension == 2 && (p_flag[2] || p_flag[3] || p_flag[4]))
    error->all(FLERR,"Invalid fix box/relax command for a 2d simulation");
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR,"Invalid fix box/relax command for a 2d simulation");

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix box/relax command pressure settings");
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR,"Invalid fix box/relax command pressure settings");
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix box/relax command pressure settings");
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix box/relax command pressure settings");
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix box/relax command pressure settings");

  // require periodicity in tensile dimension

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with non-periodic x dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with non-periodic y dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with non-periodic z dimension");

  // require periodicity in 2nd dim of off-diagonal tilt component

  if (p_flag[3] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix box/relax yz with non-periodic z dimension");
  if (p_flag[4] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix box/relax xz with non-periodic z dimension");
  if (p_flag[5] && domain->yperiodic == 0)
    error->all(FLERR, "Cannot use fix box/relax xz with non-periodic y dimension");

  if (scaleyz == 1 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with yz tilt factor scaling "
               "with non-periodic z dimension");
  if (scalexz == 1 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with xz tilt factor scaling "
               "with non-periodic z dimension");
  if (scalexy == 1 && domain->yperiodic == 0)
    error->all(FLERR,"Cannot use fix box/relax with xy tilt factor scaling "
               "with non-periodic y dimension");

  if (p_flag[3] && scaleyz == 1)
    error->all(FLERR,"Cannot use fix box/relax with both yz relaxation "
               "and scaling on yz tilt factor");
  if (p_flag[4] && scalexz == 1)
    error->all(FLERR,"Cannot use fix box/relax with both xy relaxation "
               "and scaling on xz tilt factor");
  if (p_flag[5] && scalexy == 1)
    error->all(FLERR,"Cannot use fix box/relax with both xy relaxation "
               "and scaling on xy tilt factor");

  if (!domain->triclinic && (p_flag[3] || p_flag[4] || p_flag[5]))
    error->all(FLERR,"Cannot specify Pxy/Pxz/Pyz in fix box/relax with non-triclinic box");

  if (pcouple == XYZ && dimension == 3 &&
      (p_target[0] != p_target[1] || p_target[0] != p_target[2]))
    error->all(FLERR,"Fix box/relax pressure for x, y, and z must be the same with couple xyz");
  if (pcouple == XYZ && dimension == 2 && p_target[0] != p_target[1])
    error->all(FLERR,"Fix box/relax pressure for x and y must be the same with couple xyz");
  if (pcouple == XY && p_target[0] != p_target[1])
    error->all(FLERR,"Fix box/relax pressure for x and y must be the same with couple xy");
  if (pcouple == YZ && p_target[1] != p_target[2])
    error->all(FLERR,"Fix box/relax pressure for y and z must be the same with couple yz");
  if (pcouple == XZ && p_target[0] != p_target[2])
    error->all(FLERR,"Fix box/relax pressure for x and z must be the same with couple xz");

  // pstyle = TRICLINIC if any off-diagonal term is controlled -> 6 dof
  // else pstyle = ISO if XYZ coupling or XY coupling in 2d -> 1 dof
  // else pstyle = ANISO -> 3 dof

  if (p_flag[3] || p_flag[4] || p_flag[5]) pstyle = TRICLINIC;
  else if (pcouple == XYZ || (dimension == 2 && pcouple == XY)) pstyle = ISO;
  else pstyle = ANISO;

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id) + "_temp");
  temperature = modify->add_compute(fmt::format("{} all temp",id_temp));
  tflag = 1;

  // create a new compute pressure style (virial only)
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id) + "_press");
  pressure = modify->add_compute(fmt::format("{} all pressure {} virial", id_press, id_temp));
  pflag = 1;

  dimension = domain->dimension;
  nrigid = 0;
  rfix = nullptr;

  current_lifo = 0;
}

/* ---------------------------------------------------------------------- */

FixBoxRelax::~FixBoxRelax()
{
  delete[] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete[] id_temp;
  delete[] id_press;
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

  temperature = modify->get_compute_by_id(id_temp);
  if (!temperature)
    error->all(FLERR,"Temperature compute ID {} for fix box/relax does not exist", id_temp);

  pressure = modify->get_compute_by_id(id_press);
  if (!pressure)
    error->all(FLERR,"Pressure compute ID {} for fix box/relax does not exist", id_press);

  pv2e = 1.0 / force->nktv2p;

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any rigid fixes exist so rigid bodies move when box is remapped
  // rfix[] = indices to each fix rigid

  delete[] rfix;
  nrigid = 0;
  rfix = nullptr;

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  }

  // initial box dimensions

  xprdinit = domain->xprd;
  yprdinit = domain->yprd;
  zprdinit = domain->zprd;
  if (dimension == 2) zprdinit = 1.0;
  vol0 = xprdinit * yprdinit * zprdinit;
  h0[0] = domain->h[0];
  h0[1] = domain->h[1];
  h0[2] = domain->h[2];
  h0[3] = domain->h[3];
  h0[4] = domain->h[4];
  h0[5] = domain->h[5];

  // hydrostatic target pressure and deviatoric target stress

  compute_press_target();
  if (deviatoric_flag) compute_sigma();
}

/* ----------------------------------------------------------------------
   compute energy and force due to extra degrees of freedom
------------------------------------------------------------------------- */

double FixBoxRelax::min_energy(double *fextra)
{
  double eng,scale,scalex,scaley,scalez,scalevol;

  temperature->compute_scalar();
  if (pstyle == ISO) pressure->compute_scalar();
  else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  // trigger virial computation on every iteration of minimizer

  pressure->addstep(update->ntimestep+1);

  // compute energy, forces for each extra degree of freedom
  // returned eng = PV must be in units of energy
  // returned fextra must likewise be in units of energy

  if (pstyle == ISO) {
    scale = domain->xprd/xprdinit;
    if (dimension == 3) {
      eng = pv2e * p_target[0] * (scale*scale*scale-1.0)*vol0;
      fextra[0] = pv2e * (p_current[0] - p_target[0])*3.0*scale*scale*vol0;
    } else {
      eng = pv2e * p_target[0] * (scale*scale-1.0)*vol0;
      fextra[0] = pv2e * (p_current[0] - p_target[0])*2.0*scale*vol0;
    }

  } else {
    fextra[0] = fextra[1] = fextra[2] = 0.0;
    scalex = scaley = scalez = 1.0;
    if (p_flag[0]) scalex = domain->xprd/xprdinit;
    if (p_flag[1]) scaley = domain->yprd/yprdinit;
    if (p_flag[2]) scalez = domain->zprd/zprdinit;
    scalevol = scalex*scaley*scalez;
    eng = pv2e * p_hydro * (scalevol-1.0)*vol0;
    if (p_flag[0])
      fextra[0] = pv2e * (p_current[0] - p_hydro)*scaley*scalez*vol0;
    if (p_flag[1])
      fextra[1] = pv2e * (p_current[1] - p_hydro)*scalex*scalez*vol0;
    if (p_flag[2])
      fextra[2] = pv2e * (p_current[2] - p_hydro)*scalex*scaley*vol0;

    if (pstyle == TRICLINIC) {
      fextra[3] = fextra[4] = fextra[5] = 0.0;
      if (p_flag[3])
        fextra[3] = pv2e*p_current[3]*scaley*yprdinit*scalex*xprdinit*yprdinit;
      if (p_flag[4])
        fextra[4] = pv2e*p_current[4]*scalex*xprdinit*scaley*yprdinit*xprdinit;
      if (p_flag[5])
        fextra[5] = pv2e*p_current[5]*scalex*xprdinit*scalez*zprdinit*xprdinit;
    }

    if (deviatoric_flag) {
      compute_deviatoric();
      if (p_flag[0]) fextra[0] -= fdev[0]*xprdinit;
      if (p_flag[1]) fextra[1] -= fdev[1]*yprdinit;
      if (p_flag[2]) fextra[2] -= fdev[2]*zprdinit;
      if (pstyle == TRICLINIC) {
        if (p_flag[3]) fextra[3] -= fdev[3]*yprdinit;
        if (p_flag[4]) fextra[4] -= fdev[4]*xprdinit;
        if (p_flag[5]) fextra[5] -= fdev[5]*xprdinit;
      }

      eng += compute_strain_energy();
    }
  }

  return eng;
}

/* ----------------------------------------------------------------------
   store extra dof values for minimization linesearch starting point
   boxlo0,boxhi0 = box dimensions
   box values are pushed onto a LIFO stack so nested calls can be made
   values are popped by calling min_step(0.0)
------------------------------------------------------------------------- */

void FixBoxRelax::min_store()
{
  for (int i = 0; i < 3; i++) {
    boxlo0[current_lifo][i] = domain->boxlo[i];
    boxhi0[current_lifo][i] = domain->boxhi[i];
  }
  if (pstyle == TRICLINIC) {
    boxtilt0[current_lifo][0] = domain->yz;
    boxtilt0[current_lifo][1] = domain->xz;
    boxtilt0[current_lifo][2] = domain->xy;
  }
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
    error->all(FLERR,"Attempt to push beyond stack limit in fix box/relax");
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
    error->all(FLERR,"Attempt to pop empty stack in fix box/relax");
    return;
  }
  current_lifo--;
}

/* ----------------------------------------------------------------------
   check if time to reset reference state. If so, do so.
------------------------------------------------------------------------- */

int FixBoxRelax::min_reset_ref()
{
  int itmp = 0;

  // if nreset_h0 > 0, reset reference box
  // every nreset_h0 timesteps
  // only needed for deviatoric external stress

  if (deviatoric_flag && nreset_h0 > 0) {
    int delta = update->ntimestep - update->beginstep;
    if (delta % nreset_h0 == 0) {
      compute_sigma();
      itmp = 1;
    }
  }
  return itmp;
}

/* ----------------------------------------------------------------------
   change the box dimensions by fraction ds = alpha*hextra
------------------------------------------------------------------------- */

void FixBoxRelax::min_step(double alpha, double *hextra)
{
  if (pstyle == ISO) {
    ds[0] = ds[1] = ds[2] = alpha*hextra[0];
  } else {
    ds[0] = ds[1] = ds[2] = 0.0;
    if (p_flag[0]) ds[0] = alpha*hextra[0];
    if (p_flag[1]) ds[1] = alpha*hextra[1];
    if (p_flag[2]) ds[2] = alpha*hextra[2];
    if (pstyle == TRICLINIC) {
      ds[3] = ds[4] = ds[5] = 0.0;
      if (p_flag[3]) ds[3] = alpha*hextra[3];
      if (p_flag[4]) ds[4] = alpha*hextra[4];
      if (p_flag[5]) ds[5] = alpha*hextra[5];
    }
  }
  remap();
  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   max allowed step size along hextra
------------------------------------------------------------------------- */

double FixBoxRelax::max_alpha(double *hextra)
{
  double alpha = 1.0;
  if (pstyle == ISO) alpha = vmax/fabs(hextra[0]);
  else {
    if (p_flag[0]) alpha = MIN(alpha,vmax/fabs(hextra[0]));
    if (p_flag[1]) alpha = MIN(alpha,vmax/fabs(hextra[1]));
    if (p_flag[2]) alpha = MIN(alpha,vmax/fabs(hextra[2]));
    if (pstyle == TRICLINIC) {
      if (p_flag[3]) alpha = MIN(alpha,vmax/fabs(hextra[3]));
      if (p_flag[4]) alpha = MIN(alpha,vmax/fabs(hextra[4]));
      if (p_flag[5]) alpha = MIN(alpha,vmax/fabs(hextra[5]));
    }
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   return number of degrees of freedom added by this fix
------------------------------------------------------------------------- */

int FixBoxRelax::min_dof()
{
  if (pstyle == ISO) return 1;
  if (pstyle == TRICLINIC) return 6;
  return 3;
}

/* ----------------------------------------------------------------------
   dilate the box and owned/ghost atoms around center of box
------------------------------------------------------------------------- */

void FixBoxRelax::remap()
{
  int i,n;

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
      domain->boxlo[i] = currentBoxLo0 + (currentBoxLo0 - fixedpoint[i])/domain->h[i]*ds[i]*h0[i];
      domain->boxhi[i] = currentBoxHi0 + (currentBoxHi0 - fixedpoint[i])/domain->h[i]*ds[i]*h0[i];
      if (domain->boxlo[i] >= domain->boxhi[i])
        error->all(FLERR,"Fix box/relax generated negative box length");
    }

  // scale tilt factors with cell, if set

  if (scaleyz) domain->yz = (domain->boxhi[2] - domain->boxlo[2])*h0[3]/h0[2];
  if (scalexz) domain->xz = (domain->boxhi[2] - domain->boxlo[2])*h0[4]/h0[2];
  if (scalexy) domain->xy = (domain->boxhi[1] - domain->boxlo[1])*h0[5]/h0[1];

  if (pstyle == TRICLINIC) {
    if (p_flag[3]) domain->yz = boxtilt0[current_lifo][0]+ds[3]*yprdinit;
    if (p_flag[4]) domain->xz = boxtilt0[current_lifo][1]+ds[4]*xprdinit;
    if (p_flag[5]) domain->xy = boxtilt0[current_lifo][2]+ds[5]*xprdinit;
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

  if (pstyle == ISO)
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  else if (pcouple == XYZ) {
    double ave = 1.0/3.0 * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }

  if (!std::isfinite(p_current[0]) || !std::isfinite(p_current[1]) || !std::isfinite(p_current[2]))
    error->all(FLERR,"Non-numeric pressure - simulation unstable");

  // switch order from xy-xz-yz to Voigt ordering

  if (pstyle == TRICLINIC) {
    p_current[3] = tensor[5];
    p_current[4] = tensor[4];
    p_current[5] = tensor[3];

    if (!std::isfinite(p_current[3]) || !std::isfinite(p_current[4]) || !std::isfinite(p_current[5]))
      error->all(FLERR,"Non-numeric pressure - simulation unstable");
  }
}

/* ---------------------------------------------------------------------- */

int FixBoxRelax::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR,"Could not find fix_modify {} temperature compute ID {}", id, id_temp);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify {} temperature compute {} does not compute temperature",
                 id, id_temp);
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature compute {} for fix modify {} is not for group all",
                     id, id_temp);

    // reset id_temp of pressure to new temperature ID

    pressure = modify->get_compute_by_id(id_press);
    if (!pressure)
      error->all(FLERR,"Pressure compute ID {} for fix modify {} does not exist", id_press, id);
    pressure->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete[] id_press;
    id_press = utils::strdup(arg[1]);

    pressure = modify->get_compute_by_id(id_press);
    if (!pressure)
      error->all(FLERR,"Could not find fix_modify {} compute pressure ID {}", id, id_press);

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify {} pressure compute {} does not compute pressure", id, id_press);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   compute sigma tensor (needed whenever reference box is reset)
-----------------------------------------------------------------------*/

void FixBoxRelax::compute_sigma()
{
  double pdeviatoric[3][3];
  double tmp1[3][3],sigma_tensor[3][3],h_invtmp[3][3];

  // reset reference box dimensions

  xprdinit = domain->xprd;
  yprdinit = domain->yprd;
  zprdinit = domain->zprd;
  if (dimension == 2) zprdinit = 1.0;
  vol0 = xprdinit * yprdinit * zprdinit;

  h0_inv[0] = domain->h_inv[0];
  h0_inv[1] = domain->h_inv[1];
  h0_inv[2] = domain->h_inv[2];
  h0_inv[3] = domain->h_inv[3];
  h0_inv[4] = domain->h_inv[4];
  h0_inv[5] = domain->h_inv[5];

  h_invtmp[0][0] = h0_inv[0];
  h_invtmp[1][1] = h0_inv[1];
  h_invtmp[2][2] = h0_inv[2];
  h_invtmp[1][2] = h0_inv[3];
  h_invtmp[0][2] = h0_inv[4];
  h_invtmp[0][1] = h0_inv[5];
  h_invtmp[2][1] = 0.0;
  h_invtmp[2][0] = 0.0;
  h_invtmp[1][0] = 0.0;

  // compute target deviatoric stress tensor pdevmod

  pdeviatoric[0][0] = pdeviatoric[1][1] = pdeviatoric[2][2] = 0.0;
  if (p_flag[0]) pdeviatoric[0][0] = p_target[0] - p_hydro;
  if (p_flag[1]) pdeviatoric[1][1] = p_target[1] - p_hydro;
  if (p_flag[2]) pdeviatoric[2][2] = p_target[2] - p_hydro;
  pdeviatoric[1][2] = pdeviatoric[2][1] = p_target[3];
  pdeviatoric[0][2] = pdeviatoric[2][0] = p_target[4];
  pdeviatoric[0][1] = pdeviatoric[1][0] = p_target[5];

  // Modify to account for off-diagonal terms
  // These equations come from the stationarity relation:
  //    Pdev,sys = Pdev,targ*hinv^t*hdiag
  // where:
  // Pdev,sys is the system deviatoric stress tensor,
  // Pdev,targ = pdeviatoric, effective target deviatoric stress
  // hinv^t is the transpose of the inverse h tensor
  // hdiag is the diagonal part of the h tensor

  pdeviatoric[1][1] -= pdeviatoric[1][2]*h0_inv[3]*h0[1];
  pdeviatoric[0][1] -= pdeviatoric[0][2]*h0_inv[3]*h0[1];
  pdeviatoric[1][0] = pdeviatoric[0][1];
  pdeviatoric[0][0] -= pdeviatoric[0][1]*h0_inv[5]*h0[0] +
    pdeviatoric[0][2]*h0_inv[4]*h0[0];

  // compute symmetric sigma tensor

  MathExtra::times3(h_invtmp,pdeviatoric,tmp1);
  MathExtra::times3_transpose(tmp1,h_invtmp,sigma_tensor);
  MathExtra::scalar_times3(vol0,sigma_tensor);

  sigma[0] = sigma_tensor[0][0];
  sigma[1] = sigma_tensor[1][1];
  sigma[2] = sigma_tensor[2][2];
  sigma[3] = sigma_tensor[1][2];
  sigma[4] = sigma_tensor[0][2];
  sigma[5] = sigma_tensor[0][1];
}

/* ----------------------------------------------------------------------
   compute strain energy
-----------------------------------------------------------------------*/

double FixBoxRelax::compute_strain_energy()
{
  // compute strain energy = 0.5*Tr(sigma*h*h^t) in energy units

  double* h = domain->h;
  double d0,d1,d2;

  if (dimension == 3) {
    d0 =
      sigma[0]*(h[0]*h[0]+h[5]*h[5]+h[4]*h[4]) +
      sigma[5]*(          h[1]*h[5]+h[3]*h[4]) +
      sigma[4]*(                    h[2]*h[4]);
    d1 =
      sigma[5]*(          h[5]*h[1]+h[4]*h[3]) +
      sigma[1]*(          h[1]*h[1]+h[3]*h[3]) +
      sigma[3]*(                    h[2]*h[3]);
    d2 =
      sigma[4]*(                    h[4]*h[2]) +
      sigma[3]*(                    h[3]*h[2]) +
      sigma[2]*(                    h[2]*h[2]);
  } else {
    d0 = sigma[0]*(h[0]*h[0]+h[5]*h[5]) + sigma[5]*h[1]*h[5];
    d1 = sigma[5]*h[5]*h[1] + sigma[1]*h[1]*h[1];
    d2 = 0.0;
  }

  double energy = 0.5*(d0+d1+d2)*pv2e;
  return energy;
}

/* ----------------------------------------------------------------------
   compute deviatoric barostat force = h*sigma*h^t
-----------------------------------------------------------------------*/

void FixBoxRelax::compute_deviatoric()
{
  double* h = domain->h;

  // [ 0 5 4 ]   [ 0 5 4 ] [ 0 5 4 ]
  // [ 5 1 3 ] = [ - 1 3 ] [ 5 1 3 ]
  // [ 4 3 2 ]   [ - - 2 ] [ 4 3 2 ]

  if (dimension == 3) {
    fdev[0] = pv2e*(h[0]*sigma[0]+h[5]*sigma[5]+h[4]*sigma[4]);
    fdev[1] = pv2e*(h[1]*sigma[1]+h[3]*sigma[3]);
    fdev[2] = pv2e*(h[2]*sigma[2]);
    fdev[3] = pv2e*(h[1]*sigma[3]+h[3]*sigma[2]);
    fdev[4] = pv2e*(h[0]*sigma[4]+h[5]*sigma[3]+h[4]*sigma[2]);
    fdev[5] = pv2e*(h[0]*sigma[5]+h[5]*sigma[1]+h[4]*sigma[3]);
  } else {
    fdev[0] = pv2e*(h[0]*sigma[0]+h[5]*sigma[5]);
    fdev[1] = pv2e*(h[1]*sigma[1]);
    fdev[5] = pv2e*(h[0]*sigma[5]+h[5]*sigma[1]);
  }
}

/* ----------------------------------------------------------------------
   compute hydrostatic target pressure
-----------------------------------------------------------------------*/

void FixBoxRelax::compute_press_target()
{
  pflagsum = p_flag[0] + p_flag[1] + p_flag[2];

  p_hydro = 0.0;
  for (int i = 0; i < 3; i++)
    if (p_flag[i]) p_hydro += p_target[i];
  if (pflagsum) p_hydro /= pflagsum;

  for (int i = 0; i < 3; i++) {
    if (p_flag[i] && fabs(p_hydro - p_target[i]) > 1.0e-6) deviatoric_flag = 1;
  }

  if (pstyle == TRICLINIC) {
    for (int i = 3; i < 6; i++)
      if (p_flag[i] && fabs(p_target[i]) > 1.0e-6) deviatoric_flag = 1;
  }
}

/* ----------------------------------------------------------------------
   compute PV and strain energy for access to the user
   ---------------------------------------------------------------------- */

double FixBoxRelax::compute_scalar()
{
  double ftmp[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  if (update->ntimestep == 0) return 0.0;
  return min_energy(ftmp);
}
