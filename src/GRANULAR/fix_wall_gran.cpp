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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL),
                         Dan Bolintineanu (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// XYZ PLANE need to be 0,1,2

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,GRANULAR};
enum{NONE,CONSTANT,EQUAL};

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS 1.333333333333333       // 4/3
#define THREEQUARTERS 0.75                 // 3/4
#define TWOPI 6.28318530717959             // 2*PI

#define EPSILON 1e-10

enum {NORMAL_HOOKE, NORMAL_HERTZ, HERTZ_MATERIAL, DMT, JKR};
enum {VELOCITY, VISCOELASTIC, TSUJI};
enum {TANGENTIAL_NOHISTORY, TANGENTIAL_HISTORY,
      TANGENTIAL_MINDLIN, TANGENTIAL_MINDLIN_RESCALE};
enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
enum {ROLL_NONE, ROLL_SDS};

#define BIG 1.0e20
#define EPSILON 1e-10

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(NULL), history_one(NULL),
  fix_rigid(NULL), mass_rigid(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix wall/gran command");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix wall/gran requires atom style sphere");

  create_attribute = 1;

  // set interaction style
  // disable bonded/history option for now

  if (strcmp(arg[3],"hooke") == 0) pairstyle = HOOKE;
  else if (strcmp(arg[3],"hooke/history") == 0) pairstyle = HOOKE_HISTORY;
  else if (strcmp(arg[3],"hertz/history") == 0) pairstyle = HERTZ_HISTORY;
  else if (strcmp(arg[3],"granular") == 0) pairstyle = GRANULAR;
  else error->all(FLERR,"Invalid fix wall/gran interaction style");

  use_history = restart_peratom = 1;
  if (pairstyle == HOOKE) use_history = restart_peratom = 0;

  // wall/particle coefficients

  int iarg;

  if (pairstyle != GRANULAR) {
    size_history = 3;
    if (narg < 11) error->all(FLERR,"Illegal fix wall/gran command");

    kn = force->numeric(FLERR,arg[4]);
    if (strcmp(arg[5],"NULL") == 0) kt = kn * 2.0/7.0;
    else kt = force->numeric(FLERR,arg[5]);

    gamman = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) gammat = 0.5 * gamman;
    else gammat = force->numeric(FLERR,arg[7]);

    xmu = force->numeric(FLERR,arg[8]);
    int dampflag = force->inumeric(FLERR,arg[9]);
    if (dampflag == 0) gammat = 0.0;

    if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
        xmu < 0.0 || xmu > 10000.0 || dampflag < 0 || dampflag > 1)
      error->all(FLERR,"Illegal fix wall/gran command");

    // convert Kn and Kt from pressure units to force/distance^2 if Hertzian

    if (pairstyle == HERTZ_HISTORY) {
      kn /= force->nktv2p;
      kt /= force->nktv2p;
    }
    iarg = 10;
  } else {
    iarg = 4;
    damping_model = VISCOELASTIC;
    roll_model = twist_model = NONE;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "hooke") == 0) {
        if (iarg + 2 >= narg) 
          error->all(FLERR,"Illegal fix wall/gran command, "
                     "not enough parameters provided for Hooke option");
        normal_model = NORMAL_HOOKE;
        normal_coeffs[0] = force->numeric(FLERR,arg[iarg+1]); //kn
        normal_coeffs[1] = force->numeric(FLERR,arg[iarg+2]); //damping
        iarg += 3;
      } else if (strcmp(arg[iarg], "hertz") == 0) {
        int num_coeffs = 2;
        if (iarg + num_coeffs >= narg) 
          error->all(FLERR,"Illegal fix wall/gran command, "
                     "not enough parameters provided for Hertz option");
        normal_model = NORMAL_HERTZ;
        normal_coeffs[0] = force->numeric(FLERR,arg[iarg+1]); //kn
        normal_coeffs[1] = force->numeric(FLERR,arg[iarg+2]); //damping
        iarg += num_coeffs+1;
      } else if (strcmp(arg[iarg], "hertz/material") == 0) {
        int num_coeffs = 3;
        if (iarg + num_coeffs >= narg) 
          error->all(FLERR,"Illegal fix wall/gran command, "
                     "not enough parameters provided for Hertz option");
        normal_model = HERTZ_MATERIAL;
        Emod = force->numeric(FLERR,arg[iarg+1]); //E
        normal_coeffs[1] = force->numeric(FLERR,arg[iarg+2]); //damping
        poiss = force->numeric(FLERR,arg[iarg+3]); //Poisson's ratio
        normal_coeffs[0] = Emod/(2*(1-poiss))*FOURTHIRDS;
        normal_coeffs[2] = poiss;
        iarg += num_coeffs+1;
      } else if (strcmp(arg[iarg], "dmt") == 0) {
        if (iarg + 4 >= narg) 
          error->all(FLERR,"Illegal fix wall/gran command, "
                     "not enough parameters provided for Hertz option");
        normal_model = DMT;
        Emod = force->numeric(FLERR,arg[iarg+1]); //E
        normal_coeffs[1] = force->numeric(FLERR,arg[iarg+2]); //damping
        poiss = force->numeric(FLERR,arg[iarg+3]); //Poisson's ratio
        normal_coeffs[0] = Emod/(2*(1-poiss))*FOURTHIRDS;
        normal_coeffs[2] = poiss;
        normal_coeffs[3] = force->numeric(FLERR,arg[iarg+4]); //cohesion
        iarg += 5;
      } else if (strcmp(arg[iarg], "jkr") == 0) {
        if (iarg + 4 >= narg) 
          error->all(FLERR,"Illegal wall/gran command, "
                     "not enough parameters provided for JKR option");
        beyond_contact = 1;
        normal_model = JKR;
        Emod = force->numeric(FLERR,arg[iarg+1]); //E
        normal_coeffs[1] = force->numeric(FLERR,arg[iarg+2]); //damping
        poiss = force->numeric(FLERR,arg[iarg+3]); //Poisson's ratio
        normal_coeffs[0] = Emod/(2*(1-poiss))*FOURTHIRDS;
        normal_coeffs[2] = poiss;
        normal_coeffs[3] = force->numeric(FLERR,arg[iarg+4]); //cohesion
        iarg += 5;
      } else if (strcmp(arg[iarg], "damping") == 0) {
        if (iarg+1 >= narg) 
          error->all(FLERR, "Illegal wall/gran command, "
                     "not enough parameters provided for damping model");
        if (strcmp(arg[iarg+1], "velocity") == 0) {
          damping_model = VELOCITY;
          iarg += 1;
        } else if (strcmp(arg[iarg+1], "viscoelastic") == 0) {
          damping_model = VISCOELASTIC;
          iarg += 1;
        } else if (strcmp(arg[iarg+1], "tsuji") == 0) {
          damping_model = TSUJI;
          iarg += 1;
        } else error->all(FLERR, "Illegal wall/gran command, "
                          "unrecognized damping model");
        iarg += 1;
      } else if (strcmp(arg[iarg], "tangential") == 0) {
        if (iarg + 1 >= narg) 
          error->all(FLERR,"Illegal pair_coeff command, "
                     "must specify tangential model after tangential keyword");
        if (strcmp(arg[iarg+1], "linear_nohistory") == 0) {
          if (iarg + 3 >= narg) 
            error->all(FLERR,"Illegal pair_coeff command, "
                       "not enough parameters provided for tangential model");
          tangential_model = TANGENTIAL_NOHISTORY;
          tangential_coeffs[0] = 0;
          // gammat and friction coeff
          tangential_coeffs[1] = force->numeric(FLERR,arg[iarg+2]);
          tangential_coeffs[2] = force->numeric(FLERR,arg[iarg+3]);
          iarg += 4;
        } else if ((strcmp(arg[iarg+1], "linear_history") == 0) ||
            (strcmp(arg[iarg+1], "mindlin") == 0) ||
            (strcmp(arg[iarg+1], "mindlin_rescale") == 0)) {
          if (iarg + 4 >= narg) 
            error->all(FLERR,"Illegal pair_coeff command, "
                       "not enough parameters provided for tangential model");
          if (strcmp(arg[iarg+1], "linear_history") == 0) 
            tangential_model = TANGENTIAL_HISTORY;
          else if (strcmp(arg[iarg+1], "mindlin") == 0) 
            tangential_model = TANGENTIAL_MINDLIN;
          else if (strcmp(arg[iarg+1], "mindlin_rescale") == 0) 
            tangential_model = TANGENTIAL_MINDLIN_RESCALE;
          if ((tangential_model == TANGENTIAL_MINDLIN || 
               tangential_model == TANGENTIAL_MINDLIN_RESCALE) &&
              (strcmp(arg[iarg+2], "NULL") == 0)) {
            if (normal_model == NORMAL_HERTZ || normal_model == NORMAL_HOOKE) {
              error->all(FLERR, "NULL setting for Mindlin tangential "
                         "stiffness requires a normal contact model "
                         "that specifies material properties");
            }
            tangential_coeffs[0] = 4*(2-poiss)*(1+poiss)/Emod;
          } else {
            tangential_coeffs[0] = force->numeric(FLERR,arg[iarg+2]); //kt
          }
          tangential_history = 1;
          // gammat and friction coeff
          tangential_coeffs[1] = force->numeric(FLERR,arg[iarg+3]);
          tangential_coeffs[2] = force->numeric(FLERR,arg[iarg+4]);
          iarg += 5;
        } else {
          error->all(FLERR, "Illegal pair_coeff command, "
                     "tangential model not recognized");
        }
      } else if (strcmp(arg[iarg], "rolling") == 0) {
        if (iarg + 1 >= narg) 
          error->all(FLERR, "Illegal wall/gran command, not enough parameters");
        if (strcmp(arg[iarg+1], "none") == 0) {
          roll_model = ROLL_NONE;
          iarg += 2;
        } else if (strcmp(arg[iarg+1], "sds") == 0) {
          if (iarg + 4 >= narg) 
            error->all(FLERR,"Illegal wall/gran command, "
                       "not enough parameters provided for rolling model");
          roll_model = ROLL_SDS;
          roll_history = 1;
          // kR, gammaR, rolling friction coeff
          roll_coeffs[0] = force->numeric(FLERR,arg[iarg+2]);
          roll_coeffs[1] = force->numeric(FLERR,arg[iarg+3]);
          roll_coeffs[2] = force->numeric(FLERR,arg[iarg+4]);
          iarg += 5;
        } else {
          error->all(FLERR, "Illegal wall/gran command, "
                     "rolling friction model not recognized");
        }
      } else if (strcmp(arg[iarg], "twisting") == 0) {
        if (iarg + 1 >= narg)
          error->all(FLERR, "Illegal wall/gran command, not enough parameters");
        if (strcmp(arg[iarg+1], "none") == 0) {
          twist_model = TWIST_NONE;
          iarg += 2;
        } else if (strcmp(arg[iarg+1], "marshall") == 0) {
          twist_model = TWIST_MARSHALL;
          twist_history = 1;
          iarg += 2;
        } else if (strcmp(arg[iarg+1], "sds") == 0) {
          if (iarg + 4 >= narg) 
            error->all(FLERR,"Illegal wall/gran command, "
                       "not enough parameters provided for twist model");
          twist_model = TWIST_SDS;
          twist_history = 1;
          twist_coeffs[0] = force->numeric(FLERR,arg[iarg+2]); //kt
          twist_coeffs[1] = force->numeric(FLERR,arg[iarg+3]); //gammat
          twist_coeffs[2] = force->numeric(FLERR,arg[iarg+4]); //friction coeff.
          iarg += 5;
        } else {
          error->all(FLERR, "Illegal wall/gran command, "
                     "twisting friction model not recognized");
        }
      } else if (strcmp(arg[iarg], "xplane") == 0 ||
          strcmp(arg[iarg], "yplane") == 0 ||
          strcmp(arg[iarg], "zplane") == 0 ||
          strcmp(arg[iarg], "zcylinder") == 0 ||
          strcmp(arg[iarg], "region") == 0) {
        break;
      } else {
        error->all(FLERR, "Illegal fix wall/gran command");
      }
    }
    size_history = 3*tangential_history + 3*roll_history + twist_history;
    if (normal_model == JKR) size_history += 1;
    if (tangential_model == TANGENTIAL_MINDLIN_RESCALE) size_history += 1;
  }

  // wallstyle args

  idregion = NULL;

  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
  } else if (strcmp(arg[iarg],"region") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = REGION;
    int n = strlen(arg[iarg+1]) + 1;
    idregion = new char[n];
    strcpy(idregion,arg[iarg+1]);
    iarg += 2;
  }

  // optional args

  wiggle = 0;
  wshear = 0;
  peratom_flag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      amplitude = force->numeric(FLERR,arg[iarg+2]);
      period = force->numeric(FLERR,arg[iarg+3]);
      wiggle = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"shear") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      vshear = force->numeric(FLERR,arg[iarg+2]);
      wshear = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"store_contacts") == 0) {
      peratom_flag = 1;
      size_peratom_cols = 8;
      peratom_freq = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix wall/gran command");
  }

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all(FLERR,"Cannot use wall in periodic dimension");

  if (wiggle && wshear)
    error->all(FLERR,"Cannot wiggle and shear fix wall/gran");
  if (wiggle && wallstyle == ZCYLINDER && axis != 2)
    error->all(FLERR,"Invalid wiggle direction for fix wall/gran");
  if (wshear && wallstyle == XPLANE && axis == 0)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == YPLANE && axis == 1)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == ZPLANE && axis == 2)
    error->all(FLERR,"Invalid shear direction for fix wall/gran");
  if ((wiggle || wshear) && wallstyle == REGION)
    error->all(FLERR,"Cannot wiggle or shear with fix wall/gran/region");

  // setup oscillations

  if (wiggle) omega = 2.0*MY_PI / period;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  history_one = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  nmax = 0;
  mass_rigid = NULL;

  // initialize history as if particle is not touching region
  // history_one will be NULL for wallstyle = REGION

  if (use_history && history_one) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < size_history; j++)
        history_one[i][j] = 0.0;
  }

  if (peratom_flag) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < size_peratom_cols; m++)
        array_atom[i][m] = 0.0;
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete local storage

  delete [] idregion;
  memory->destroy(history_one);
  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
  int i;

  dt = update->dt;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = NULL;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];

  tangential_history_index = 0;
  if (roll_history) {
    if (tangential_history) roll_history_index = 3;
    else roll_history_index = 0;
  }
  if (twist_history) {
    if (tangential_history) {
      if (roll_history) twist_history_index = 6;
      else twist_history_index = 3;
    }
    else{
      if (roll_history) twist_history_index = 3;
      else twist_history_index = 0;
    }
  }
  if (normal_model == JKR) {
    tangential_history_index += 1;
    roll_history_index += 1;
    twist_history_index += 1;
  }
  if (tangential_model == TANGENTIAL_MINDLIN_RESCALE) {
    roll_history_index += 1;
    twist_history_index += 1;
  }

  if (damping_model == TSUJI) {
    double cor = normal_coeffs[1];
    normal_coeffs[1] = 1.2728-4.2783*cor+11.087*pow(cor,2)-22.348*pow(cor,3)+
        27.467*pow(cor,4)-18.022*pow(cor,5)+
        4.8218*pow(cor,6);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force(int /*vflag*/)
{
  int i,j;
  double dx,dy,dz,del1,del2,delxy,delr,rsq,rwall,meff;
  double vwall[3];

  // do not update history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"wall/gran:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    }
  }

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set history if pair potential stores history

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  rwall = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > radius[i]) {
          dz = cylradius;
          rwall = 0.0;
        } else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
          // rwall = -2r_c if inside cylinder, 2r_c outside
          rwall = (delxy < cylradius) ? -2*cylradius : 2*cylradius;
          if (wshear && axis != 2) {
            vwall[0] += vshear * x[i][1]/delxy;
            vwall[1] += -vshear * x[i][0]/delxy;
            vwall[2] = 0.0;
          }
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;

      double rad;
      if (pairstyle == GRANULAR && normal_model == JKR) {
        rad = radius[i] + pulloff_distance(radius[i]);
      }
      else
        rad = radius[i];

      if (rsq > rad*rad) {
        if (use_history)
          for (j = 0; j < size_history; j++)
            history_one[i][j] = 0.0;
      }
      else {
        if (pairstyle == GRANULAR && normal_model == JKR && use_history) {
          if ((history_one[i][0] == 0) && (rsq > radius[i]*radius[i])) {
            // Particles have not contacted yet,
            // and are outside of contact distance
            for (j = 0; j < size_history; j++)
              history_one[i][j] = 0.0;
            continue;
          }
        }

        // meff = effective mass of sphere
        // if I is part of rigid body, use body mass

        meff = rmass[i];
        if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

        // store contact info
        if (peratom_flag) {
          array_atom[i][0] = (double)atom->tag[i];
          array_atom[i][4] = x[i][0] - dx;
          array_atom[i][5] = x[i][1] - dy;
          array_atom[i][6] = x[i][2] - dz;
          array_atom[i][7] = radius[i];
        }

        // invoke sphere/wall interaction
        double *contact;
        if (peratom_flag)
          contact = array_atom[i];
        else
          contact = NULL;

        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,vwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff, contact);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,vwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,history_one[i],
              contact);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,history_one[i],
              contact);
        else if (pairstyle == GRANULAR)
          granular(rsq,dx,dy,dz,vwall,rwall,v[i],f[i],
              omega[i],torque[i],radius[i],meff,history_one[i],
              contact);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke(double rsq, double dx, double dy, double dz,
    double *vwall, double *v,
    double *f, double *omega, double *torque,
    double radius, double meff, double* contact)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = meff*gammat*vrel;
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // tangential force due to tangential velocity damping

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  if (peratom_flag) {
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hooke_history(double rsq, double dx, double dy, double dz,
    double *vwall, double *v,
    double *f, double *omega, double *torque,
    double radius, double meff, double *history,
    double *contact)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (history_update) {
    history[0] += vtr1*dt;
    history[1] += vtr2*dt;
    history[2] += vtr3*dt;
  }
  shrmag = sqrt(history[0]*history[0] + history[1]*history[1] + 
                history[2]*history[2]);

  // rotate shear displacements

  rsht = history[0]*dx + history[1]*dy + history[2]*dz;
  rsht = rsht*rsqinv;
  if (history_update) {
    history[0] -= rsht*dx;
    history[1] -= rsht*dy;
    history[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = - (kt*history[0] + meff*gammat*vtr1);
  fs2 = - (kt*history[1] + meff*gammat*vtr2);
  fs3 = - (kt*history[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      history[0] = (fn/fs) * (history[0] + meff*gammat*vtr1/kt) -
          meff*gammat*vtr1/kt;
      history[1] = (fn/fs) * (history[1] + meff*gammat*vtr2/kt) -
          meff*gammat*vtr2/kt;
      history[2] = (fn/fs) * (history[2] + meff*gammat*vtr3/kt) -
          meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  if (peratom_flag) {
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::hertz_history(double rsq, double dx, double dy, double dz,
    double *vwall, double rwall, double *v,
    double *f, double *omega, double *torque,
    double radius, double meff, double *history,
    double *contact)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,polyhertz,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr / rsq;
  vn2 = dy*vnnr / rsq;
  vn3 = dz*vnnr / rsq;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hertzian contact + normal velocity damping
  // rwall = 0 is flat wall case
  // rwall positive or negative is curved wall
  //   will break (as it should) if rwall is negative and
  //   its absolute value < radius of particle

  damp = meff*gamman*vnnr*rsqinv;
  ccel = kn*(radius-r)*rinv - damp;
  if (rwall == 0.0) polyhertz = sqrt((radius-r)*radius);
  else polyhertz = sqrt((radius-r)*radius*rwall/(rwall+radius));
  ccel *= polyhertz;

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  if (history_update) {
    history[0] += vtr1*dt;
    history[1] += vtr2*dt;
    history[2] += vtr3*dt;
  }
  shrmag = sqrt(history[0]*history[0] + history[1]*history[1] + 
                history[2]*history[2]);

  // rotate history displacements

  rsht = history[0]*dx + history[1]*dy + history[2]*dz;
  rsht = rsht*rsqinv;
  if (history_update) {
    history[0] -= rsht*dx;
    history[1] -= rsht*dy;
    history[2] -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt*history[0] + meff*gammat*vtr1);
  fs2 = -polyhertz * (kt*history[1] + meff*gammat*vtr2);
  fs3 = -polyhertz * (kt*history[2] + meff*gammat*vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      history[0] = (fn/fs) * (history[0] + meff*gammat*vtr1/kt) -
          meff*gammat*vtr1/kt;
      history[1] = (fn/fs) * (history[1] + meff*gammat*vtr2/kt) -
          meff*gammat*vtr2/kt;
      history[2] = (fn/fs) * (history[2] + meff*gammat*vtr3/kt) -
          meff*gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  if (peratom_flag) {
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::granular(double rsq, double dx, double dy, double dz,
                           double *vwall, double rwall, double *v,
                           double *f, double *omega, double *torque,
                           double radius, double meff, double *history,
                           double *contact)
{
  double fx,fy,fz,nx,ny,nz;
  double radsum,r,rinv;
  double Reff, delta, dR, dR2;

  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;

  double knfac, damp_normal, damp_normal_prefactor;
  double k_tangential, damp_tangential;
  double Fne, Ft, Fdamp, Fntot, Fncrit, Fscrit, Frcrit;
  double fs, fs1, fs2, fs3;

  double tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3;

  // for JKR
  double R2, coh, F_pulloff, a, a2, E;
  double t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  // rolling
  double k_roll, damp_roll;
  double torroll1, torroll2, torroll3;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;

  // twisting
  double k_twist, damp_twist, mu_twist;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double tortwist1, tortwist2, tortwist3;

  double shrmag,rsht;

  r = sqrt(rsq);
  radsum = rwall + radius;

  E = normal_coeffs[0];

  radsum = radius + rwall;
  if (rwall == 0) Reff = radius;
  else Reff = radius*rwall/(radius+rwall);

  rinv = 1.0/r;

  nx = dx*rinv;
  ny = dy*rinv;
  nz = dz*rinv;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*nx + vr2*ny + vr3*nz; //v_R . n
  vn1 = nx*vnnr;
  vn2 = ny*vnnr;
  vn3 = nz*vnnr;

  delta = radsum - r;
  dR = delta*Reff;
  if (normal_model == JKR) {
    history[0] = 1.0;
    E *= THREEQUARTERS;
    R2=Reff*Reff;
    coh = normal_coeffs[3];
    dR2 = dR*dR;
    t0 = coh*coh*R2*R2*E;
    t1 = PI27SQ*t0;
    t2 = 8*dR*dR2*E*E*E;
    t3 = 4*dR2*E;
    sqrt1 = MAX(0, t0*(t1+2*t2)); // in case sqrt(0) < 0 due to precision issues
    t4 = cbrt(t1+t2+THREEROOT3*M_PI*sqrt(sqrt1));
    t5 = t3/t4 + t4/E;
    sqrt2 = MAX(0, 2*dR + t5);
    t6 = sqrt(sqrt2);
    sqrt3 = MAX(0, 4*dR - t5 + SIXROOT6*coh*M_PI*R2/(E*t6));
    a = INVROOT6*(t6 + sqrt(sqrt3));
    a2 = a*a;
    knfac = normal_coeffs[0]*a;
    Fne = knfac*a2/Reff - TWOPI*a2*sqrt(4*coh*E/(M_PI*a));
  }
  else{
    knfac = E; //Hooke
    a = sqrt(dR);
    if (normal_model != HOOKE) {
      Fne *= a;
      knfac *= a;
    }
    Fne = knfac*delta;
    if (normal_model == DMT)
      Fne -= 4*MY_PI*normal_coeffs[3]*Reff;
  }

  if (damping_model == VELOCITY) {
    damp_normal = 1;
  }
  else if (damping_model == VISCOELASTIC) {
    damp_normal = a*meff;
  }
  else if (damping_model == TSUJI) {
    damp_normal = sqrt(meff*knfac);
  }

  damp_normal_prefactor = normal_coeffs[1]*damp_normal;
  Fdamp = -damp_normal_prefactor*vnnr;

  Fntot = Fne + Fdamp;

  //****************************************
  // tangential force, including history effects
  //****************************************

  // tangential component
  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity
  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // relative tangential velocities
  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  if (normal_model == JKR) {
    F_pulloff = 3*M_PI*coh*Reff;
    Fncrit = fabs(Fne + 2*F_pulloff);
  }
  else if (normal_model == DMT) {
    F_pulloff = 4*M_PI*coh*Reff;
    Fncrit = fabs(Fne + 2*F_pulloff);
  }
  else{
    Fncrit = fabs(Fntot);
  }

  //------------------------------
  // tangential forces
  //------------------------------

  k_tangential = tangential_coeffs[0];
  damp_tangential = tangential_coeffs[1]*damp_normal_prefactor;

  int thist0 = tangential_history_index;
  int thist1 = thist0 + 1;
  int thist2 = thist1 + 1;

  if (tangential_history) {
    if (tangential_model == TANGENTIAL_MINDLIN) {
      k_tangential *= a;
    }
    else if (tangential_model == TANGENTIAL_MINDLIN_RESCALE) {
      k_tangential *= a;
      if (a < history[3]) { //On unloading, rescale the shear displacements
        double factor = a/history[thist2+1];
        history[thist0] *= factor;
        history[thist1] *= factor;
        history[thist2] *= factor;
      }
    }
    shrmag = sqrt(history[thist0]*history[thist0] + 
                  history[thist1]*history[thist1] +
                  history[thist2]*history[thist2]);

    // rotate and update displacements.
    // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
    if (history_update) {
      rsht = history[thist0]*nx + history[thist1]*ny + history[thist2]*nz;
      if (fabs(rsht) < EPSILON) rsht = 0;
      if (rsht > 0) {
        // if rhst == shrmag, contacting pair has rotated 90 deg in one step,
        // in which case you deserve a crash!
        scalefac = shrmag/(shrmag - rsht); 
        history[thist0] -= rsht*nx;
        history[thist1] -= rsht*ny;
        history[thist2] -= rsht*nz;
        // also rescale to preserve magnitude
        history[thist0] *= scalefac;
        history[thist1] *= scalefac;
        history[thist2] *= scalefac;
      }
      // update history
      history[thist0] += vtr1*dt;
      history[thist1] += vtr2*dt;
      history[thist2] += vtr3*dt;
    }

    // tangential forces = history + tangential velocity damping
    fs1 = -k_tangential*history[thist0] - damp_tangential*vtr1;
    fs2 = -k_tangential*history[thist1] - damp_tangential*vtr2;
    fs3 = -k_tangential*history[thist2] - damp_tangential*vtr3;

    // rescale frictional displacements and forces if needed
    Fscrit = tangential_coeffs[2] * Fncrit;
    fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    if (fs > Fscrit) {
      if (shrmag != 0.0) {
        history[thist0] = -1.0/k_tangential*(Fscrit*fs1/fs + 
                                             damp_tangential*vtr1);
        history[thist1] = -1.0/k_tangential*(Fscrit*fs2/fs + 
                                             damp_tangential*vtr2);
        history[thist2] = -1.0/k_tangential*(Fscrit*fs3/fs + 
                                             damp_tangential*vtr3);
        fs1 *= Fscrit/fs;
        fs2 *= Fscrit/fs;
        fs3 *= Fscrit/fs;
      } else fs1 = fs2 = fs3 = 0.0;
    }
  } else { // classic pair gran/hooke (no history)
    fs = meff*damp_tangential*vrel;
    if (vrel != 0.0) Ft = MIN(Fne,fs) / vrel;
    else Ft = 0.0;
    fs1 = -Ft*vtr1;
    fs2 = -Ft*vtr2;
    fs3 = -Ft*vtr3;
  }

  //****************************************
  // rolling resistance
  //****************************************

  if (roll_model != ROLL_NONE) {
    relrot1 = omega[0];
    relrot2 = omega[1];
    relrot3 = omega[2];

    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // This is different from the Marshall papers,
    // which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl1 = Reff*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
    vrl2 = Reff*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
    vrl3 = Reff*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;

    int rhist0 = roll_history_index;
    int rhist1 = rhist0 + 1;
    int rhist2 = rhist1 + 1;

    // rolling displacement
    rollmag = sqrt(history[rhist0]*history[rhist0] +
                   history[rhist1]*history[rhist1] +
                   history[rhist2]*history[rhist2]);

    rolldotn = history[rhist0]*nx + history[rhist1]*ny + history[rhist2]*nz;

    if (history_update) {
      if (fabs(rolldotn) < EPSILON) rolldotn = 0;
      if (rolldotn > 0) { // rotate into tangential plane
        scalefac = rollmag/(rollmag - rolldotn);
        history[rhist0] -= rolldotn*nx;
        history[rhist1] -= rolldotn*ny;
        history[rhist2] -= rolldotn*nz;
        // also rescale to preserve magnitude
        history[rhist0] *= scalefac;
        history[rhist1] *= scalefac;
        history[rhist2] *= scalefac;
      }
      history[rhist0] += vrl1*dt;
      history[rhist1] += vrl2*dt;
      history[rhist2] += vrl3*dt;
    }

    k_roll = roll_coeffs[0];
    damp_roll = roll_coeffs[1];
    fr1 = -k_roll*history[rhist0] - damp_roll*vrl1;
    fr2 = -k_roll*history[rhist1] - damp_roll*vrl2;
    fr3 = -k_roll*history[rhist2] - damp_roll*vrl3;

    // rescale frictional displacements and forces if needed
    Frcrit = roll_coeffs[2] * Fncrit;

    fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
    if (fr > Frcrit) {
      if (rollmag != 0.0) {
        history[rhist0] = -1.0/k_roll*(Frcrit*fr1/fr + damp_roll*vrl1);
        history[rhist1] = -1.0/k_roll*(Frcrit*fr2/fr + damp_roll*vrl2);
        history[rhist2] = -1.0/k_roll*(Frcrit*fr3/fr + damp_roll*vrl3);
        fr1 *= Frcrit/fr;
        fr2 *= Frcrit/fr;
        fr3 *= Frcrit/fr;
      } else fr1 = fr2 = fr3 = 0.0;
    }
  }

  //****************************************
  // twisting torque, including history effects
  //****************************************

  if (twist_model != TWIST_NONE) {
    magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
    if (twist_model == TWIST_MARSHALL) {
      k_twist = 0.5*k_tangential*a*a;; // eq 32 of Marshall paper
      damp_twist = 0.5*damp_tangential*a*a;
      mu_twist = TWOTHIRDS*a*tangential_coeffs[2];
    }
    else{
      k_twist = twist_coeffs[0];
      damp_twist = twist_coeffs[1];
      mu_twist = twist_coeffs[2];
    }
    if (history_update) {
      history[twist_history_index] += magtwist*dt;
    }
    // M_t torque (eq 30)
    magtortwist = -k_twist*history[twist_history_index] - damp_twist*magtwist;
    signtwist = (magtwist > 0) - (magtwist < 0);
    Mtcrit = mu_twist*Fncrit; // critical torque (eq 44)
    if (fabs(magtortwist) > Mtcrit) {
      history[twist_history_index] = 1.0/k_twist*(Mtcrit*signtwist - 
                                                  damp_twist*magtwist);
      magtortwist = -Mtcrit * signtwist; // eq 34
    }
  }

  // apply forces & torques

  fx = nx*Fntot + fs1;
  fy = ny*Fntot + fs2;
  fz = nz*Fntot + fs3;

  if (peratom_flag) {
    contact[1] = fx;
    contact[2] = fy;
    contact[3] = fz;
  }

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = ny*fs3 - nz*fs2;
  tor2 = nz*fs1 - nx*fs3;
  tor3 = nx*fs2 - ny*fs1;

  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;

  if (twist_model != TWIST_NONE) {
    tortwist1 = magtortwist * nx;
    tortwist2 = magtortwist * ny;
    tortwist3 = magtortwist * nz;

    torque[0] += tortwist1;
    torque[1] += tortwist2;
    torque[2] += tortwist3;
  }

  if (roll_model != ROLL_NONE) {
    torroll1 = Reff*(ny*fr3 - nz*fr2); //n cross fr
    torroll2 = Reff*(nz*fr1 - nx*fr3);
    torroll3 = Reff*(nx*fr2 - ny*fr1);

    torque[0] += torroll1;
    torque[1] += torroll2;
    torque[2] += torroll3;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  if (use_history) bytes += nmax*size_history * sizeof(double);  // shear history
  if (fix_rigid) bytes += nmax * sizeof(int);                    // mass_rigid
  // store contacts
  if (peratom_flag) bytes += nmax*size_peratom_cols*sizeof(double); 
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  if (use_history) memory->grow(history_one,nmax,size_history,
                                "fix_wall_gran:history_one");
  if (peratom_flag) {
    memory->grow(array_atom,nmax,size_peratom_cols,"fix_wall_gran:array_atom");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::copy_arrays(int i, int j, int /*delflag*/)
{
  if (use_history)
    for (int m = 0; m < size_history; m++)
      history_one[j][m] = history_one[i][m];
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[j][m] = array_atom[i][m];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGran::set_arrays(int i)
{
  if (use_history)
    for (int m = 0; m < size_history; m++)
      history_one[i][m] = 0;
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[i][m] = 0;
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGran::pack_exchange(int i, double *buf)
{
  int n = 0;
  if (use_history) {
    for (int m = 0; m < size_history; m++)
      buf[n++] = history_one[i][m];
  }
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++)
      buf[n++] = array_atom[i][m];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGran::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  if (use_history) {
    for (int m = 0; m < size_history; m++)
      history_one[nlocal][m] = buf[n++];
  }
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++)
      array_atom[nlocal][m] = buf[n++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGran::pack_restart(int i, double *buf)
{
  if (!use_history) return 0;

  int n = 0;
  buf[n++] = size_history + 1;
  for (int m = 0; m < size_history; m++)
    buf[n++] = history_one[i][m];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGran::unpack_restart(int nlocal, int nth)
{
  if (!use_history) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int i = 0; i < size_history; i++)
    history_one[nlocal][i] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGran::maxsize_restart()
{
  if (!use_history) return 0;
  return 1 + size_history;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGran::size_restart(int /*nlocal*/)
{
  if (!use_history) return 0;
  return 1 + size_history;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::reset_dt()
{
  dt = update->dt;
}

double FixWallGran::pulloff_distance(double radius)
{
  double coh, E, a, dist;
  coh = normal_coeffs[3];
  E = normal_coeffs[0]*THREEQUARTERS;
  a = cbrt(9*M_PI*coh*radius/(4*E));
  dist = a*a/radius - 2*sqrt(M_PI*coh*a/E);
  return dist;
}

