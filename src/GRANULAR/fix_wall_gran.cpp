// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL),
                         Dan Bolintineanu (SNL), Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_wall_gran.h"

#include "atom.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

static constexpr double BIG = 1.0e20;

// XYZ PLANE need to be 0,1,2

enum {NOSTYLE=-1,XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum {NONE,CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(nullptr), tstr(nullptr), history_one(nullptr),
  fix_rigid(nullptr), mass_rigid(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR,"fix wall/gran", error);

  if (!atom->omega_flag) error->all(FLERR,"Fix wall/gran requires atom attribute omega");
  if (!atom->radius_flag) error->all(FLERR,"Fix wall/gran requires atom attribute radius");

  create_attribute = 1;

  // set interaction style
  // disable bonded/history option for now
  model = new GranularModel(lmp);
  model->contact_type = WALL;

  heat_flag = 0;
  int classic_flag = 1;
  if (strcmp(arg[3],"granular") == 0)  classic_flag = 0;

  // wall/particle coefficients

  int iarg;
  if (classic_flag) {
    iarg = model->define_classic_model(arg, 3, narg);

    if (iarg < narg) {
      if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      }
    }

  } else {
    iarg = 4;
    iarg = model->add_sub_model(arg, iarg, narg, NORMAL);

    while (iarg < narg) {
      if (strcmp(arg[iarg], "damping") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
      } else if (strcmp(arg[iarg], "tangential") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
      } else if (strcmp(arg[iarg], "rolling") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, ROLLING);
      } else if (strcmp(arg[iarg], "twisting") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TWISTING);
      } else if (strcmp(arg[iarg], "heat") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, HEAT);
        heat_flag = 1;
      } else if (strcmp(arg[iarg], "xplane") == 0 ||
          strcmp(arg[iarg], "yplane") == 0 ||
          strcmp(arg[iarg], "zplane") == 0 ||
          strcmp(arg[iarg], "zcylinder") == 0 ||
          strcmp(arg[iarg], "region") == 0) {
        break;
      } else if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      } else {
        error->all(FLERR, "Unknown fix wall/gran keyword {}", arg[iarg]);
      }
    }
  }

  // Define default damping sub model if unspecified, takes no args
  if (!model->damping_model)
    model->construct_sub_model("viscoelastic", DAMPING);
  model->init();

  size_history = model->size_history;
  if (model->beyond_contact) size_history += 1; //Need to track if particle is touching
  if (size_history == 0) use_history = restart_peratom = 0;
  else use_history = restart_peratom = 1;

  // wallstyle args

  if (iarg >= narg) error->all(FLERR, "Illegal fix wall/gran command");

  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    iarg += 2;
  } else if (strcmp(arg[iarg],"region") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/gran command");
    wallstyle = REGION;
    idregion = utils::strdup(arg[iarg+1]);
    iarg += 2;
  } else wallstyle = NOSTYLE;

  // optional args

  wiggle = 0;
  wshear = 0;
  peratom_flag = 0;
  int Twall_defined = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      amplitude = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      period = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      wiggle = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"shear") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/gran command");
      vshear = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      wshear = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"contacts") == 0) {
      peratom_flag = 1;
      size_peratom_cols = 8;
      peratom_freq = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/gran command");
      if (utils::strmatch(arg[iarg+1], "^v_")) {
        tstr = utils::strdup(arg[iarg+1] + 2);
      } else {
        Twall = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      Twall_defined = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/gran command");
  }

  if (heat_flag != Twall_defined)
    error->all(FLERR, "Must define wall temperature with heat model");

  if (wallstyle == NOSTYLE)
    error->all(FLERR,"No wall style defined");
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

  FixWallGran::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  nmax = 0;

  // initialize history as if particle is not touching region
  // history_one will be a null pointer for wallstyle = REGION

  if (use_history && history_one) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < size_history; j++)
        history_one[i][j] = 0.0;
  }

  if (peratom_flag) {
    clear_stored_contacts();
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
  if (copymode) return;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  // delete local storage

  delete model;
  delete[] tstr;
  delete[] idregion;
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
  model->dt = dt;

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // check for compatible heat conduction atom style

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR,"Heat conduction in fix wall/gran requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,"Heat conduction in fix wall/gran requires atom style with heatflow property");
  }

  // check for FixRigid so can extract rigid body masses

  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];

  // Define history indices

  int next_index = 0;
  if (model->beyond_contact) next_index = 1;

  for (i = 0; i < NSUBMODELS; i++) {
    model->sub_models[i]->history_index = next_index;
    next_index += model->sub_models[i]->size_history;
  }

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR, "Variable {} for fix wall/gran does not exist", tstr);
    if (! input->variable->equalstyle(tvar))
      error->all(FLERR, "Variable {} for fix wall/gran must be an equal style variable", tstr);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force(int /*vflag*/)
{
  int i,j;
  double dx,dy,dz,del1,del2,delxy,delr,rwall,meff;
  double *forces, *torquesi;
  double vwall[3];
  double w0[3] = {0.0};
  bool touchflag = false;

  // do not update history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;
  model->history_update = history_update;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal",tmp);
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
      wlo = lo + amplitude - amplitude * cos(arg);
      whi = hi + amplitude - amplitude * cos(arg);
    }
    vwall[axis] = amplitude * omega * sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to a null pointer, it's skipped since lo/hi are infinity
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
  double *temperature, *heatflow;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  rwall = 0.0;

  if (peratom_flag) {
    clear_stored_contacts();
  }

  // Define constant wall properties (atom j)
  model->radj = 0.0;
  model->vj = vwall;
  model->omegaj = w0;
  if (heat_flag) {
    temperature = atom->temperature;
    heatflow = atom->heatflow;
    if (tstr)
      Twall = input->variable->compute_equal(tvar);
    model->Tj = Twall;
  }

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

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
      delxy = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
      delr = cylradius - delxy;
      if (delr > radius[i]) {
        dz = cylradius;
        rwall = 0.0;
      } else {
        dx = -delr / delxy * x[i][0];
        dy = -delr / delxy * x[i][1];
        // rwall = -2r_c if inside cylinder, 2r_c outside
        rwall = (delxy < cylradius) ? -2 * cylradius : 2 * cylradius;
        if (wshear && axis != 2) {
          vwall[0] += vshear * x[i][1] / delxy;
          vwall[1] += -vshear * x[i][0] / delxy;
          vwall[2] = 0.0;
        }
      }
    }

    // Reset model and copy initial geometric data
    model->dx[0] = dx;
    model->dx[1] = dy;
    model->dx[2] = dz;
    model->radi = radius[i];
    model->radj = rwall;
    if (model->beyond_contact) model->touch = history_one[i][0];

    touchflag = model->check_contact();

    if (!touchflag) {
      if (use_history)
        for (j = 0; j < size_history; j++)
          history_one[i][j] = 0.0;
      continue;
    }

    if (model->beyond_contact)
      history_one[i][0] = 1;

    // meff = effective mass of sphere
    // if I is part of rigid body, use body mass

    meff = rmass[i];
    if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

    // Copy additional information and prepare force calculations
    model->meff = meff;
    model->vi = v[i];
    model->omegai = omega[i];
    if (use_history) model->history = history_one[i];
    if (heat_flag) model->Ti = temperature[i];

    model->calculate_forces();

    forces = model->forces;
    torquesi = model->torquesi;

    // apply forces & torques
    add3(f[i], forces, f[i]);

    add3(torque[i], torquesi, torque[i]);
    if (heat_flag) heatflow[i] += model->dq;

    // store contact info
    if (peratom_flag) {
      array_atom[i][0] = 1.0;
      array_atom[i][1] = forces[0];
      array_atom[i][2] = forces[1];
      array_atom[i][3] = forces[2];
      array_atom[i][4] = x[i][0] - dx;
      array_atom[i][5] = x[i][1] - dy;
      array_atom[i][6] = x[i][2] - dz;
      array_atom[i][7] = radius[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::clear_stored_contacts()
{
  const int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    for (int m = 0; m < size_peratom_cols; m++) {
      array_atom[i][m] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  if (use_history) bytes += (double)nmax * size_history * sizeof(double);  // shear history
  if (fix_rigid) bytes += (double)nmax * sizeof(int);                    // mass_rigid
  // store contacts
  if (peratom_flag) bytes += (double)nmax * size_peratom_cols * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGran::grow_arrays(int nmax)
{
  if (use_history) memory->grow(history_one,nmax,size_history,"fix_wall_gran:history_one");
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
  // pack buf[0] this way because other fixes unpack it
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
  // unpack the Nth first values this way because other fixes pack them

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
  model->dt = dt;
}
