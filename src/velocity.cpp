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

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "velocity.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "compute_temp.h"
#include "random_park.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{CREATE,SET,SCALE,RAMP,ZERO};
enum{ALL,LOCAL,GEOM};
enum{NONE,CONSTANT,EQUAL,ATOM};

#define WARMUP 100
#define SMALL  0.001

/* ---------------------------------------------------------------------- */

Velocity::Velocity(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Velocity::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal velocity command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Velocity command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Velocity command with no atoms existing");

  // atom masses must all be set

  atom->check_mass();

  // identify group

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find velocity group ID");
  groupbit = group->bitmask[igroup];

  // identify style

  if (strcmp(arg[1],"create") == 0) style = CREATE;
  else if (strcmp(arg[1],"set") == 0) style = SET;
  else if (strcmp(arg[1],"scale") == 0) style = SCALE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"zero") == 0) style = ZERO;
  else error->all(FLERR,"Illegal velocity command");

  // set defaults

  temperature = NULL;
  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 1;
  rotation_flag = 0;
  bias_flag = 0;
  loop_flag = ALL;
  scale_flag = 1;
  rfix = -1;

  // read options from end of input line
  // change defaults as options specify

  if (style == CREATE) options(narg-4,&arg[4]);
  else if (style == SET) options(narg-5,&arg[5]);
  else if (style == SCALE) options(narg-3,&arg[3]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == ZERO) options(narg-3,&arg[3]);

  // initialize velocities based on style
  // create() invoked differently, so can be called externally

  if (style == CREATE) {
    double t_desired = force->numeric(FLERR,arg[2]);
    int seed = force->inumeric(FLERR,arg[3]);
    create(t_desired,seed);
  }
  else if (style == SET) set(narg-2,&arg[2]);
  else if (style == SCALE) scale(narg-2,&arg[2]);
  else if (style == RAMP) ramp(narg-2,&arg[2]);
  else if (style == ZERO) zero(narg-2,&arg[2]);
}

/* ----------------------------------------------------------------------
   initialization of defaults before calling velocity methods externaly
------------------------------------------------------------------------- */

void Velocity::init_external(const char *extgroup)
{
  igroup = group->find(extgroup);
  if (igroup == -1) error->all(FLERR,"Could not find velocity group ID");
  groupbit = group->bitmask[igroup];

  temperature = NULL;
  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 1;
  rotation_flag = 0;
  loop_flag = ALL;
  scale_flag = 1;
}

/* ---------------------------------------------------------------------- */

void Velocity::create(double t_desired, int seed)
{
  int i;
  double **vhold;

  if (seed <= 0) error->all(FLERR,"Illegal velocity create command");

  // if sum_flag set, store a copy of current velocities

  if (sum_flag) {
    double **v = atom->v;
    int nlocal = atom->nlocal;
    memory->create(vhold,nlocal,3,"velocity:vhold");
    for (i = 0; i < nlocal; i++) {
      vhold[i][0] = v[i][0];
      vhold[i][1] = v[i][1];
      vhold[i][2] = v[i][2];
    }
  }

  // if temperature = NULL or bias_flag set,
  // create a new ComputeTemp with the velocity group

  int tcreate_flag = 0;
  Compute *temperature_nobias = NULL;

  if (temperature == NULL || bias_flag) {
    char **arg = new char*[3];
    arg[0] = (char *) "velocity_temp";
    arg[1] = group->names[igroup];
    arg[2] = (char *) "temp";
    if (temperature == NULL) {
      temperature = new ComputeTemp(lmp,3,arg);
      tcreate_flag = 1;
    } else temperature_nobias = new ComputeTemp(lmp,3,arg);
    delete [] arg;
  }

  // initialize temperature computation(s)
  // warn if groups don't match

  if (igroup != temperature->igroup && comm->me == 0)
    error->warning(FLERR,"Mismatch between velocity and compute groups");
  temperature->init();
  temperature->setup();
  if (temperature_nobias) {
    temperature_nobias->init();
    temperature_nobias->setup();
  }

  // if bias_flag set, remove bias velocity from all atoms
  // for some temperature computes, must first calculate temp to do that

  if (bias_flag) {
    temperature->compute_scalar();
    temperature->remove_bias_all();
  }

  // create new velocities, in uniform or gaussian distribution
  // loop option determines looping style, ALL is default
  //   ALL = loop over all natoms, only set those I own via atom->map
  //    cannot do this if atom IDs do not span 1-Natoms (some were deleted)
  //    will produce same V, independent of P, if atoms were read-in
  //    will NOT produce same V, independent of P, if used create_atoms
  //   LOCAL = only loop over my atoms, adjust RNG to be proc-specific
  //    will never produce same V, independent of P
  //   GEOM = only loop over my atoms
  //    choose RNG for each atom based on its xyz coord (geometry)
  //      via random->reset()
  //    will always produce same V, independent of P
  // adjust by factor for atom mass
  // set xdim,ydim,zdim = 1/0 for whether to create velocity in those dims
  //   zdim = 0 for 2d
  //   any dims can be 0 if bias temperature compute turns them off
  //     currently only temp/partial does

  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  int m;
  double vx,vy,vz,factor;
  RanPark *random;

  if (loop_flag == ALL) {

    // create an atom map if one doesn't exist already

    int mapflag = 0;
    if (atom->map_style == 0) {
      mapflag = 1;
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }

    // error check

    if (atom->natoms > MAXSMALLINT)
      error->all(FLERR,"Too big a problem to use velocity create loop all");
    if (atom->tag_enable == 0)
      error->all(FLERR,
                 "Cannot use velocity create loop all unless atoms have IDs");
    if (atom->tag_consecutive() == 0)
      error->all(FLERR,
                 "Atom IDs must be consecutive for velocity create loop all");

    // loop over all atoms in system
    // generate RNGs for all atoms, only assign to ones I own
    // use either per-type mass or per-atom rmass

    random = new RanPark(lmp,seed);
    int natoms = static_cast<int> (atom->natoms);

    for (i = 1; i <= natoms; i++) {
      if (dist_flag == 0) {
        vx = random->uniform() - 0.5;
        vy = random->uniform() - 0.5;
        vz = random->uniform() - 0.5;
      } else {
        vx = random->gaussian();
        vy = random->gaussian();
        vz = random->gaussian();
      }
      m = atom->map(i);
      if (m >= 0 && m < nlocal) {
        if (mask[m] & groupbit) {
          if (rmass) factor = 1.0/sqrt(rmass[m]);
          else factor = 1.0/sqrt(mass[type[m]]);
          v[m][0] = vx * factor;
          v[m][1] = vy * factor;
          if (dim == 3) v[m][2] = vy * factor;
          else v[m][2] = 0.0;
        }
      }
    }

    // delete temporary atom map

    if (mapflag) {
      atom->map_delete();
      atom->map_style = 0;
    }

  } else if (loop_flag == LOCAL) {
    random = new RanPark(lmp,seed + comm->me);
    for (i = 0; i < WARMUP; i++) random->uniform();

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (dist_flag == 0) {
          vx = random->uniform() - 0.5;
          vy = random->uniform() - 0.5;
          vz = random->uniform() - 0.5;
        } else {
          vx = random->gaussian();
          vy = random->gaussian();
          vz = random->gaussian();
        }
        if (rmass) factor = 1.0/sqrt(rmass[i]);
        else factor = 1.0/sqrt(mass[type[i]]);
        v[i][0] = vx * factor;
        v[i][1] = vy * factor;
        if (dim == 3) v[i][2] = vz * factor;
        else v[i][2] = 0.0;
      }
    }

  } else if (loop_flag == GEOM) {
    random = new RanPark(lmp,1);
    double **x = atom->x;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        random->reset(seed,x[i]);
        if (dist_flag == 0) {
          vx = random->uniform() - 0.5;
          vy = random->uniform() - 0.5;
          vz = random->uniform() - 0.5;
        } else {
          vx = random->gaussian();
          vy = random->gaussian();
          vz = random->gaussian();
        }

        if (rmass) factor = 1.0/sqrt(rmass[i]);
        else factor = 1.0/sqrt(mass[type[i]]);
        v[i][0] = vx * factor;
        v[i][1] = vy * factor;
        if (dim == 3) v[i][2] = vz * factor;
        else v[i][2] = 0.0;
      }
    }
  }

  // apply momentum and rotation zeroing

  if (momentum_flag) zero_momentum();
  if (rotation_flag) zero_rotation();

  // scale temp to desired value
  // if bias flag is set, bias velocities have already been removed:
  //   no-bias compute calculates temp only for new thermal velocities

  double t;
  if (bias_flag == 0) t = temperature->compute_scalar();
  else t = temperature_nobias->compute_scalar();
  rescale(t,t_desired);

  // if bias_flag set, restore bias velocity to all atoms
  // reapply for bias = compute temp/partial to reset v dims to 0.0

  if (bias_flag) {
    temperature->reapply_bias_all();
    temperature->restore_bias_all();
  }

  // if sum_flag set, add back in previous velocities

  if (sum_flag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] += vhold[i][0];
        v[i][1] += vhold[i][1];
        v[i][2] += vhold[i][2];
      }
    }
    memory->destroy(vhold);
  }

  // free local memory
  // if temperature compute was created, delete it

  delete random;
  if (tcreate_flag) delete temperature;
  if (temperature_nobias) delete temperature_nobias;
}

/* ---------------------------------------------------------------------- */

void Velocity::set(int narg, char **arg)
{
  int xstyle,ystyle,zstyle,varflag;
  double vx,vy,vz;
  char *xstr,*ystr,*zstr;
  int xvar,yvar,zvar;

  // parse 3 args

  xstyle = ystyle = zstyle = CONSTANT;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[0],"v_") == arg[0]) {
    int n = strlen(&arg[0][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[0][2]);
  } else if (strcmp(arg[0],"NULL") == 0) xstyle = NONE;
  else vx = force->numeric(FLERR,arg[0]);

  if (strstr(arg[1],"v_") == arg[1]) {
    int n = strlen(&arg[1][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[1][2]);
  } else if (strcmp(arg[1],"NULL") == 0) ystyle = NONE;
  else vy = force->numeric(FLERR,arg[1]);

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[2][2]);
  } else if (strcmp(arg[2],"NULL") == 0) zstyle = NONE;
  else vz = force->numeric(FLERR,arg[2]);

  // set and apply scale factors

  xscale = yscale = zscale = 1.0;

  if (xstyle && !xstr) {
    if (scale_flag) xscale = domain->lattice->xlattice;
    vx *= xscale;
  }
  if (ystyle && !ystr) {
    if (scale_flag) yscale = domain->lattice->ylattice;
    vy *= yscale;
  }
  if (zstyle && !zstr) {
    if (scale_flag) zscale = domain->lattice->zlattice;
    vz *= zscale;
  }

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  // error check for 2d models

  if (domain->dimension == 2) {
    if (zstyle == CONSTANT && vz != 0.0)
      error->all(FLERR,"Cannot set non-zero z velocity for 2d simulation");
    if (zstyle == EQUAL || zstyle == ATOM)
      error->all(FLERR,"Cannot set variable z velocity for 2d simulation");
  }

  // allocate vfield array if necessary

  double **vfield = NULL;
  if (varflag == ATOM) memory->create(vfield,atom->nlocal,3,"velocity:vfield");

  // set velocities via constants

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle) v[i][0] = vx;
          if (ystyle) v[i][1] = vy;
          if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle) v[i][0] += vx;
          if (ystyle) v[i][1] += vy;
          if (zstyle) v[i][2] += vz;
        }
      }
    }

  // set velocities via variables

  } else {
    if (xstyle == EQUAL) vx = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM) {
      if (vfield) input->variable->compute_atom(xvar,igroup,&vfield[0][0],3,0);
      else input->variable->compute_atom(xvar,igroup,NULL,3,0);
    }
    if (ystyle == EQUAL) vy = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM) {
      if (vfield) input->variable->compute_atom(yvar,igroup,&vfield[0][1],3,0);
      else input->variable->compute_atom(yvar,igroup,NULL,3,0);
    }
    if (zstyle == EQUAL) vz = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM) {
      if (vfield) input->variable->compute_atom(zvar,igroup,&vfield[0][2],3,0);
      else input->variable->compute_atom(zvar,igroup,NULL,3,0);
    }

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle == ATOM) v[i][0] = vfield[i][0];
          else if (xstyle) v[i][0] = vx;
          if (ystyle == ATOM) v[i][1] = vfield[i][1];
          else if (ystyle) v[i][1] = vy;
          if (zstyle == ATOM) v[i][2] = vfield[i][2];
          else if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle == ATOM) v[i][0] += vfield[i][0];
          else if (xstyle) v[i][0] += vx;
          if (ystyle == ATOM) v[i][1] += vfield[i][1];
          else if (ystyle) v[i][1] += vy;
          if (zstyle == ATOM) v[i][2] += vfield[i][2];
          else if (zstyle) v[i][2] += vz;
        }
      }
  }

  // clean up

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy(vfield);
}

/* ----------------------------------------------------------------------
   rescale velocities of a group after computing its temperature
------------------------------------------------------------------------- */

void Velocity::scale(int narg, char **arg)
{
  double t_desired = force->numeric(FLERR,arg[0]);

  // if temperature = NULL, create a new ComputeTemp with the velocity group

  int tflag = 0;
  if (temperature == NULL) {
    char **arg = new char*[3];
    arg[0] = (char *) "velocity_temp";
    arg[1] = group->names[igroup];
    arg[2] = (char *) "temp";
    temperature = new ComputeTemp(lmp,3,arg);
    tflag = 1;
    delete [] arg;
  }

  // initialize temperature computation
  // warn if groups don't match

  if (igroup != temperature->igroup && comm->me == 0)
    error->warning(FLERR,"Mismatch between velocity and compute groups");
  temperature->init();
  temperature->setup();

  // scale temp to desired value
  // if bias flag is set:
  //   temperature calculation will be done accounting for bias
  //   remove/restore bias velocities before/after rescale

  if (bias_flag == 0) {
    double t = temperature->compute_scalar();
    rescale(t,t_desired);
  } else {
    double t = temperature->compute_scalar();
    temperature->remove_bias_all();
    rescale(t,t_desired);
    temperature->restore_bias_all();
  }

  // if temperature was created, delete it

  if (tflag) delete temperature;
}

/* ----------------------------------------------------------------------
   apply a ramped set of velocities
------------------------------------------------------------------------- */

void Velocity::ramp(int narg, char **arg)
{
  // set scale factors

  if (scale_flag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // parse args

  int v_dim;
  if (strcmp(arg[0],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[0],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[0],"vz") == 0) v_dim = 2;
  else error->all(FLERR,"Illegal velocity command");

  if (v_dim == 2 && domain->dimension == 2)
    error->all(FLERR,"Velocity ramp in z for a 2d problem");

  double v_lo,v_hi;
  if (v_dim == 0) {
    v_lo = xscale*force->numeric(FLERR,arg[1]);
    v_hi = xscale*force->numeric(FLERR,arg[2]);
  } else if (v_dim == 1) {
    v_lo = yscale*force->numeric(FLERR,arg[1]);
    v_hi = yscale*force->numeric(FLERR,arg[2]);
  } else if (v_dim == 2) {
    v_lo = zscale*force->numeric(FLERR,arg[1]);
    v_hi = zscale*force->numeric(FLERR,arg[2]);
  }

  int coord_dim;
  if (strcmp(arg[3],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[3],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[3],"z") == 0) coord_dim = 2;
  else error->all(FLERR,"Illegal velocity command");

  double coord_lo,coord_hi;
  if (coord_dim == 0) {
    coord_lo = xscale*force->numeric(FLERR,arg[4]);
    coord_hi = xscale*force->numeric(FLERR,arg[5]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*force->numeric(FLERR,arg[4]);
    coord_hi = yscale*force->numeric(FLERR,arg[5]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*force->numeric(FLERR,arg[4]);
    coord_hi = zscale*force->numeric(FLERR,arg[5]);
  }

  // vramp = ramped velocity component for v_dim
  // add or set based on sum_flag

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fraction,vramp;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      if (sum_flag) v[i][v_dim] += vramp;
      else v[i][v_dim] = vramp;
    }
}

/* ----------------------------------------------------------------------
   zero linear or angular momentum of a group
   if using rigid/small requires init of entire system since
      its methods perform forward/reverse comm,
      comm::init needs neighbor::init needs pair::init needs kspace::init, etc
      also requires setup_pre_neighbor call to setup bodies
------------------------------------------------------------------------- */

void Velocity::zero(int narg, char **arg)
{
  if (strcmp(arg[0],"linear") == 0) {
    if (rfix < 0) zero_momentum();
    else {
      if (strcmp(modify->fix[rfix]->style,"rigid/small") == 0) {
        lmp->init();
        modify->fix[rfix]->setup_pre_neighbor();
        modify->fix[rfix]->zero_momentum();
      } else if (strstr(modify->fix[rfix]->style,"rigid")) {
        modify->fix[rfix]->zero_momentum();
      } else error->all(FLERR,"Velocity rigid used with non-rigid fix-ID");
    }

  } else if (strcmp(arg[0],"angular") == 0) {
    if (rfix < 0) zero_rotation();
    else {
      if (strcmp(modify->fix[rfix]->style,"rigid/small") == 0) {
        lmp->init();
        modify->fix[rfix]->setup_pre_neighbor();
        modify->fix[rfix]->zero_rotation();
      } else if (strstr(modify->fix[rfix]->style,"rigid")) {
        modify->fix[rfix]->zero_rotation();
      } else error->all(FLERR,"Velocity rigid used with non-rigid fix-ID");
    }

  } else error->all(FLERR,"Illegal velocity command");
}

/* ----------------------------------------------------------------------
   rescale velocities of group atoms to t_new from t_old
   no bias applied here, since done in create() and scale()
------------------------------------------------------------------------- */

void Velocity::rescale(double t_old, double t_new)
{
  if (t_old == 0.0) error->all(FLERR,"Attempting to rescale a 0.0 temperature");

  double factor = sqrt(t_new/t_old);

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] *= factor;
      v[i][1] *= factor;
      v[i][2] *= factor;
    }
}

/* ----------------------------------------------------------------------
   zero the linear momentum of a group of atoms by adjusting v by -Vcm
------------------------------------------------------------------------- */

void Velocity::zero_momentum()
{
  // cannot have no atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Cannot zero momentum of no atoms");

  // compute velocity of center-of-mass of group

  double masstotal = group->mass(igroup);
  double vcm[3];
  group->vcm(igroup,masstotal,vcm);

  // adjust velocities by vcm to zero linear momentum

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vcm[0];
      v[i][1] -= vcm[1];
      v[i][2] -= vcm[2];
    }
}

/* ----------------------------------------------------------------------
   zero the angular momentum of a group of atoms by adjusting v by -(w x r)
------------------------------------------------------------------------- */

void Velocity::zero_rotation()
{
  int i;

  // cannot have no atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Cannot zero momentum of no atoms");

  // compute omega (angular velocity) of group around center-of-mass

  double xcm[3],angmom[3],inertia[3][3],omega[3];
  double masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  group->angmom(igroup,xcm,angmom);
  group->inertia(igroup,xcm,inertia);
  group->omega(angmom,inertia,omega);

  // adjust velocities to zero omega
  // vnew_i = v_i - w x r_i
  // must use unwrapped coords to compute r_i correctly

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      v[i][0] -= omega[1]*dz - omega[2]*dy;
      v[i][1] -= omega[2]*dx - omega[0]*dz;
      v[i][2] -= omega[0]*dy - omega[1]*dx;
    }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of velocity input line
------------------------------------------------------------------------- */

void Velocity::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal velocity command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dist") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"uniform") == 0) dist_flag = 0;
      else if (strcmp(arg[iarg+1],"gaussian") == 0) dist_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) sum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) sum_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) momentum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) momentum_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rot") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) rotation_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) rotation_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      int icompute;
      for (icompute = 0; icompute < modify->ncompute; icompute++)
        if (strcmp(arg[iarg+1],modify->compute[icompute]->id) == 0) break;
      if (icompute == modify->ncompute)
        error->all(FLERR,"Could not find velocity temperature ID");
      temperature = modify->compute[icompute];
      if (temperature->tempflag == 0)
        error->all(FLERR,
                   "Velocity temperature ID does not compute temperature");
      iarg += 2;
    } else if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) bias_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) bias_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"all") == 0) loop_flag = ALL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop_flag = LOCAL;
      else if (strcmp(arg[iarg+1],"geom") == 0) loop_flag = GEOM;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rigid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      rfix = modify->find_fix(arg[iarg+1]);
      if (rfix < 0) error->all(FLERR,"Fix ID for velocity does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"box") == 0) scale_flag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scale_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else error->all(FLERR,"Illegal velocity command");
  }

  // error check

  if (bias_flag && temperature == NULL)
    error->all(FLERR,"Cannot use velocity bias command without temp keyword");
  if (bias_flag && temperature->tempbias == 0)
    error->all(FLERR,"Velocity temperature ID does calculate a velocity bias");
}
