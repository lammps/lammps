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
#include "fix_move.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{LINEAR,WIGGLE,ROTATE,VARIABLE};
enum{EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixMove::FixMove(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix move command");

  restart_global = 1;
  restart_peratom = 1;
  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;
  time_integrate = 1;
  create_attribute = 1;

  // parse args

  int iarg;
  xvarstr = yvarstr = zvarstr = NULL;
  vxvarstr = vyvarstr = vzvarstr = NULL;

  if (strcmp(arg[3],"linear") == 0) {
    if (narg < 7) error->all(FLERR,"Illegal fix move command");
    iarg = 7;
    mstyle = LINEAR;
    if (strcmp(arg[4],"NULL") == 0) vxflag = 0;
    else {
      vxflag = 1;
      vx = force->numeric(FLERR,arg[4]);
    }
    if (strcmp(arg[5],"NULL") == 0) vyflag = 0;
    else {
      vyflag = 1;
      vy = force->numeric(FLERR,arg[5]);
    }
    if (strcmp(arg[6],"NULL") == 0) vzflag = 0;
    else {
      vzflag = 1;
      vz = force->numeric(FLERR,arg[6]);
    }

  } else if (strcmp(arg[3],"wiggle") == 0) {
    if (narg < 8) error->all(FLERR,"Illegal fix move command");
    iarg = 8;
    mstyle = WIGGLE;
    if (strcmp(arg[4],"NULL") == 0) axflag = 0;
    else {
      axflag = 1;
      ax = force->numeric(FLERR,arg[4]);
    }
    if (strcmp(arg[5],"NULL") == 0) ayflag = 0;
    else {
      ayflag = 1;
      ay = force->numeric(FLERR,arg[5]);
    }
    if (strcmp(arg[6],"NULL") == 0) azflag = 0;
    else {
      azflag = 1;
      az = force->numeric(FLERR,arg[6]);
    }
    period = force->numeric(FLERR,arg[7]);

  } else if (strcmp(arg[3],"rotate") == 0) {
    if (narg < 11) error->all(FLERR,"Illegal fix move command");
    iarg = 11;
    mstyle = ROTATE;
    point[0] = force->numeric(FLERR,arg[4]);
    point[1] = force->numeric(FLERR,arg[5]);
    point[2] = force->numeric(FLERR,arg[6]);
    axis[0] = force->numeric(FLERR,arg[7]);
    axis[1] = force->numeric(FLERR,arg[8]);
    axis[2] = force->numeric(FLERR,arg[9]);
    period = force->numeric(FLERR,arg[10]);

  } else if (strcmp(arg[3],"variable") == 0) {
    if (narg < 10) error->all(FLERR,"Illegal fix move command");
    iarg = 10;
    mstyle = VARIABLE;
    if (strcmp(arg[4],"NULL") == 0) xvarstr = NULL;
    else if (strstr(arg[4],"v_") == arg[4]) {
      int n = strlen(&arg[4][2]) + 1;
      xvarstr = new char[n];
      strcpy(xvarstr,&arg[4][2]);
    } else error->all(FLERR,"Illegal fix move command");
    if (strcmp(arg[5],"NULL") == 0) yvarstr = NULL;
    else if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      yvarstr = new char[n];
      strcpy(yvarstr,&arg[5][2]);
    } else error->all(FLERR,"Illegal fix move command");
    if (strcmp(arg[6],"NULL") == 0) zvarstr = NULL;
    else if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      zvarstr = new char[n];
      strcpy(zvarstr,&arg[6][2]);
    } else error->all(FLERR,"Illegal fix move command");
    if (strcmp(arg[7],"NULL") == 0) vxvarstr = NULL;
    else if (strstr(arg[7],"v_") == arg[7]) {
      int n = strlen(&arg[7][2]) + 1;
      vxvarstr = new char[n];
      strcpy(vxvarstr,&arg[7][2]);
    } else error->all(FLERR,"Illegal fix move command");
    if (strcmp(arg[8],"NULL") == 0) vyvarstr = NULL;
    else if (strstr(arg[8],"v_") == arg[8]) {
      int n = strlen(&arg[8][2]) + 1;
      vyvarstr = new char[n];
      strcpy(vyvarstr,&arg[8][2]);
    } else error->all(FLERR,"Illegal fix move command");
    if (strcmp(arg[9],"NULL") == 0) vzvarstr = NULL;
    else if (strstr(arg[9],"v_") == arg[9]) {
      int n = strlen(&arg[9][2]) + 1;
      vzvarstr = new char[n];
      strcpy(vzvarstr,&arg[9][2]);
    } else error->all(FLERR,"Illegal fix move command");

  } else error->all(FLERR,"Illegal fix move command");

  // optional args

  int scaleflag = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix move command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix move command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix move command");
  }

  // error checks and warnings

  if (domain->dimension == 2) {
    if (mstyle == LINEAR && vzflag && vz != 0.0)
      error->all(FLERR,"Fix move cannot set linear z motion for 2d problem");
    if (mstyle == WIGGLE && azflag && az != 0.0)
      error->all(FLERR,"Fix move cannot set wiggle z motion for 2d problem");
    if (mstyle == ROTATE && (axis[0] != 0.0 || axis[1] != 0.0))
      error->all(FLERR,
                 "Fix move cannot rotate aroung non z-axis for 2d problem");
    if (mstyle == VARIABLE && (zvarstr || vzvarstr))
      error->all(FLERR,
                 "Fix move cannot define z or vz variable for 2d problem");
  }

  if (atom->angmom_flag && comm->me == 0)
    error->warning(FLERR,"Fix move does not update angular momentum");
  if (atom->ellipsoid_flag && comm->me == 0)
    error->warning(FLERR,"Fix move does not update quaternions");

  // setup scaling and apply scaling factors to velocity & amplitude

  if ((mstyle == LINEAR || mstyle == WIGGLE || mstyle == ROTATE) &&
      scaleflag) {
    double xscale,yscale,zscale;
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    if (mstyle == LINEAR) {
      if (vxflag) vx *= xscale;
      if (vyflag) vy *= yscale;
      if (vzflag) vz *= zscale;
    } else if (mstyle == WIGGLE) {
      if (axflag) ax *= xscale;
      if (ayflag) ay *= yscale;
      if (azflag) az *= zscale;
    } else if (mstyle == ROTATE) {
      point[0] *= xscale;
      point[1] *= yscale;
      point[2] *= zscale;
    }
  }

  // set omega_rotate from period

  if (mstyle == WIGGLE || mstyle == ROTATE) omega_rotate = 2.0*MY_PI / period;

  // runit = unit vector along rotation axis

  if (mstyle == ROTATE) {
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Zero length rotation vector with fix move");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }

  // set omega_flag if particles store omega

  omega_flag = atom->omega_flag;

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  maxatom = 0;
  displace = velocity = NULL;

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixMove::~FixMove()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(xoriginal);
  memory->destroy(displace);
  memory->destroy(velocity);

  delete [] xvarstr;
  delete [] yvarstr;
  delete [] zvarstr;
  delete [] vxvarstr;
  delete [] vyvarstr;
  delete [] vzvarstr;
}

/* ---------------------------------------------------------------------- */

int FixMove::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMove::init()
{
  dt = update->dt;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // set indices and style of all variables

  displaceflag = velocityflag = 0;

  if (mstyle == VARIABLE) {
    if (xvarstr) {
      xvar = input->variable->find(xvarstr);
      if (xvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(xvar)) xvarstyle = EQUAL;
      else if (input->variable->atomstyle(xvar)) xvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (yvarstr) {
      yvar = input->variable->find(yvarstr);
      if (yvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(yvar)) yvarstyle = EQUAL;
      else if (input->variable->atomstyle(yvar)) yvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (zvarstr) {
      zvar = input->variable->find(zvarstr);
      if (zvar < 0) error->all(FLERR,
                               "Variable name for fix move does not exist");
      if (input->variable->equalstyle(zvar)) zvarstyle = EQUAL;
      else if (input->variable->atomstyle(zvar)) zvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (vxvarstr) {
      vxvar = input->variable->find(vxvarstr);
      if (vxvar < 0) error->all(FLERR,
                                "Variable name for fix move does not exist");
      if (input->variable->equalstyle(vxvar)) vxvarstyle = EQUAL;
      else if (input->variable->atomstyle(vxvar)) vxvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (vyvarstr) {
      vyvar = input->variable->find(vyvarstr);
      if (vyvar < 0) error->all(FLERR,
                                "Variable name for fix move does not exist");
      if (input->variable->equalstyle(vyvar)) vyvarstyle = EQUAL;
      else if (input->variable->atomstyle(vyvar)) vyvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }
    if (vzvarstr) {
      vzvar = input->variable->find(vzvarstr);
      if (vzvar < 0) error->all(FLERR,
                                "Variable name for fix move does not exist");
      if (input->variable->equalstyle(vzvar)) vzvarstyle = EQUAL;
      else if (input->variable->atomstyle(vzvar)) vzvarstyle = ATOM;
      else error->all(FLERR,"Variable for fix move is invalid style");
    }

    if (xvarstr && xvarstyle == ATOM) displaceflag = 1;
    if (yvarstr && yvarstyle == ATOM) displaceflag = 1;
    if (zvarstr && zvarstyle == ATOM) displaceflag = 1;
    if (vxvarstr && vxvarstyle == ATOM) velocityflag = 1;
    if (vyvarstr && vyvarstyle == ATOM) velocityflag = 1;
    if (vzvarstr && vzvarstyle == ATOM) velocityflag = 1;
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ----------------------------------------------------------------------
   set x,v of particles
------------------------------------------------------------------------- */

void FixMove::initial_integrate(int vflag)
{
  double dtfm;
  double xold[3],a[3],b[3],c[3],d[3],disp[3];
  double ddotr,dx,dy,dz;

  double delta = (update->ntimestep - time_origin) * dt;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // for linear: X = X0 + V*dt

  if (mstyle == LINEAR) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        xold[0] = x[i][0];
        xold[1] = x[i][1];
        xold[2] = x[i][2];

        if (vxflag) {
          v[i][0] = vx;
          x[i][0] = xoriginal[i][0] + vx*delta;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          x[i][0] += dtv * v[i][0];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          x[i][0] += dtv * v[i][0];
        }

        if (vyflag) {
          v[i][1] = vy;
          x[i][1] = xoriginal[i][1] + vy*delta;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][1] += dtfm * f[i][1];
          x[i][1] += dtv * v[i][1];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][1] += dtfm * f[i][1];
          x[i][1] += dtv * v[i][1];
        }

        if (vzflag) {
          v[i][2] = vz;
          x[i][2] = xoriginal[i][2] + vz*delta;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][2] += dtfm * f[i][2];
          x[i][2] += dtv * v[i][2];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][2] += dtfm * f[i][2];
          x[i][2] += dtv * v[i][2];
        }

        domain->remap_near(x[i],xold);
      }
    }

  // for wiggle: X = X0 + A sin(w*dt)

  } else if (mstyle == WIGGLE) {
    double arg = omega_rotate * delta;
    double sine = sin(arg);
    double cosine = cos(arg);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        xold[0] = x[i][0];
        xold[1] = x[i][1];
        xold[2] = x[i][2];

        if (axflag) {
          v[i][0] = ax*omega_rotate*cosine;
          x[i][0] = xoriginal[i][0] + ax*sine;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          x[i][0] += dtv * v[i][0];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          x[i][0] += dtv * v[i][0];
        }

        if (ayflag) {
          v[i][1] = ay*omega_rotate*cosine;
          x[i][1] = xoriginal[i][1] + ay*sine;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][1] += dtfm * f[i][1];
          x[i][1] += dtv * v[i][1];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][1] += dtfm * f[i][1];
          x[i][1] += dtv * v[i][1];
        }

        if (azflag) {
          v[i][2] = az*omega_rotate*cosine;
          x[i][2] = xoriginal[i][2] + az*sine;
        } else if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][2] += dtfm * f[i][2];
          x[i][2] += dtv * v[i][2];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][2] += dtfm * f[i][2];
          x[i][2] += dtv * v[i][2];
        }

        domain->remap_near(x[i],xold);
      }
    }

  // for rotate by right-hand rule around omega:
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // w = omega of rotation (from period)
  // X0 = xoriginal = initial coord of atom
  // R0 = runit = unit vector for R
  // D = X0 - P = vector from P to X0
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X0
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(w*dt) + B sin(w*dt)
  // V = w R0 cross (A cos(w*dt) + B sin(w*dt))

  } else if (mstyle == ROTATE) {
    double arg = omega_rotate * delta;
    double sine = sin(arg);
    double cosine = cos(arg);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        xold[0] = x[i][0];
        xold[1] = x[i][1];
        xold[2] = x[i][2];

        d[0] = xoriginal[i][0] - point[0];
        d[1] = xoriginal[i][1] - point[1];
        d[2] = xoriginal[i][2] - point[2];
        ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
        c[0] = ddotr*runit[0];
        c[1] = ddotr*runit[1];
        c[2] = ddotr*runit[2];
        a[0] = d[0] - c[0];
        a[1] = d[1] - c[1];
        a[2] = d[2] - c[2];
        b[0] = runit[1]*a[2] - runit[2]*a[1];
        b[1] = runit[2]*a[0] - runit[0]*a[2];
        b[2] = runit[0]*a[1] - runit[1]*a[0];
        disp[0] = a[0]*cosine  + b[0]*sine;
        disp[1] = a[1]*cosine  + b[1]*sine;
        disp[2] = a[2]*cosine  + b[2]*sine;

        x[i][0] = point[0] + c[0] + disp[0];
        x[i][1] = point[1] + c[1] + disp[1];
        x[i][2] = point[2] + c[2] + disp[2];
        v[i][0] = omega_rotate * (runit[1]*disp[2] - runit[2]*disp[1]);
        v[i][1] = omega_rotate * (runit[2]*disp[0] - runit[0]*disp[2]);
        v[i][2] = omega_rotate * (runit[0]*disp[1] - runit[1]*disp[0]);
        if (omega_flag) {
          omega[i][0] = omega_rotate*runit[0];
          omega[i][1] = omega_rotate*runit[1];
          omega[i][2] = omega_rotate*runit[2];
        }

        domain->remap_near(x[i],xold);
      }
    }

  // for variable: compute x,v from variables

  } else if (mstyle == VARIABLE) {

    // reallocate displace and velocity arrays as necessary

    if ((displaceflag || velocityflag) && nlocal > maxatom) {
      maxatom = atom->nmax;
      if (displaceflag) {
        memory->destroy(displace);
        memory->create(displace,maxatom,3,"move:displace");
      }
      if (velocityflag) {
        memory->destroy(velocity);
        memory->create(velocity,maxatom,3,"move:velocity");
      }
    }

    // pre-compute variable values, wrap with clear/add

    modify->clearstep_compute();

    if (xvarstr) {
      if (xvarstyle == EQUAL) dx = input->variable->compute_equal(xvar);
      else if (displace)
        input->variable->compute_atom(xvar,igroup,&displace[0][0],3,0);
    }
    if (yvarstr) {
      if (yvarstyle == EQUAL) dy = input->variable->compute_equal(yvar);
      else if (displace)
        input->variable->compute_atom(yvar,igroup,&displace[0][1],3,0);
    }
    if (zvarstr) {
      if (zvarstyle == EQUAL) dz = input->variable->compute_equal(zvar);
      else if (displace)
        input->variable->compute_atom(zvar,igroup,&displace[0][2],3,0);
    }
    if (vxvarstr) {
      if (vxvarstyle == EQUAL) vx = input->variable->compute_equal(vxvar);
      else if (velocity)
        input->variable->compute_atom(vxvar,igroup,&velocity[0][0],3,0);
    }
    if (vyvarstr) {
      if (vyvarstyle == EQUAL) vy = input->variable->compute_equal(vyvar);
      else if (velocity)
        input->variable->compute_atom(vyvar,igroup,&velocity[0][1],3,0);
    }
    if (vzvarstr) {
      if (vzvarstyle == EQUAL) vz = input->variable->compute_equal(vzvar);
      else if (velocity)
        input->variable->compute_atom(vzvar,igroup,&velocity[0][2],3,0);
    }

    modify->addstep_compute(update->ntimestep + 1);

    // update x,v

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        xold[0] = x[i][0];
        xold[1] = x[i][1];
        xold[2] = x[i][2];

        if (xvarstr && vxvarstr) {
          if (vxvarstyle == EQUAL) v[i][0] = vx;
          else v[i][0] = velocity[i][0];
          if (xvarstyle == EQUAL) x[i][0] = xoriginal[i][0] + dx;
          else x[i][0] = xoriginal[i][0] + displace[i][0];
        } else if (xvarstr) {
          if (xvarstyle == EQUAL) x[i][0] = xoriginal[i][0] + dx;
          else x[i][0] = xoriginal[i][0] + displace[i][0];
        } else if (vxvarstr) {
          if (vxvarstyle == EQUAL) v[i][0] = vx;
          else v[i][0] = velocity[i][0];
          if (rmass) {
            dtfm = dtf / rmass[i];
            x[i][0] += dtv * v[i][0];
          } else {
            dtfm = dtf / mass[type[i]];
            x[i][0] += dtv * v[i][0];
          }
        } else {
          if (rmass) {
            dtfm = dtf / rmass[i];
            v[i][0] += dtfm * f[i][0];
            x[i][0] += dtv * v[i][0];
          } else {
            dtfm = dtf / mass[type[i]];
            v[i][0] += dtfm * f[i][0];
            x[i][0] += dtv * v[i][0];
          }
        }

        if (yvarstr && vyvarstr) {
          if (vyvarstyle == EQUAL) v[i][1] = vy;
          else v[i][1] = velocity[i][1];
          if (yvarstyle == EQUAL) x[i][1] = xoriginal[i][1] + dy;
          else x[i][1] = xoriginal[i][1] + displace[i][1];
        } else if (yvarstr) {
          if (yvarstyle == EQUAL) x[i][1] = xoriginal[i][1] + dy;
          else x[i][1] = xoriginal[i][1] + displace[i][1];
        } else if (vyvarstr) {
          if (vyvarstyle == EQUAL) v[i][1] = vy;
          else v[i][1] = velocity[i][1];
          if (rmass) {
            dtfm = dtf / rmass[i];
            x[i][1] += dtv * v[i][1];
          } else {
            dtfm = dtf / mass[type[i]];
            x[i][1] += dtv * v[i][1];
          }
        } else {
          if (rmass) {
            dtfm = dtf / rmass[i];
            v[i][1] += dtfm * f[i][1];
            x[i][1] += dtv * v[i][1];
          } else {
            dtfm = dtf / mass[type[i]];
            v[i][1] += dtfm * f[i][1];
            x[i][1] += dtv * v[i][1];
          }
        }

        if (zvarstr && vzvarstr) {
          if (vzvarstyle == EQUAL) v[i][2] = vz;
          else v[i][2] = velocity[i][2];
          if (zvarstyle == EQUAL) x[i][2] = xoriginal[i][2] + dz;
          else x[i][2] = xoriginal[i][2] + displace[i][2];
        } else if (zvarstr) {
          if (zvarstyle == EQUAL) x[i][2] = xoriginal[i][2] + dz;
          else x[i][2] = xoriginal[i][2] + displace[i][2];
        } else if (vzvarstr) {
          if (vzvarstyle == EQUAL) v[i][2] = vz;
          else v[i][2] = velocity[i][2];
          if (rmass) {
            dtfm = dtf / rmass[i];
            x[i][2] += dtv * v[i][2];
          } else {
            dtfm = dtf / mass[type[i]];
            x[i][2] += dtv * v[i][2];
          }
        } else {
          if (rmass) {
            dtfm = dtf / rmass[i];
            v[i][2] += dtfm * f[i][2];
            x[i][2] += dtv * v[i][2];
          } else {
            dtfm = dtf / mass[type[i]];
            v[i][2] += dtfm * f[i][2];
            x[i][2] += dtv * v[i][2];
          }
        }

        domain->remap_near(x[i],xold);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   final NVE of particles with NULL components
------------------------------------------------------------------------- */

void FixMove::final_integrate()
{
  double dtfm;

  int xflag = 1;
  if (mstyle == LINEAR && vxflag) xflag = 0;
  else if (mstyle == WIGGLE && axflag) xflag = 0;
  else if (mstyle == ROTATE) xflag = 0;
  else if (mstyle == VARIABLE && (xvarstr || vxvarstr)) xflag = 0;

  int yflag = 1;
  if (mstyle == LINEAR && vyflag) yflag = 0;
  else if (mstyle == WIGGLE && ayflag) yflag = 0;
  else if (mstyle == ROTATE) yflag = 0;
  else if (mstyle == VARIABLE && (yvarstr || vyvarstr)) yflag = 0;

  int zflag = 1;
  if (mstyle == LINEAR && vzflag) zflag = 0;
  else if (mstyle == WIGGLE && azflag) zflag = 0;
  else if (mstyle == ROTATE) zflag = 0;
  else if (mstyle == VARIABLE && (zvarstr || vzvarstr)) zflag = 0;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (xflag)
        if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
        }

      if (yflag)
        if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][1] += dtfm * f[i][1];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][1] += dtfm * f[i][1];
        }

      if (zflag)
        if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][2] += dtfm * f[i][2];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][2] += dtfm * f[i][2];
        }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMove::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  // outermost level - update v and x
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMove::final_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) final_integrate();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMove::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  if (displaceflag) bytes += atom->nmax*3 * sizeof(double);
  if (velocityflag) bytes += atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMove::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = time_origin;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMove::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  time_origin = static_cast<int> (list[n++]);
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixMove::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"move:xoriginal");
  array_atom = xoriginal;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixMove::copy_arrays(int i, int j, int delflag)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixMove::set_arrays(int i)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;

  // particle not in group

  if (!(mask[i] & groupbit)) {
    xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
    return;
  }

  // current time still equal fix creation time

  if (update->ntimestep == time_origin) {
    domain->unmap(x[i],image[i],xoriginal[i]);
    return;
  }

  // backup particle to time_origin

  if (mstyle == VARIABLE)
    error->all(FLERR,"Cannot add atoms to fix move variable");

  domain->unmap(x[i],image[i],xoriginal[i]);
  double delta = (update->ntimestep - time_origin) * update->dt;

  if (mstyle == LINEAR) {
    if (vxflag) xoriginal[i][0] -= vx * delta;
    if (vyflag) xoriginal[i][1] -= vy * delta;
    if (vzflag) xoriginal[i][2] -= vz * delta;
  } else if (mstyle == WIGGLE) {
    double arg = omega_rotate * delta;
    double sine = sin(arg);
    if (axflag) xoriginal[i][0] -= ax*sine;
    if (ayflag) xoriginal[i][1] -= ay*sine;
    if (azflag) xoriginal[i][2] -= az*sine;
  } else if (mstyle == ROTATE) {
    double a[3],b[3],c[3],d[3],disp[3],ddotr;
    double arg = - omega_rotate * delta;
    double sine = sin(arg);
    double cosine = cos(arg);
    d[0] = x[i][0] - point[0];
    d[1] = x[i][1] - point[1];
    d[2] = x[i][2] - point[2];
    ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
    c[0] = ddotr*runit[0];
    c[1] = ddotr*runit[1];
    c[2] = ddotr*runit[2];

    a[0] = d[0] - c[0];
    a[1] = d[1] - c[1];
    a[2] = d[2] - c[2];
    b[0] = runit[1]*a[2] - runit[2]*a[1];
    b[1] = runit[2]*a[0] - runit[0]*a[2];
    b[2] = runit[0]*a[1] - runit[1]*a[0];
    disp[0] = a[0]*cosine  + b[0]*sine;
    disp[1] = a[1]*cosine  + b[1]*sine;
    disp[2] = a[2]*cosine  + b[2]*sine;

    xoriginal[i][0] = point[0] + c[0] + disp[0];
    xoriginal[i][1] = point[1] + c[1] + disp[1];
    xoriginal[i][2] = point[2] + c[2] + disp[2];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixMove::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixMove::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixMove::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixMove::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixMove::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixMove::size_restart(int nlocal)
{
  return 4;
}

/* ---------------------------------------------------------------------- */

void FixMove::reset_dt()
{
  error->all(FLERR,"Resetting timestep is not allowed with fix move");
}
