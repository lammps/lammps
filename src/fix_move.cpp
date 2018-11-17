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
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "atom_vec_body.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{LINEAR,WIGGLE,ROTATE,VARIABLE};
enum{EQUAL,ATOM};

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixMove::FixMove(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xvarstr(NULL), yvarstr(NULL), zvarstr(NULL), vxvarstr(NULL),
  vyvarstr(NULL), vzvarstr(NULL),
  xoriginal(NULL), toriginal(NULL), qoriginal(NULL),
  displace(NULL), velocity(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix move command");

  restart_global = 1;
  restart_peratom = 1;
  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;
  time_integrate = 1;
  create_attribute = 1;
  displaceflag = 0;
  velocityflag = 0;
  maxatom = 0;

  // parse args

  int iarg = 0;

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
    if (period <= 0.0) error->all(FLERR,"Illegal fix move command");

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
    if (period <= 0.0) error->all(FLERR,"Illegal fix move command");

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
                 "Fix move cannot rotate around non z-axis for 2d problem");
    if (mstyle == VARIABLE && (zvarstr || vzvarstr))
      error->all(FLERR,
                 "Fix move cannot define z or vz variable for 2d problem");
  }

  // setup scaling and apply scaling factors to velocity & amplitude

  if ((mstyle == LINEAR || mstyle == WIGGLE || mstyle == ROTATE) &&
      scaleflag) {
    double xscale = domain->lattice->xlattice;
    double yscale = domain->lattice->ylattice;
    double zscale = domain->lattice->zlattice;

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

  if (mstyle == WIGGLE || mstyle == ROTATE) omega_rotate = MY_2PI / period;

  // runit = unit vector along rotation axis

  if (mstyle == ROTATE) {
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Zero length rotation vector with fix move");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }

  // set flags for extra attributes particles may store
  // relevant extra attributes = omega, angmom, theta, quat

  omega_flag = atom->omega_flag;
  angmom_flag = atom->angmom_flag;

  radius_flag = atom->radius_flag;
  ellipsoid_flag = atom->ellipsoid_flag;
  line_flag = atom->line_flag;
  tri_flag = atom->tri_flag;
  body_flag = atom->body_flag;

  theta_flag = quat_flag = 0;
  if (line_flag) theta_flag = 1;
  if (ellipsoid_flag || tri_flag || body_flag) quat_flag = 1;

  extra_flag = 0;
  if (omega_flag || angmom_flag || theta_flag || quat_flag) extra_flag = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  displace = velocity = NULL;

  // AtomVec pointers to retrieve per-atom storage of extra quantities

  avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  avec_line = (AtomVecLine *) atom->style_match("line");
  avec_tri = (AtomVecTri *) atom->style_match("tri");
  avec_body = (AtomVecBody *) atom->style_match("body");

  // xoriginal = initial unwrapped positions of atoms
  // toriginal = initial theta of lines
  // qoriginal = initial quat of extended particles

  double **x = atom->x;
  imageint *image = atom->image;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  if (theta_flag) {
    for (int i = 0; i < nlocal; i++) {
      if ((mask[i] & groupbit) && line[i] >= 0)
        toriginal[i] = avec_line->bonus[line[i]].theta;
      else toriginal[i] = 0.0;
    }
  }

  if (quat_flag) {
    double *quat;
    for (int i = 0; i < nlocal; i++) {
      quat = NULL;
      if (mask[i] & groupbit) {
        if (ellipsoid_flag && ellipsoid[i] >= 0)
          quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
        else if (tri_flag && tri[i] >= 0)
          quat = avec_tri->bonus[tri[i]].quat;
        else if (body_flag && body[i] >= 0)
          quat = avec_body->bonus[body[i]].quat;
      }
      if (quat) {
        qoriginal[i][0] = quat[0];
        qoriginal[i][1] = quat[1];
        qoriginal[i][2] = quat[2];
        qoriginal[i][3] = quat[3];
      } else qoriginal[i][0] = qoriginal[i][1] =
               qoriginal[i][2] = qoriginal[i][3] = 0.0;
    }
  }

  // nrestart = size of per-atom restart data
  // nrestart = 1 + xorig + torig + qorig

  nrestart = 4;
  if (theta_flag) nrestart++;
  if (quat_flag) nrestart += 4;

  // time origin for movement = current timestep

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
  memory->destroy(toriginal);
  memory->destroy(qoriginal);
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

  maxatom = atom->nmax;
  memory->destroy(displace);
  memory->destroy(velocity);
  if (displaceflag) memory->create(displace,maxatom,3,"move:displace");
  else displace = NULL;
  if (velocityflag) memory->create(velocity,maxatom,3,"move:velocity");
  else velocity = NULL;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ----------------------------------------------------------------------
   set x,v of particles
------------------------------------------------------------------------- */

void FixMove::initial_integrate(int /*vflag*/)
{
  int flag;
  double ddotr,dx,dy,dz;
  double dtfm,theta_new;
  double xold[3],a[3],b[3],c[3],d[3],disp[3],w[3],ex[3],ey[3],ez[3];
  double inertia_ellipsoid[3],qrotate[4];
  double *quat,*inertia,*shape;

  double delta = (update->ntimestep - time_origin) * dt;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *body = atom->body;
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
    double cosine = cos(arg);
    double sine = sin(arg);

    double qcosine = cos(0.5*arg);
    double qsine = sin(0.5*arg);
    qrotate[0] = qcosine;
    qrotate[1] = runit[0]*qsine;
    qrotate[2] = runit[1]*qsine;
    qrotate[3] = runit[2]*qsine;

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

        // set any extra attributes affected by rotation

        if (extra_flag) {

          // omega for spheres, lines, tris

          if (omega_flag) {
            flag = 0;
            if (radius_flag && radius[i] > 0.0) flag = 1;
            if (line_flag && line[i] >= 0.0) flag = 1;
            if (tri_flag && tri[i] >= 0.0) flag = 1;
            if (flag) {
              omega[i][0] = omega_rotate*runit[0];
              omega[i][1] = omega_rotate*runit[1];
              omega[i][2] = omega_rotate*runit[2];
            }
          }

          // angmom for ellipsoids, tris, and bodies

          if (angmom_flag) {
            quat = inertia = NULL;
            if (ellipsoid_flag && ellipsoid[i] >= 0) {
              quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
              shape = avec_ellipsoid->bonus[ellipsoid[i]].shape;
              inertia_ellipsoid[0] =
                INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
              inertia_ellipsoid[1] =
                INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
              inertia_ellipsoid[2] =
                INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
              inertia = inertia_ellipsoid;
            } else if (tri_flag && tri[i] >= 0) {
              quat = avec_tri->bonus[tri[i]].quat;
              inertia = avec_tri->bonus[tri[i]].inertia;
            } else if (body_flag && body[i] >= 0) {
              quat = avec_body->bonus[body[i]].quat;
              inertia = avec_body->bonus[body[i]].inertia;
            }
            if (quat) {
              w[0] = omega_rotate*runit[0];
              w[1] = omega_rotate*runit[1];
              w[2] = omega_rotate*runit[2];
              MathExtra::q_to_exyz(quat,ex,ey,ez);
              MathExtra::omega_to_angmom(w,ex,ey,ez,inertia,angmom[i]);
            }
          }

          // theta for lines

          if (theta_flag && line[i] >= 0.0) {
            theta_new = fmod(toriginal[i]+arg,MY_2PI);
            avec_line->bonus[atom->line[i]].theta = theta_new;
          }

          // quats for ellipsoids, tris, and bodies

          if (quat_flag) {
            quat = NULL;
            if (ellipsoid_flag && ellipsoid[i] >= 0)
              quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
            else if (tri_flag && tri[i] >= 0)
              quat = avec_tri->bonus[tri[i]].quat;
            else if (body_flag && body[i] >= 0)
              quat = avec_body->bonus[body[i]].quat;
            if (quat) MathExtra::quatquat(qrotate,qoriginal[i],quat);
          }
        }

        domain->remap_near(x[i],xold);
      }
    }

  // for variable: compute x,v from variables
  // NOTE: also allow for changes to extra attributes?
  //       omega, angmom, theta, quat
  //       only necessary if prescribed motion involves rotation

  } else if (mstyle == VARIABLE) {

    // reallocate displace and velocity arrays as necessary

    if ((displaceflag || velocityflag) && atom->nmax > maxatom) {
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
      else input->variable->compute_atom(xvar,igroup,&displace[0][0],3,0);
    }
    if (yvarstr) {
      if (yvarstyle == EQUAL) dy = input->variable->compute_equal(yvar);
      else input->variable->compute_atom(yvar,igroup,&displace[0][1],3,0);
    }
    if (zvarstr) {
      if (zvarstyle == EQUAL) dz = input->variable->compute_equal(zvar);
      else input->variable->compute_atom(zvar,igroup,&displace[0][2],3,0);
    }
    if (vxvarstr) {
      if (vxvarstyle == EQUAL) vx = input->variable->compute_equal(vxvar);
      else input->variable->compute_atom(vxvar,igroup,&velocity[0][0],3,0);
    }
    if (vyvarstr) {
      if (vyvarstyle == EQUAL) vy = input->variable->compute_equal(vyvar);
      else input->variable->compute_atom(vyvar,igroup,&velocity[0][1],3,0);
    }
    if (vzvarstr) {
      if (vzvarstyle == EQUAL) vz = input->variable->compute_equal(vzvar);
      else input->variable->compute_atom(vzvar,igroup,&velocity[0][2],3,0);
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
      if (xflag) {
        if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
        }
      }

      if (yflag) {
        if (rmass) {
          dtfm = dtf / rmass[i];
          v[i][1] += dtfm * f[i][1];
        } else {
          dtfm = dtf / mass[type[i]];
          v[i][1] += dtfm * f[i][1];
        }
      }

      if (zflag) {
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
}

/* ---------------------------------------------------------------------- */

void FixMove::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  // outermost level - update v and x
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMove::final_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) final_integrate();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMove::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  if (theta_flag) bytes += atom->nmax * sizeof(double);
  if (quat_flag) bytes += atom->nmax*4 * sizeof(double);
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
  if (theta_flag) memory->grow(toriginal,nmax,"move:toriginal");
  if (quat_flag) memory->grow(qoriginal,nmax,4,"move:qoriginal");
  array_atom = xoriginal;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixMove::copy_arrays(int i, int j, int /*delflag*/)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
  if (theta_flag) toriginal[j] = toriginal[i];
  if (quat_flag) {
    qoriginal[j][0] = qoriginal[i][0];
    qoriginal[j][1] = qoriginal[i][1];
    qoriginal[j][2] = qoriginal[i][2];
    qoriginal[j][3] = qoriginal[i][3];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixMove::set_arrays(int i)
{
  double theta;
  double *quat;

  double **x = atom->x;
  imageint *image = atom->image;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *body = atom->body;
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

    // set theta and quat extra attributes affected by rotation

    if (extra_flag) {

      // theta for lines

      if (theta_flag && line[i] >= 0.0) {
        theta = avec_line->bonus[atom->line[i]].theta;
        toriginal[i] = theta - 0.0;  // NOTE: edit this line
      }

      // quats for ellipsoids, tris, and bodies

      if (quat_flag) {
        quat = NULL;
        if (ellipsoid_flag && ellipsoid[i] >= 0)
          quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
        else if (tri_flag && tri[i] >= 0)
          quat = avec_tri->bonus[tri[i]].quat;
        else if (body_flag && body[i] >= 0)
          quat = avec_body->bonus[body[i]].quat;
        if (quat) {
          // qoriginal = f(quat,-delta);   // NOTE: edit this line
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixMove::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = xoriginal[i][0];
  buf[n++] = xoriginal[i][1];
  buf[n++] = xoriginal[i][2];
  if (theta_flag) buf[n++] = toriginal[i];
  if (quat_flag) {
    buf[n++] = qoriginal[i][0];
    buf[n++] = qoriginal[i][1];
    buf[n++] = qoriginal[i][2];
    buf[n++] = qoriginal[i][3];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixMove::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  xoriginal[nlocal][0] = buf[n++];
  xoriginal[nlocal][1] = buf[n++];
  xoriginal[nlocal][2] = buf[n++];
  if (theta_flag) toriginal[nlocal] = buf[n++];
  if (quat_flag) {
    qoriginal[nlocal][0] = buf[n++];
    qoriginal[nlocal][1] = buf[n++];
    qoriginal[nlocal][2] = buf[n++];
    qoriginal[nlocal][3] = buf[n++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixMove::pack_restart(int i, double *buf)
{
  int n = 1;
  buf[n++] = xoriginal[i][0];
  buf[n++] = xoriginal[i][1];
  buf[n++] = xoriginal[i][2];
  if (theta_flag) buf[n++] = toriginal[i];
  if (quat_flag) {
    buf[n++] = qoriginal[i][0];
    buf[n++] = qoriginal[i][1];
    buf[n++] = qoriginal[i][2];
    buf[n++] = qoriginal[i][3];
  }
  buf[0] = n;
  return n;
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
  if (theta_flag) toriginal[nlocal] = extra[nlocal][m++];
  if (quat_flag) {
    qoriginal[nlocal][0] = extra[nlocal][m++];
    qoriginal[nlocal][1] = extra[nlocal][m++];
    qoriginal[nlocal][2] = extra[nlocal][m++];
    qoriginal[nlocal][3] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixMove::maxsize_restart()
{
  return nrestart;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixMove::size_restart(int /*nlocal*/)
{
  return nrestart;
}

/* ---------------------------------------------------------------------- */

void FixMove::reset_dt()
{
  error->all(FLERR,"Resetting timestep size is not allowed with fix move");
}
