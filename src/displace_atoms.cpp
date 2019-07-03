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

#include "displace_atoms.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "irregular.h"
#include "group.h"
#include "math_const.h"
#include "random_park.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "atom_vec_body.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{MOVE,RAMP,RANDOM,ROTATE};

/* ---------------------------------------------------------------------- */

DisplaceAtoms::DisplaceAtoms(LAMMPS *lmp) : Pointers(lmp)
{
  mvec = NULL;
}

/* ---------------------------------------------------------------------- */

DisplaceAtoms::~DisplaceAtoms()
{
  memory->destroy(mvec);
}

/* ---------------------------------------------------------------------- */

void DisplaceAtoms::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace_atoms command");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot displace_atoms after "
               "reading restart file with per-atom info");

  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_atoms group ID");
  groupbit = group->bitmask[igroup];

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Attempting to displace atoms in rigid bodies");

  int style = -1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  else if (strcmp(arg[1],"rotate") == 0) style = ROTATE;
  else error->all(FLERR,"Illegal displace_atoms command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == RANDOM) options(narg-6,&arg[6]);
  else if (style == ROTATE) options(narg-9,&arg[9]);

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector or specified variable(s)

  if (style == MOVE) {
    move(0,arg[2],xscale);
    move(1,arg[3],yscale);
    move(2,arg[4],zscale);
  }

  // move atoms in ramped fashion

  if (style == RAMP) {

    int d_dim = 0;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*force->numeric(FLERR,arg[3]);
      d_hi = xscale*force->numeric(FLERR,arg[4]);
    } else if (d_dim == 1) {
      d_lo = yscale*force->numeric(FLERR,arg[3]);
      d_hi = yscale*force->numeric(FLERR,arg[4]);
    } else if (d_dim == 2) {
      d_lo = zscale*force->numeric(FLERR,arg[3]);
      d_hi = zscale*force->numeric(FLERR,arg[4]);
    }

    int coord_dim = 0;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*force->numeric(FLERR,arg[6]);
      coord_hi = xscale*force->numeric(FLERR,arg[7]);
    } else if (coord_dim == 1) {
      coord_lo = yscale*force->numeric(FLERR,arg[6]);
      coord_hi = yscale*force->numeric(FLERR,arg[7]);
    } else if (coord_dim == 2) {
      coord_lo = zscale*force->numeric(FLERR,arg[6]);
      coord_hi = zscale*force->numeric(FLERR,arg[7]);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double fraction,dramp;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
        fraction = MAX(fraction,0.0);
        fraction = MIN(fraction,1.0);
        dramp = d_lo + fraction*(d_hi - d_lo);
        x[i][d_dim] += dramp;
      }
    }
  }

  // move atoms randomly
  // makes atom result independent of what proc owns it via random->reset()

  if (style == RANDOM) {
    RanPark *random = new RanPark(lmp,1);

    double dx = xscale*force->numeric(FLERR,arg[2]);
    double dy = yscale*force->numeric(FLERR,arg[3]);
    double dz = zscale*force->numeric(FLERR,arg[4]);
    int seed = force->inumeric(FLERR,arg[5]);
    if (seed <= 0) error->all(FLERR,"Illegal displace_atoms random command");

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        random->reset(seed,x[i]);
        x[i][0] += dx * 2.0*(random->uniform()-0.5);
        x[i][1] += dy * 2.0*(random->uniform()-0.5);
        x[i][2] += dz * 2.0*(random->uniform()-0.5);
      }
    }

    delete random;
  }

  // rotate atoms by right-hand rule by theta around R
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // R0 = runit = unit vector for R
  // D = X - P = vector from P to X
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(theta) + B sin(theta)

  if (style == ROTATE) {
    double theta_new;
    double axis[3],point[3],qrotate[4],qnew[4];
    double a[3],b[3],c[3],d[3],disp[3],runit[3];
    double *quat;

    int dim = domain->dimension;
    point[0] = xscale*force->numeric(FLERR,arg[2]);
    point[1] = yscale*force->numeric(FLERR,arg[3]);
    point[2] = zscale*force->numeric(FLERR,arg[4]);
    axis[0] = force->numeric(FLERR,arg[5]);
    axis[1] = force->numeric(FLERR,arg[6]);
    axis[2] = force->numeric(FLERR,arg[7]);
    double theta = force->numeric(FLERR,arg[8]);
    if (dim == 2 && (axis[0] != 0.0 || axis[1] != 0.0))
      error->all(FLERR,"Invalid displace_atoms rotate axis for 2d");

    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Zero length rotation vector with displace_atoms");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;

    double angle = MY_PI*theta/180.0;
    double cosine = cos(angle);
    double sine = sin(angle);

    double qcosine = cos(0.5*angle);
    double qsine = sin(0.5*angle);
    qrotate[0] = qcosine;
    qrotate[1] = runit[0]*qsine;
    qrotate[2] = runit[1]*qsine;
    qrotate[3] = runit[2]*qsine;

    double ddotr;

    // flags for additional orientation info stored by some atom styles

    int ellipsoid_flag = atom->ellipsoid_flag;
    int line_flag = atom->line_flag;
    int tri_flag = atom->tri_flag;
    int body_flag = atom->body_flag;

    int theta_flag = 0;
    int quat_flag = 0;
    if (line_flag) theta_flag = 1;
    if (ellipsoid_flag || tri_flag || body_flag) quat_flag = 1;

    // AtomVec pointers to retrieve per-atom storage of extra quantities

    AtomVecEllipsoid *avec_ellipsoid =
      (AtomVecEllipsoid *) atom->style_match("ellipsoid");
    AtomVecLine *avec_line = (AtomVecLine *) atom->style_match("line");
    AtomVecTri *avec_tri = (AtomVecTri *) atom->style_match("tri");
    AtomVecBody *avec_body = (AtomVecBody *) atom->style_match("body");

    double **x = atom->x;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;
    int *body = atom->body;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    imageint *image = atom->image;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // unwrap coordinate and reset image flags accordingly
        domain->unmap(x[i],image[i]);
        image[i] = ((imageint) IMGMAX << IMG2BITS) |
          ((imageint) IMGMAX << IMGBITS) | IMGMAX;

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
        x[i][0] = point[0] + c[0] + disp[0];
        x[i][1] = point[1] + c[1] + disp[1];
        if (dim == 3) x[i][2] = point[2] + c[2] + disp[2];

        // theta for lines

        if (theta_flag && line[i] >= 0.0) {
          theta_new = fmod(avec_line->bonus[line[i]].theta+angle,MY_2PI);
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
          if (quat) {
            MathExtra::quatquat(qrotate,quat,qnew);
            quat[0] = qnew[0];
            quat[1] = qnew[1];
            quat[2] = qnew[2];
            quat[3] = qnew[3];
          }
        }
      }
    }
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0) {
    char str[128];
    sprintf(str,"Lost atoms via displace_atoms: original " BIGINT_FORMAT
            " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   move atoms either by specified numeric displacement or variable evaluation
------------------------------------------------------------------------- */

void DisplaceAtoms::move(int idim, char *arg, double scale)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (strstr(arg,"v_") != arg) {
    double delta = scale*force->numeric(FLERR,arg);
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) x[i][idim] += delta;

  } else {
    int ivar = input->variable->find(arg+2);
    if (ivar < 0)
      error->all(FLERR,"Variable name for displace_atoms does not exist");

    if (input->variable->equalstyle(ivar)) {
      double delta = scale * input->variable->compute_equal(ivar);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) x[i][idim] += delta;
    } else if (input->variable->atomstyle(ivar)) {
      if (mvec == NULL) memory->create(mvec,nlocal,"displace_atoms:mvec");
      input->variable->compute_atom(ivar,igroup,mvec,1,0);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) x[i][idim] += scale*mvec[i];
    } else error->all(FLERR,"Variable for displace_atoms is invalid style");
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_atoms input line
------------------------------------------------------------------------- */

void DisplaceAtoms::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace_atoms command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace_atoms command");
  }
}
