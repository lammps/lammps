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

#include "lmptype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "displace_atoms.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "irregular.h"
#include "group.h"
#include "math_const.h"
#include "random_park.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{MOVE,RAMP,RANDOM,ROTATE};

/* ---------------------------------------------------------------------- */

DisplaceAtoms::DisplaceAtoms(LAMMPS *lmp) : Pointers(lmp) {}

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

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_atoms group ID");
  int groupbit = group->bitmask[igroup];

  int style;
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

  // move atoms by 3-vector

  if (style == MOVE) {

    double delx = xscale*force->numeric(FLERR,arg[2]);
    double dely = yscale*force->numeric(FLERR,arg[3]);
    double delz = zscale*force->numeric(FLERR,arg[4]);

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][0] += delx;
        x[i][1] += dely;
        x[i][2] += delz;
      }
    }
  }

  // move atoms in ramped fashion

  if (style == RAMP) {

    int d_dim;
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

    int coord_dim;
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
    double axis[3],point[3];
    double a[3],b[3],c[3],d[3],disp[3],runit[3];
    
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

    double sine = sin(MY_PI*theta/180.0);
    double cosine = cos(MY_PI*theta/180.0);
    double ddotr;

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
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
      }
    }
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
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
