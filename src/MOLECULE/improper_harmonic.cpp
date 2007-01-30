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
#include "stdlib.h"
#include "improper_harmonic.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperHarmonic::ImproperHarmonic(LAMMPS *lmp) : Improper(lmp) {}

/* ---------------------------------------------------------------------- */

ImproperHarmonic::~ImproperHarmonic()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(k);
    memory->sfree(chi);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonic::compute(int eflag, int vflag)
{
  int n,i1,i2,i3,i4,type,factor;
  double rfactor;
  double v1x,v1y,v1z,v2x,v2y,v2z,v3x;
  double v3y,v3z,ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c,s,domega,a,a11,a22,a33,a12,a13,a23,sx1;
  double sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12;

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {

    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    if (newton_bond) factor = 4;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
      if (i3 < nlocal) factor++;
      if (i4 < nlocal) factor++;
      }
    rfactor = 0.25 * factor;

    // geometry of 4-body

    v1x = x[i2][0] - x[i1][0];
    v1y = x[i2][1] - x[i1][1];
    v1z = x[i2][2] - x[i1][2];
    domain->minimum_image(&v1x,&v1y,&v1z);

    v2x = x[i3][0] - x[i2][0];
    v2y = x[i3][1] - x[i2][1];
    v2z = x[i3][2] - x[i2][2];
    domain->minimum_image(&v2x,&v2y,&v2z);

    v3x = x[i4][0] - x[i3][0];
    v3y = x[i4][1] - x[i3][1];
    v3z = x[i4][2] - x[i3][2];
    domain->minimum_image(&v3x,&v3y,&v3z);

    ss1 = 1.0 / (v1x*v1x + v1y*v1y + v1z*v1z);
    ss2 = 1.0 / (v2x*v2x + v2y*v2y + v2z*v2z);
    ss3 = 1.0 / (v3x*v3x + v3y*v3y + v3z*v3z);
        
    r1 = sqrt(ss1);
    r2 = sqrt(ss2);
    r3 = sqrt(ss3);
        
    // sin and cos of angle
        
    c0 = -(v1x * v3x + v1y * v3y + v1z * v3z) * r1 * r3;
    c1 = -(v1x * v2x + v1y * v2y + v1z * v2z) * r1 * r2;
    c2 = -(v3x * v2x + v3y * v2y + v3z * v2z) * r3 * r2;

    s1 = 1.0 - c1*c1;
    if (s1 < SMALL) s1 = SMALL;
    s1 = 1.0 / s1;

    s2 = 1.0 - c2*c2;
    if (s2 < SMALL) s2 = SMALL;
    s2 = 1.0 / s2;

    s12 = sqrt(s1*s2);
    c = (c1*c2 + c0) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
	fprintf(screen,"Improper problem: %d %d %d %d %d %d\n",
		me,update->ntimestep,
		atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	fprintf(screen,"  1st atom: %d %g %g %g\n",
		me,x[i1][0],x[i1][1],x[i1][2]);
	fprintf(screen,"  2nd atom: %d %g %g %g\n",
		me,x[i2][0],x[i2][1],x[i2][2]);
	fprintf(screen,"  3rd atom: %d %g %g %g\n",
		me,x[i3][0],x[i3][1],x[i3][2]);
	fprintf(screen,"  4th atom: %d %g %g %g\n",
		me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;

    // force & energy

    domega = acos(c) - chi[type];
    a = k[type] * domega;

    if (eflag) energy += rfactor * a * domega;

    a = -a * 2.0/s;
    c = c * a;

    s12 = s12 * a;
    a11 = (-c * ss1 * s1);
    a22 = ss2 * (2.0 * c0 * s12 - c * (s1 + s2));
    a33 = (-c * ss3 * s2);
    a12 = r1 * r2 * (c1 * c * s1 + c2 * s12);
    a13 = r1 * r3 * s12;
    a23 = r2 * r3 * (-c2 * c * s2 - c1 * s12);

    sx1  = a12*v2x + a13*v3x - a11*v1x;
    sx2  = a22*v2x + a23*v3x - a12*v1x;
    sx12 = a23*v2x + a33*v3x - a13*v1x;
    sy1  = a12*v2y + a13*v3y - a11*v1y;
    sy2  = a22*v2y + a23*v3y - a12*v1y;
    sy12 = a23*v2y + a33*v3y - a13*v1y;
    sz1  = a12*v2z + a13*v3z - a11*v1z;
    sz2  = a22*v2z + a23*v3z - a12*v1z;
    sz12 = a23*v2z + a33*v3z - a13*v1z;

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= sx1;
      f[i1][1] -= sy1;
      f[i1][2] -= sz1;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += sx2 + sx1;
      f[i2][1] += sy2 + sy1;
      f[i2][2] += sz2 + sz1;
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += sx12 - sx2;
      f[i3][1] += sy12 - sy2;
      f[i3][2] += sz12 - sz2;
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] -= sx12;
      f[i4][1] -= sy12;
      f[i4][2] -= sz12;
    }

    // virial contribution

    if (vflag) {
      virial[0] += rfactor * (v1x*sx1 - v2x*sx2 - v3x*sx12);
      virial[1] += rfactor * (v1y*sy1 - v2y*sy2 - v3y*sy12);
      virial[2] += rfactor * (v1z*sz1 - v2z*sz2 - v3z*sz12);
      virial[3] += rfactor * (v1x*sy1 - v2x*sy2 - v3x*sy12);
      virial[4] += rfactor * (v1x*sz1 - v2x*sz2 - v3x*sz12);
      virial[5] += rfactor * (v1y*sz1 - v2y*sz2 - v3y*sz12);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  k = (double *) memory->smalloc((n+1)*sizeof(double),"improper:k");
  chi = (double *) memory->smalloc((n+1)*sizeof(double),"improper:chi");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperHarmonic::coeff(int which, int narg, char **arg)
{
  if (which != 0) error->all("Invalid coeffs for this improper style");
  if (narg != 3) error->all("Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = atof(arg[1]);
  double chi_one = atof(arg[2]);

  // convert chi from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    chi[i] = chi_one/180.0 * PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void ImproperHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&chi[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void ImproperHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&chi[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&k[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&chi[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}
