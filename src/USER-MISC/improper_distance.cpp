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
   Contributing author: Paolo Raiteri (Curtin University)
------------------------------------------------------------------------- */

#include "improper_distance.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperDistance::ImproperDistance(LAMMPS *lmp) : Improper(lmp) {}

/* ---------------------------------------------------------------------- */

ImproperDistance::~ImproperDistance()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(chi);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperDistance::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double xab, yab, zab; // bond 1-2
  double xac, yac, zac; // bond 1-3
  double xad, yad, zad; // bond 1-4
  double xbc, ybc, zbc; // bond 2-3
  double xbd, ybd, zbd; // bond 2-4
  double xna, yna, zna, rna; // normal
  double da;

  double eimproper,f1[3],f2[3],f3[3],f4[3];
//  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
//  double s12,c,s,domega,a,a11,a22,a33,a12,a13,a23;
  double domega,a;

  eimproper = 0.0;
  ev_init(eflag,vflag);

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

    // geometry of 4-body
    // 1 is the central atom
    // 2-3-4 are ment to be equivalent
    // I need the bonds between 2-3 and 2-4 to get the plane normal
    // Then I need the bond 1-2 to project it onto the normal to the plane

    // bond 1->2
    xab = x[i2][0] - x[i1][0];
    yab = x[i2][1] - x[i1][1];
    zab = x[i2][2] - x[i1][2];
    domain->minimum_image(xab,yab,zab);

    // bond 1->3
    xac = x[i3][0] - x[i1][0];
    yac = x[i3][1] - x[i1][1];
    zac = x[i3][2] - x[i1][2];
    domain->minimum_image(xac,yac,zac);

    // bond 1->4
    xad = x[i4][0] - x[i1][0];
    yad = x[i4][1] - x[i1][1];
    zad = x[i4][2] - x[i1][2];
    domain->minimum_image(xad,yad,zad);

    // bond 2-3
    xbc = x[i3][0] - x[i2][0];
    ybc = x[i3][1] - x[i2][1];
    zbc = x[i3][2] - x[i2][2];
    domain->minimum_image(xbc,ybc,zbc);

    // bond 2-4
    xbd = x[i4][0] - x[i2][0];
    ybd = x[i4][1] - x[i2][1];
    zbd = x[i4][2] - x[i2][2];
    domain->minimum_image(xbd,ybd,zbd);

    xna =   ybc*zbd - zbc*ybd;
    yna = -(xbc*zbd - zbc*xbd);
    zna =   xbc*ybd - ybc*xbd;
    rna = 1.0 / sqrt(xna*xna+yna*yna+zna*zna);
    xna *= rna;
    yna *= rna;
    zna *= rna;

    da = xna*xab + yna*yab + zna*zab;

    domega = k[type]*da*da + chi[type]*da*da*da*da;
    //printf("%3i %3i %3i %3i %10.5f %10.5f \n",i1,i2,i3,i4,da,domega);
    a =  2.0* (k[type]*da + 2.0*chi[type]*da*da*da);

    if (eflag) eimproper = domega;

    f1[0] = a*( xna);
    f1[1] = a*( yna);
    f1[2] = a*( zna);

    f2[0] = a*( -xna               -yab*(zbd-zbc)*rna +zab*(ybd-ybc)*rna -da*( -yna*(zbd-zbc) + zna*(ybd-ybc) )*rna);
    f2[1] = a*( +xab*(zbd-zbc)*rna -yna               +zab*(xbc-xbd)*rna -da*( +xna*(zbd-zbc) + zna*(xbc-xbd) )*rna);
    f2[2] = a*( -xab*(ybd-ybc)*rna -yab*(xbc-xbd)*rna -zna               -da*( +xna*(ybc-ybd) - yna*(xbc-xbd) )*rna);

    f3[0] = a*( (           yab*zbd -zab*ybd ) *rna +da*( -yna*zbd +zna*ybd )*rna);
    f3[1] = a*( ( -xab*zbd          +zab*xbd ) *rna +da*( +xna*zbd -zna*xbd )*rna);
    f3[2] = a*( ( +xab*ybd -yab*xbd          ) *rna +da*( -xna*ybd +yna*xbd )*rna);

    f4[0] = a*( (          -yab*zbc +zab*ybc ) *rna -da*( -yna*zbc +zna*ybc )*rna);
    f4[1] = a*( ( +xab*zbc          -zab*xbc ) *rna -da*( +xna*zbc -zna*xbc )*rna);
    f4[2] = a*( ( -xab*ybc +yab*xbc          ) *rna -da*( -xna*ybc +yna*xbc )*rna);
    //printf("%10.5f %10.5f %10.5f \n",f1[0],f1[1],f1[2]);
    //printf("%10.5f %10.5f %10.5f \n",f2[0],f2[1],f2[2]);
    //printf("%10.5f %10.5f %10.5f \n",f3[0],f3[1],f3[2]);
    //printf("%10.5f %10.5f %10.5f \n",f4[0],f4[1],f4[2]);

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f2,f3,f4,
       xab,yab,zab,xac,yac,zac,xad-xac,yad-yac,zad-zac);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperDistance::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");
  memory->create(chi,n+1,"improper:chi");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperDistance::coeff(int narg, char **arg)
{
//  if (which > 0) return;
  if (narg != 3) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double chi_one = force->numeric(FLERR,arg[2]);

  // convert chi from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    chi[i] = chi_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperDistance::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&chi[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperDistance::read_restart(FILE *fp)
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

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperDistance::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],chi[i]);
}
