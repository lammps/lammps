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
   Contributing author: Hendrik Heenen (Technical University of Munich)
                        and Rochus Schmid (Ruhr-Universitaet Bochum)
   [ based on improper_fourier.cpp Loukas D. Peristeras (Scienomics SARL) ]
   [ based on improper_umbrella.cpp Tod A Pascal (Caltech) ]
   [ abbreviated from and verified via DLPOLY2.0 ]
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "improper_inversion_harmonic.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperInversionHarmonic::ImproperInversionHarmonic(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperInversionHarmonic::~ImproperInversionHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kw);
    memory->destroy(w0);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperInversionHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double rrvb1,rrvb2,rrvb3,rr2vb1,rr2vb2,rr2vb3;

  ev_init(eflag,vflag);

  double **x = atom->x;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    // 1st bond - IJ

    vb1x = x[i2][0] - x[i1][0];
    vb1y = x[i2][1] - x[i1][1];
    vb1z = x[i2][2] - x[i1][2];
    rrvb1 = 1.0/sqrt(vb1x*vb1x+vb1y*vb1y+vb1z*vb1z);
    rr2vb1 = rrvb1*rrvb1;

    // 2nd bond - IK

    vb2x = x[i3][0] - x[i1][0];
    vb2y = x[i3][1] - x[i1][1];
    vb2z = x[i3][2] - x[i1][2];
    rrvb2 = 1.0/sqrt(vb2x*vb2x+vb2y*vb2y+vb2z*vb2z);
    rr2vb2 = rrvb2*rrvb2;

    // 3rd bond - IL

    vb3x = x[i4][0] - x[i1][0];
    vb3y = x[i4][1] - x[i1][1];
    vb3z = x[i4][2] - x[i1][2];
    rrvb3 = 1.0/sqrt(vb3x*vb3x+vb3y*vb3y+vb3z*vb3z);
    rr2vb3 = rrvb3*rrvb3;

    // compute all three inversion angles
    invang(i1,i2,i3,i4, type,evflag,eflag,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1);
    invang(i1,i3,i4,i2, type,evflag,eflag,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2);
    invang(i1,i4,i2,i3, type,evflag,eflag,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3);
  }
}

/* ----------------------------------------------------------------------
   compute inversion angles + energy and forces
------------------------------------------------------------------------- */

void ImproperInversionHarmonic::invang(const int &i1,const int &i2,
          const int &i3,const int &i4,
          const int &type,const int &evflag,const int &eflag,
          const double &vb1x, const double &vb1y, const double &vb1z,
          const double &rrvb1, const double &rr2vb1,
          const double &vb2x, const double &vb2y, const double &vb2z,
          const double &rrvb2, const double &rr2vb2,
          const double &vb3x, const double &vb3y, const double &vb3z,
          const double &rrvb3, const double &rr2vb3)
{
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double omega,cosomega,domega,gomega,rjk,rjl;
  double upx,upy,upz,upn,rup,umx,umy,umz,umn,rum,wwr;
  double rucb,rudb,rvcb,rvdb,rupupn,rumumn;

  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  eimproper = 0.0;

  // scalar products of IJ*IK and IJ*IL
  rjk=vb3x*vb2x+vb3y*vb2y+vb3z*vb2z;
  rjl=vb1x*vb3x+vb1y*vb3y+vb1z*vb3z;

  // unit-vector: IK+IL
  upx=vb2x*rrvb2+vb1x*rrvb1;
  upy=vb2y*rrvb2+vb1y*rrvb1;
  upz=vb2z*rrvb2+vb1z*rrvb1;
  upn=1.0/sqrt(upx*upx+upy*upy+upz*upz);
  upx=upx*upn;
  upy=upy*upn;
  upz=upz*upn;
  rup=vb3x*upx+vb3y*upy+vb3z*upz;

  // unit-vector: IK-IL
  umx=vb2x*rrvb2-vb1x*rrvb1;
  umy=vb2y*rrvb2-vb1y*rrvb1;
  umz=vb2z*rrvb2-vb1z*rrvb1;
  umn=1.0/sqrt(umx*umx+umy*umy+umz*umz);
  umx=umx*umn;
  umy=umy*umn;
  umz=umz*umn;
  rum=vb3x*umx+vb3y*umy+vb3z*umz;

  // angle theta
  wwr=sqrt(rup*rup+rum*rum);

  cosomega=wwr*rrvb3;

  if (cosomega > 1.0) cosomega = 1.0;

  omega=acos(cosomega);

  domega = acos(cosomega)-w0[type];
  if (eflag) eimproper = kw[type]*(domega*domega);

  // kw[type] is divided by 3 -> threefold contribution
  gomega=0.0;
  if (omega*omega > 1.0e-24)  gomega=2.0*kw[type]*(domega)/(sin(omega));

  // projection IK and IL on unit vectors and contribution on IK and IL
  rucb = rjk-rup*(vb2x*upx+vb2y*upy+vb2z*upz);
  rudb = rjl-rup*(vb1x*upx+vb1y*upy+vb1z*upz);
  rvcb = rjk-rum*(vb2x*umx+vb2y*umy+vb2z*umz);
  rvdb = rjl-rum*(vb1x*umx+vb1y*umy+vb1z*umz);

  rupupn = rup*upn;
  rumumn = rum*umn;

  // force contributions of angle
  f2[0]=gomega*(-cosomega*vb3x*rr2vb3+rrvb3*(rup*upx+rum*umx)/wwr);
  f2[1]=gomega*(-cosomega*vb3y*rr2vb3+rrvb3*(rup*upy+rum*umy)/wwr);
  f2[2]=gomega*(-cosomega*vb3z*rr2vb3+rrvb3*(rup*upz+rum*umz)/wwr);

  f3[0]=gomega*rrvb3*(rupupn*rrvb2*(vb3x-rup*upx-rucb*vb2x*rr2vb2) +
        rumumn*rrvb2*(vb3x-rum*umx-rvcb*vb2x*rr2vb2))/wwr;
  f3[1]=gomega*rrvb3*(rupupn*rrvb2*(vb3y-rup*upy-rucb*vb2y*rr2vb2) +
        rumumn*rrvb2*(vb3y-rum*umy-rvcb*vb2y*rr2vb2))/wwr;
  f3[2]=gomega*rrvb3*(rupupn*rrvb2*(vb3z-rup*upz-rucb*vb2z*rr2vb2) +
        rumumn*rrvb2*(vb3z-rum*umz-rvcb*vb2z*rr2vb2))/wwr;

  f4[0]=gomega*rrvb3*(rupupn*rrvb1*(vb3x-rup*upx-rudb*vb1x*rr2vb1) -
        rumumn*rrvb1*(vb3x-rum*umx-rvdb*vb1x*rr2vb1))/wwr;
  f4[1]=gomega*rrvb3*(rupupn*rrvb1*(vb3y-rup*upy-rudb*vb1y*rr2vb1) -
        rumumn*rrvb1*(vb3y-rum*umy-rvdb*vb1y*rr2vb1))/wwr;
  f4[2]=gomega*rrvb3*(rupupn*rrvb1*(vb3z-rup*upz-rudb*vb1z*rr2vb1) -
        rumumn*rrvb1*(vb3z-rum*umz-rvdb*vb1z*rr2vb1))/wwr;

  f1[0] = -(f2[0] + f3[0] + f4[0]);
  f1[1] = -(f2[1] + f3[1] + f4[1]);
  f1[2] = -(f2[2] + f3[2] + f4[2]);

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

  if (evflag) {
    double rb3x, rb3y, rb3z;

    rb3x = vb1x - vb2x;
    rb3y = vb1y - vb2y;
    rb3z = vb1z - vb2z;

    ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f2,f3,f4,
             vb3x,vb3y,vb3z,
             vb2x,vb2y,vb2z,
             rb3x,rb3y,rb3z);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperInversionHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(kw,n+1,"improper:kw");
  memory->create(w0,n+1,"improper:w0");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperInversionHarmonic::coeff(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double w_one = force->numeric(FLERR,arg[2]);

  // convert w0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kw[i] = k_one/3.0; // parameter division due to 3 vector averaging
    w0[i] = w_one/180.0 * MY_PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperInversionHarmonic::write_restart(FILE *fp)
{
  fwrite(&kw[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&w0[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperInversionHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kw[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&w0[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&kw[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&w0[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperInversionHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,kw[i],w0[i]/MY_PI*180.0);
}
