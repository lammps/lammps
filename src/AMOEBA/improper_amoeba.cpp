// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "improper_amoeba.h"

#include <cmath>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
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

ImproperAmoeba::ImproperAmoeba(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperAmoeba::~ImproperAmoeba()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperAmoeba::compute(int eflag, int vflag)
{
  int ia,ib,ic,id,n,type;
  double xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid;
  double xab,yab,zab,xcb,ycb,zcb,xdb,ydb,zdb,xad,yad,zad,xcd,ycd,zcd;
  double rad2,rcd2,rdb2,dot,cc,ee;
  double sine,angle;
  double eimproper,f1[3],f2[3],f3[3],f4[3];

  eimproper = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {

    // in Tinker code, atom1 = D, atom2 = B, atom3 = A, atom4 = C
    // for Alligner angle:
    //   atoms A,C,D form a plane, B is out-of-plane
    //   angle is between plane and the vector from D to B

    id = improperlist[n][0];
    ib = improperlist[n][1];
    ia = improperlist[n][2];
    ic = improperlist[n][3];
    type = improperlist[n][4];

    // coordinates of the atoms at trigonal center

    xia = x[ia][0];
    yia = x[ia][1];
    zia = x[ia][2];
    xib = x[ib][0];
    yib = x[ib][1];
    zib = x[ib][2];
    xic = x[ic][0];
    yic = x[ic][1];
    zic = x[ic][2];
    xid = x[id][0];
    yid = x[id][1];
    zid = x[id][2];

    // compute the out-of-plane bending angle

    xab = xia - xib;
    yab = yia - yib;
    zab = zia - zib;
    xcb = xic - xib;
    ycb = yic - yib;
    zcb = zic - zib;
    xdb = xid - xib;
    ydb = yid - yib;
    zdb = zid - zib;
    xad = xia - xid;
    yad = yia - yid;
    zad = zia - zid;
    xcd = xic - xid;
    ycd = yic - yid;
    zcd = zic - zid;
    
    // Allinger angle between A-C-D plane and D-B vector for D-B < AC
      
    rad2 = xad*xad + yad*yad + zad*zad;
    rcd2 = xcd*xcd + ycd*ycd + zcd*zcd;
    dot = xad*xcd + yad*ycd + zad*zcd;
    cc = rad2*rcd2 - dot*dot;

    // find the out-of-plane angle bending energy

    ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb) + zdb*(xab*ycb-yab*xcb);
    rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
    if (rdb2 == 0.0 || cc == 0.0) continue;

    sine = fabs(ee) / sqrt(cc*rdb2);
    sine = MIN(1.0,sine);
    angle = asin(sine) * MY_PI/180.0;

    //dt = angle;
    //dt2 = dt * dt;
    //dt3 = dt2 * dt;
    //dt4 = dt2 * dt2;
    //e = opbunit * force * dt2 * (1.0d0+copb*dt+qopb*dt2+popb*dt3+sopb*dt4);
    //deddt = opbunit * force * dt * radian * 
    // (2.0d0 + 3.0d0*copb*dt + 4.0d0*qopb*dt2 + 5.0d0*popb*dt3 + 6.0d0*sopb*dt4);
    //dedcos = -deddt * sign(1.0d0,ee) / sqrt(cc*rdb2-ee*ee);

    /*
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
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
    */
  }
}

/* ---------------------------------------------------------------------- */

void ImproperAmoeba::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperAmoeba::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nimpropertypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);

  // convert chi from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperAmoeba::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperAmoeba::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nimpropertypes,fp,nullptr,error);
  }
  MPI_Bcast(&k[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperAmoeba::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g\n",i,k[i]);
}
