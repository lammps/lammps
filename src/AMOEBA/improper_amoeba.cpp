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

#include "improper_amoeba.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"

#include <cmath>

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
  if (disable) return;

  int ia,ib,ic,id,n,type;
  double xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid;
  double xab,yab,zab,xcb,ycb,zcb,xdb,ydb,zdb,xad,yad,zad,xcd,ycd,zcd;
  double rad2,rcd2,rdb2,dot,cc,ee;
  double sine,angle;
  double dt,dt2,dt3,dt4,e;
  double deddt,sign,dedcos,term;
  double dccdxia,dccdyia,dccdzia,dccdxic,dccdyic,dccdzic;
  double dccdxid,dccdyid,dccdzid;
  double deedxia,deedyia,deedzia,deedxic,deedyic,deedzic;
  double deedxid,deedyid,deedzid;
  double fa[3],fb[3],fc[3],fd[3];

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  // conversion factors for radians to degrees and vice versa

  double rad2degree = 180.0/MY_PI;
  double eprefactor = 1.0 / (rad2degree*rad2degree);
  double fprefactor = 1.0 / rad2degree;

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

    // angle needs to be in degrees for Tinker formulas
    // b/c opbend_3456 coeffs are in mixed units

    angle = rad2degree * asin(sine);
    dt = angle;
    dt2 = dt * dt;
    dt3 = dt2 * dt;
    dt4 = dt2 * dt2;
    e = eprefactor * k[type] * dt2 *
      (1.0 + opbend_cubic*dt + opbend_quartic*dt2 +
       opbend_pentic*dt3 + opbend_sextic*dt4);

    deddt = fprefactor * k[type] * dt *
      (2.0 + 3.0*opbend_cubic*dt + 4.0*opbend_quartic*dt2 +
       5.0*opbend_pentic*dt3 + 6.0*opbend_sextic*dt4);
    sign = (ee >= 0.0) ? 1.0 : -1.0;
    dedcos = -deddt * sign / sqrt(cc*rdb2 - ee*ee);

    // chain rule terms for first derivative components

    term = ee / cc;
    dccdxia = (xad*rcd2-xcd*dot) * term;
    dccdyia = (yad*rcd2-ycd*dot) * term;
    dccdzia = (zad*rcd2-zcd*dot) * term;
    dccdxic = (xcd*rad2-xad*dot) * term;
    dccdyic = (ycd*rad2-yad*dot) * term;
    dccdzic = (zcd*rad2-zad*dot) * term;
    dccdxid = -dccdxia - dccdxic;
    dccdyid = -dccdyia - dccdyic;
    dccdzid = -dccdzia - dccdzic;

    term = ee / rdb2;
    deedxia = ydb*zcb - zdb*ycb;
    deedyia = zdb*xcb - xdb*zcb;
    deedzia = xdb*ycb - ydb*xcb;
    deedxic = yab*zdb - zab*ydb;
    deedyic = zab*xdb - xab*zdb;
    deedzic = xab*ydb - yab*xdb;
    deedxid = ycb*zab - zcb*yab + xdb*term;
    deedyid = zcb*xab - xcb*zab + ydb*term;
    deedzid = xcb*yab - ycb*xab + zdb*term;

    //  compute first derivative components for this angle

    fa[0] = dedcos * (dccdxia+deedxia);
    fa[1] = dedcos * (dccdyia+deedyia);
    fa[2] = dedcos * (dccdzia+deedzia);
    fc[0] = dedcos * (dccdxic+deedxic);
    fc[1] = dedcos * (dccdyic+deedyic);
    fc[2] = dedcos * (dccdzic+deedzic);
    fd[0] = dedcos * (dccdxid+deedxid);
    fd[1] = dedcos * (dccdyid+deedyid);
    fd[2] = dedcos * (dccdzid+deedzid);
    fb[0] = -fa[0] - fc[0] - fd[0];
    fb[1] = -fa[1] - fc[1] - fd[1];
    fb[2] = -fa[2] - fc[2] - fd[2];

    // apply force to each of 4 atoms

    if (newton_bond || id < nlocal) {
      f[id][0] -= fd[0];
      f[id][1] -= fd[1];
      f[id][2] -= fd[2];
    }

    if (newton_bond || ib < nlocal) {
      f[ib][0] -= fb[0];
      f[ib][1] -= fb[1];
      f[ib][2] -= fb[2];
    }

    if (newton_bond || ia < nlocal) {
      f[ia][0] -= fa[0];
      f[ia][1] -= fa[1];
      f[ia][2] -= fa[2];
    }

    if (newton_bond || ic < nlocal) {
      f[ic][0] -= fc[0];
      f[ic][1] -= fc[1];
      f[ic][2] -= fc[2];
    }

    if (evflag) {
      fd[0] = -fd[0]; fd[1] = -fd[1]; fd[2] = -fd[2];
      fa[0] = -fa[0]; fa[1] = -fa[1]; fa[2] = -fa[2];
      fc[0] = -fc[0]; fc[1] = -fc[1]; fc[2] = -fc[2];
      ev_tally(id,ib,ia,ic,nlocal,newton_bond,e,fd,fa,fc,
               xdb,ydb,zdb,xab,yab,zab,xic-xia,yic-yia,zic-zia);
    }
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
   set opbend higher-order term weights from PairAmoeba
------------------------------------------------------------------------- */

void ImproperAmoeba::init_style()
{
  // check if PairAmoeba disabled improper terms

  Pair *pair = nullptr;
  pair = force->pair_match("^amoeba",0,0);
  if (!pair) pair = force->pair_match("^hippo",0,0);

  if (!pair) error->all(FLERR,"Improper amoeba could not find pair amoeba/hippo");

  int tmp;
  int flag = *((int *) pair->extract("improper_flag",tmp));
  disable = flag ? 0 : 1;

  // also extract opbend params

  int dim;
  opbend_cubic = *(double *) pair->extract("opbend_cubic",dim);
  opbend_quartic = *(double *) pair->extract("opbend_quartic",dim);
  opbend_pentic = *(double *) pair->extract("opbend_pentic",dim);
  opbend_sextic = *(double *) pair->extract("opbend_sextic",dim);
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

  if (comm->me == 0)
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nimpropertypes,fp,nullptr,error);
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
