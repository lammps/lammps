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
   Contributing author: Carsten Svaneborg, science@zqex.dk
------------------------------------------------------------------------- */

#include "dihedral_cosine_shift_exp.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

DihedralCosineShiftExp::DihedralCosineShiftExp(LAMMPS *lmp) : Dihedral(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralCosineShiftExp::~DihedralCosineShiftExp()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(umin);
    memory->destroy(a);
    memory->destroy(opt1);
    memory->destroy(cost);
    memory->destroy(sint);
    memory->destroy(theta);
    memory->destroy(doExpansion);
   }
}

/* ---------------------------------------------------------------------- */

void DihedralCosineShiftExp::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
  double c,s,sx2,sy2,sz2;
  double cccpsss,cssmscc,exp2;

  edihedral = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    // c,s calculation

    ax = vb1y*vb2zm - vb1z*vb2ym;
    ay = vb1z*vb2xm - vb1x*vb2zm;
    az = vb1x*vb2ym - vb1y*vb2xm;
    bx = vb3y*vb2zm - vb3z*vb2ym;
    by = vb3z*vb2xm - vb3x*vb2zm;
    bz = vb3x*vb2ym - vb3y*vb2xm;

    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
    rg = sqrt(rgsq);

    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
        char str[128];
        sprintf(str,"Dihedral problem: %d " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
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

    double aa=a[type];
    double uumin=umin[type];

    cccpsss = c*cost[type]+s*sint[type];
    cssmscc = c*sint[type]-s*cost[type];

  //  eflag=1;

    if (doExpansion[type])
       {  //  |a|<0.001 so use expansions relative precision <1e-5
            if (eflag) edihedral = -0.125*(1+cccpsss)*(4+aa*(cccpsss-1))*uumin;
            df=0.5*uumin*( cssmscc + 0.5*aa*cccpsss);
       }
     else
       {
            exp2=exp(0.5*aa*(1+cccpsss));
            if (eflag) edihedral = opt1[type]*(1-exp2);
            df= 0.5*opt1[type]*aa* ( exp2*cssmscc );
       }

    fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
    hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
    fga = fg*ra2inv*rginv;
    hgb = hg*rb2inv*rginv;
    gaa = -ra2inv*rg;
    gbb = rb2inv*rg;

    dtfx = gaa*ax;
    dtfy = gaa*ay;
    dtfz = gaa*az;
    dtgx = fga*ax - hgb*bx;
    dtgy = fga*ay - hgb*by;
    dtgz = fga*az - hgb*bz;
    dthx = gbb*bx;
    dthy = gbb*by;
    dthz = gbb*bz;

    sx2 = df*dtgx;
    sy2 = df*dtgy;
    sz2 = df*dtgz;

    f1[0] = df*dtfx;
    f1[1] = df*dtfy;
    f1[2] = df*dtfz;

    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];

    f4[0] = df*dthx;
    f4[1] = df*dthy;
    f4[2] = df*dthz;

    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];

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
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCosineShiftExp::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(doExpansion, n+1,  "dihedral:doExpansion");
  memory->create(umin,n+1,"dihedral:umin");
  memory->create(a,n+1,"dihedral:a");
  memory->create(sint,n+1,"dihedral:sind");
  memory->create(cost,n+1,"dihedral:cosd");
  memory->create(opt1,n+1,"dihedral:opt1");
  memory->create(theta,n+1,"dihedral:opt1");
  memory->create(setflag, n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralCosineShiftExp::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  double umin_   = force->numeric(FLERR,arg[1]);
  double theta0_ = force->numeric(FLERR,arg[2]);
  double a_      = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    doExpansion[i]=(fabs(a_)<0.001);
    umin[i]  = umin_;
    a[i]     = a_;
    cost[i]  = cos(theta0_*MathConst::MY_PI/180.0);
    sint[i]  = sin(theta0_*MathConst::MY_PI/180.0);
    theta[i] = theta0_*MathConst::MY_PI/180.0;

    if (!doExpansion[i]) opt1[i]=umin_/(exp(a_)-1);

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralCosineShiftExp::write_restart(FILE *fp)
{
  fwrite(&umin[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&cost[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&sint[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&theta[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralCosineShiftExp::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&umin[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&a[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&cost[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&sint[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&theta[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&umin[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&cost[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sint[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
    doExpansion[i]=(fabs(a[i])<0.01);
    if (!doExpansion[i]) opt1[i]=umin[i]/(exp(a[i])-1);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralCosineShiftExp::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,umin[i],
            theta[i]*180.0/MathConst::MY_PI,a[i]);
}
