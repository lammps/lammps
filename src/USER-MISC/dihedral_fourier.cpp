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
   Contributing author: Loukas D. Peristeras (Scienomics SARL)
   [ based on dihedral_charmm.cpp Paul Crozier (SNL) ]
------------------------------------------------------------------------- */

#include "dihedral_fourier.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05

/* ---------------------------------------------------------------------- */

DihedralFourier::DihedralFourier(LAMMPS *lmp) : Dihedral(lmp)
{
   writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralFourier::~DihedralFourier()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(nterms);

    for (int i=1; i<= atom->ndihedraltypes; i++) {
      if ( k[i] ) delete [] k[i];
      if ( multiplicity[i] ) delete [] multiplicity[i];
      if ( shift[i] ) delete [] shift[i];
      if ( cos_shift[i] ) delete [] cos_shift[i];
      if ( sin_shift[i] ) delete [] sin_shift[i];
    }
    delete [] k;
    delete [] multiplicity;
    delete [] shift;
    delete [] cos_shift;
    delete [] sin_shift;

  }
}

/* ---------------------------------------------------------------------- */

void DihedralFourier::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,i,j,m,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1_,ddf1_,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
  double c,s,p_,sx2,sy2,sz2;

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

    // force and energy
    // p = sum(i=1,nterms) k_i*(1+cos(n_i*phi-d_i)
    // dp = dp / dphi
    edihedral = 0.0;
    df = 0.0;
    for (j=0; j<nterms[type]; j++)
    {
      m = multiplicity[type][j];
      p_ = 1.0;
      ddf1_ = df1_ = 0.0;

      for (i = 0; i < m; i++) {
        ddf1_ = p_*c - df1_*s;
        df1_ = p_*s + df1_*c;
        p_ = ddf1_;
      }

      p_ = p_*cos_shift[type][j] + df1_*sin_shift[type][j];
      df1_ = df1_*cos_shift[type][j] - ddf1_*sin_shift[type][j];
      df1_ *= -m;
      p_ += 1.0;

      if (m == 0) {
        p_ = 1.0 + cos_shift[type][j];
        df1_ = 0.0;
      }

      if (eflag) edihedral += k[type][j] * p_;

      df += (-k[type][j] * df1_);
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

void DihedralFourier::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(nterms,n+1,"dihedral:nterms");
  k = new double * [n+1];
  multiplicity = new int * [n+1];
  shift = new double * [n+1];
  cos_shift = new double * [n+1];
  sin_shift = new double * [n+1];
  for (int i = 1; i <= n; i++) {
    k[i] = shift[i] = cos_shift[i] = sin_shift[i] = 0;
    multiplicity[i] = 0;
  }

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralFourier::coeff(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  // require integer values of shift for backwards compatibility
  // arbitrary phase angle shift could be allowed, but would break
  //   backwards compatibility and is probably not needed

  double k_one;
  int multiplicity_one;
  double shift_one;
  int nterms_one = force->inumeric(FLERR,arg[1]);

  if (nterms_one < 1)
    error->all(FLERR,"Incorrect number of terms arg for dihedral coefficients");

  if (2+3*nterms_one < narg)
    error->all(FLERR,"Incorrect number of arguments for dihedral coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    nterms[i] = nterms_one;
    k[i] = new double [nterms_one];
    multiplicity[i] = new int [nterms_one];
    shift[i] = new double [nterms_one];
    cos_shift[i] = new double [nterms_one];
    sin_shift[i] = new double [nterms_one];
    for (int j = 0; j<nterms_one; j++) {
      int offset = 1+3*j;
      k_one = force->numeric(FLERR,arg[offset+1]);
      multiplicity_one = force->inumeric(FLERR,arg[offset+2]);
      shift_one = force->numeric(FLERR,arg[offset+3]);
      k[i][j] = k_one;
      multiplicity[i][j] = multiplicity_one;
      shift[i][j] = shift_one;
      cos_shift[i][j] = cos(MY_PI*shift_one/180.0);
      sin_shift[i][j] = sin(MY_PI*shift_one/180.0);
    }
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralFourier::write_restart(FILE *fp)
{

  fwrite(&nterms[1],sizeof(int),atom->ndihedraltypes,fp);
  for(int i = 1; i <= atom->ndihedraltypes; i++) {
    fwrite(k[i],sizeof(double),nterms[i],fp);
    fwrite(multiplicity[i],sizeof(int),nterms[i],fp);
    fwrite(shift[i],sizeof(double),nterms[i],fp);
  }

}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralFourier::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0)
    utils::sfread(FLERR,&nterms[1],sizeof(int),atom->ndihedraltypes,fp,NULL,error);

  MPI_Bcast(&nterms[1],atom->ndihedraltypes,MPI_INT,0,world);

  // allocate
  for (int i=1; i<=atom->ndihedraltypes; i++) {
    k[i] = new double [nterms[i]];
    multiplicity[i] = new int [nterms[i]];
    shift[i] = new double [nterms[i]];
    cos_shift[i] = new double [nterms[i]];
    sin_shift[i] = new double [nterms[i]];
  }

  if (comm->me == 0) {
    for (int i=1; i<=atom->ndihedraltypes; i++) {
      utils::sfread(FLERR,k[i],sizeof(double),nterms[i],fp,NULL,error);
      utils::sfread(FLERR,multiplicity[i],sizeof(int),nterms[i],fp,NULL,error);
      utils::sfread(FLERR,shift[i],sizeof(double),nterms[i],fp,NULL,error);
    }
  }

  for (int i=1; i<=atom->ndihedraltypes; i++) {
    MPI_Bcast(k[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(multiplicity[i],nterms[i],MPI_INT,0,world);
    MPI_Bcast(shift[i],nterms[i],MPI_DOUBLE,0,world);
  }

  for (int i=1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
    for (int j = 0; j < nterms[i]; j++) {
      cos_shift[i][j] = cos(MY_PI*shift[i][j]/180.0);
      sin_shift[i][j] = sin(MY_PI*shift[i][j]/180.0);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralFourier::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
  {
    fprintf(fp,"%d %d",i,nterms[i]);
    for(int j = 0; j < nterms[i]; j++)
       fprintf(fp," %g %d %g",k[i][j],multiplicity[i][j],shift[i][j]);
    fprintf(fp,"\n");
  }
}

