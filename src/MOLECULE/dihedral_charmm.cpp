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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "dihedral_charmm.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

DihedralCharmm::DihedralCharmm(LAMMPS *lmp) : Dihedral(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralCharmm::~DihedralCharmm()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(k);
    memory->sfree(multiplicity);
    memory->sfree(shift);
    memory->sfree(cos_shift);
    memory->sfree(sin_shift);
    memory->sfree(weight);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCharmm::compute(int eflag, int vflag)
{
  int i,m,n,i1,i2,i3,i4,type,factor;
  double rfactor;
  double vb1x,vb1y,vb1z,vb2x,vb2y;
  double vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z;
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;  
  double c,s,p,sx1,sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12;
  int itype,jtype;
  double delx,dely,delz,rsq,r2inv,r6inv;
  double fforce,forcecoul,forcelj,phicoul,philj;

  energy = 0.0;
  eng_coul = eng_vdwl = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *atomtype = atom->type;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double qqrd2e = force->qqrd2e;

  for (n = 0; n < ndihedrallist; n++) {

    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    if (newton_bond) factor = 4;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
      if (i3 < nlocal) factor++;
      if (i4 < nlocal) factor++;
      }
    rfactor = 0.25 * factor;

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
    domain->minimum_image(vb1x,vb1y,vb1z);

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    domain->minimum_image(vb2x,vb2y,vb2z);

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;
    domain->minimum_image(vb2xm,vb2ym,vb2zm);

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    domain->minimum_image(vb3x,vb3y,vb3z);
    
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
	fprintf(screen,"Dihedral problem: %d %d %d %d %d %d\n",
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
         
    m = multiplicity[type];
    p = 1.0;
    df1 = 0.0;
    
    for (i = 0; i < m; i++) {
      ddf1 = p*c - df1*s;
      df1 = p*s + df1*c;
      p = ddf1;
    }

    p = p*cos_shift[type] + df1*sin_shift[type];
    df1 = df1*cos_shift[type] - ddf1*sin_shift[type];
    df1 *= -m;
    p += 1.0;
 
    if (m == 0) {
      p = 1.0 + cos_shift[type];
      df1 = 0.0;
    }

    if (eflag) energy += rfactor * k[type] * p; 
       
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
    
    df = k[type] * df1;
    
    sx1 = df*dtfx;
    sy1 = df*dtfy;
    sz1 = df*dtfz;
    sx2 = -df*dtgx;
    sy2 = -df*dtgy;
    sz2 = -df*dtgz;
    sx12 = df*dthx;
    sy12 = df*dthy;
    sz12 = df*dthz; 
    
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
      virial[0] -= rfactor * (vb1x*sx1 + vb2x*sx2 + vb3x*sx12);
      virial[1] -= rfactor * (vb1y*sy1 + vb2y*sy2 + vb3y*sy12);
      virial[2] -= rfactor * (vb1z*sz1 + vb2z*sz2 + vb3z*sz12);
      virial[3] -= rfactor * (vb1x*sy1 + vb2x*sy2 + vb3x*sy12);
      virial[4] -= rfactor * (vb1x*sz1 + vb2x*sz2 + vb3x*sz12);
      virial[5] -= rfactor * (vb1y*sz1 + vb2y*sz2 + vb3y*sz12);
    }

    // 1-4 LJ and Coulomb interactions
    // force, energy, and virial

    if (weight[type] > 0.0) {

      itype = atomtype[i1];
      jtype = atomtype[i4];

      delx = x[i1][0] - x[i4][0];
      dely = x[i1][1] - x[i4][1];
      delz = x[i1][2] - x[i4][2];
      domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;

      if (implicit) forcecoul = qqrd2e * q[i1]*q[i4]*r2inv;
      else forcecoul = qqrd2e * q[i1]*q[i4]*sqrt(r2inv);
      forcelj = r6inv * (lj14_1[itype][jtype]*r6inv - lj14_2[itype][jtype]);
      fforce = weight[type] * (forcelj+forcecoul)*r2inv;

      if (eflag) {
	phicoul = weight[type] * rfactor * forcecoul;
	philj = r6inv * (lj14_3[itype][jtype]*r6inv - lj14_4[itype][jtype]);
	philj = weight[type] * rfactor * philj;
	eng_coul += phicoul;
	eng_vdwl += philj;
      }

      if (newton_bond || i1 < nlocal) {
	f[i1][0] += delx*fforce;
	f[i1][1] += dely*fforce;
	f[i1][2] += delz*fforce;
      }
      if (newton_bond || i4 < nlocal) {
	f[i4][0] -= delx*fforce;
	f[i4][1] -= dely*fforce;
	f[i4][2] -= delz*fforce;
      }

      if (vflag) {
	virial[0] += rfactor * delx*delx*fforce;
	virial[1] += rfactor * dely*dely*fforce;
	virial[2] += rfactor * delz*delz*fforce;
	virial[3] += rfactor * delx*dely*fforce;
	virial[4] += rfactor * delx*delz*fforce;
	virial[5] += rfactor * dely*delz*fforce;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCharmm::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  k = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:k");
  multiplicity = (int *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:multiplicity");
  shift = (int *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:shift");
  cos_shift = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:cos_shift");
  sin_shift = (double *) 
    memory->smalloc((n+1)*sizeof(double),"dihedral:sin_shift");
  weight = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:weight");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralCharmm::coeff(int which, int narg, char **arg)
{
  if (which > 0) return;
  if (narg != 5) error->all("Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);
  
  // require integer values of shift for backwards compatibility
  // arbitrary phase angle shift could be allowed, but would break
  //   backwards compatibility and is probably not needed
  
  double k_one = atof(arg[1]);
  int multiplicity_one = atoi(arg[2]);
  int shift_one = atoi(arg[3]);
  double weight_one = atof(arg[4]);

  if (multiplicity_one < 0)
    error->all("Incorrect multiplicity arg for dihedral coefficients");
  if (weight_one < 0.0 || weight_one > 1.0) 
    error->all("Incorrect weight arg for dihedral coefficients");

  double PI = 4.0*atan(1.0);
                       
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    shift[i] = shift_one;
    cos_shift[i] = cos(PI*shift_one/180.0);
    sin_shift[i] = sin(PI*shift_one/180.0);
    multiplicity[i] = multiplicity_one;
    weight[i] = weight_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   error check and initialize all values needed for force computation 
------------------------------------------------------------------------- */

void DihedralCharmm::init_style()
{
  // insure use of CHARMM pair_style if any weight factors are non-zero
  // set local ptrs to LJ 14 arrays setup by Pair

  int weightflag = 0;
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    if (weight[i] > 0.0) weightflag = 1;

  if (weightflag) {
    if (force->pair == NULL)
      error->all("Dihedral charmm is incompatible with Pair style");
    lj14_1 = (double **) force->pair->extract("lj14_1");
    lj14_2 = (double **) force->pair->extract("lj14_2");
    lj14_3 = (double **) force->pair->extract("lj14_3");
    lj14_4 = (double **) force->pair->extract("lj14_4");
    int *ptr = (int *) force->pair->extract("implicit");
    if (!lj14_1 || !lj14_2 || !lj14_3 || !lj14_4 || !ptr)
      error->all("Dihedral charmm is incompatible with Pair style");
    implicit = *ptr;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralCharmm::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&multiplicity[1],sizeof(int),atom->ndihedraltypes,fp);
  fwrite(&shift[1],sizeof(int),atom->ndihedraltypes,fp);
  fwrite(&weight[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralCharmm::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&multiplicity[1],sizeof(int),atom->ndihedraltypes,fp);
    fread(&shift[1],sizeof(int),atom->ndihedraltypes,fp);
    fread(&weight[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&multiplicity[1],atom->ndihedraltypes,MPI_INT,0,world);
  MPI_Bcast(&shift[1],atom->ndihedraltypes,MPI_INT,0,world);
  MPI_Bcast(&weight[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  double PI = 4.0*atan(1.0);
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
    cos_shift[i] = cos(PI*shift[i]/180.0);
    sin_shift[i] = sin(PI*shift[i]/180.0);
  }
}
