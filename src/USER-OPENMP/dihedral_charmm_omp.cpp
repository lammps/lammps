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
#include "dihedral_charmm_omp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair_omp.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

DihedralCharmmOMP::DihedralCharmmOMP(LAMMPS *lmp) : DihedralOMP(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralCharmmOMP::~DihedralCharmmOMP()
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

void DihedralCharmmOMP::compute(int eflag, int vflag)
 {
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = 0;

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (weightflag && vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  if (evflag) {
    if (eflag) {
      if (force->newton_bond) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_bond) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_bond) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void DihedralCharmmOMP::eval()
{
 
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i1,i2,i3,i4,i,m,n,type,tid;
    double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
    double edihedral,f1[3],f2[3],f3[3],f4[3];
    double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
    double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
    double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
    double c,s,p,sx2,sy2,sz2;
    int itype,jtype;
    double delx,dely,delz,rsq,r2inv,r6inv;
    double forcecoul,forcelj,fpair,ecoul,evdwl;

    edihedral = evdwl = ecoul = 0.0;

    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;
    const int nthreads = comm->nthreads;
    const double qqrd2e = force->qqrd2e;

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *atomtype = atom->type;
    int **dihedrallist = neighbor->dihedrallist;
    int ndihedrallist = neighbor->ndihedrallist;

    // loop over list of dihedrals

    int nfrom, nto;
    f = loop_setup_thr(f,nfrom, nto, tid, ndihedrallist, nall, nthreads);
    for (n = nfrom; n < nto; ++n) {
      i1 = dihedrallist[n][0];
      i2 = dihedrallist[n][1];
      i3 = dihedrallist[n][2];
      i4 = dihedrallist[n][3];
      type = dihedrallist[n][4];

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
	  sprintf(str,"Dihedral problem: %d/%d " BIGINT_FORMAT " %d %d %d %d",
		  me,tid,update->ntimestep,
		  atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	  error->warning(str,0);
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

      if (EFLAG) edihedral = k[type] * p;
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

      df = -k[type] * df1;

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

      if (NEWTON_BOND || i1 < nlocal) {
	f[i1][0] += f1[0];
	f[i1][1] += f1[1];
	f[i1][2] += f1[2];
      }

      if (NEWTON_BOND || i2 < nlocal) {
	f[i2][0] += f2[0];
	f[i2][1] += f2[1];
	f[i2][2] += f2[2];
      }

      if (NEWTON_BOND || i3 < nlocal) {
	f[i3][0] += f3[0];
	f[i3][1] += f3[1];
	f[i3][2] += f3[2];
      }

      if (NEWTON_BOND || i4 < nlocal) {
	f[i4][0] += f4[0];
	f[i4][1] += f4[1];
	f[i4][2] += f4[2];
      }

      if (EVFLAG) ev_tally_thr(i1,i2,i3,i4,nlocal,NEWTON_BOND,edihedral,f1,f3,f4,
			       vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,tid);

      // 1-4 LJ and Coulomb interactions
      // tally energy/virial in pair, using newton_bond as newton flag

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
	fpair = weight[type] * (forcelj+forcecoul)*r2inv;

	if (EFLAG) {
	  ecoul = weight[type] * forcecoul;
	  evdwl = r6inv * (lj14_3[itype][jtype]*r6inv - lj14_4[itype][jtype]);
	  evdwl *= weight[type];
	}

	if (NEWTON_BOND || i1 < nlocal) {
	  f[i1][0] += delx*fpair;
	  f[i1][1] += dely*fpair;
	  f[i1][2] += delz*fpair;
	}
	if (NEWTON_BOND || i4 < nlocal) {
	  f[i4][0] -= delx*fpair;
	  f[i4][1] -= dely*fpair;
	  f[i4][2] -= delz*fpair;
	}
	// NOTE: we need thread safety here.
	if (EVFLAG) pair->ev_tally_thr(i1,i4,nlocal,NEWTON_BOND,evdwl,ecoul,
				       fpair,delx,dely,delz,tid);
      }
    }
    // reduce per thread forces into global force array.
    force_reduce_thr(atom->f, nall, nthreads, tid);
  }
  if (EVFLAG) ev_reduce_thr();
}

/* ---------------------------------------------------------------------- */

void DihedralCharmmOMP::allocate()
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

void DihedralCharmmOMP::coeff(int narg, char **arg)
{
  if (narg != 5) error->all("Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);
  
  // require integer values of shift for backwards compatibility
  // arbitrary phase angle shift could be allowed, but would break
  //   backwards compatibility and is probably not needed
  
  double k_one = force->numeric(arg[1]);
  int multiplicity_one = force->inumeric(arg[2]);
  int shift_one = force->inumeric(arg[3]);
  double weight_one = force->numeric(arg[4]);

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

void DihedralCharmmOMP::init_style()
{
  // insure use of CHARMM pair_style if any weight factors are non-zero
  // set local ptrs to LJ 14 arrays setup by Pair

  weightflag = 0;
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    if (weight[i] > 0.0) weightflag = 1;

  if (weightflag) {
    int itmp;
    if (force->pair == NULL)
      error->all("Dihedral charmm is incompatible with Pair style");
    lj14_1 = (double **) force->pair->extract("lj14_1",itmp);
    lj14_2 = (double **) force->pair->extract("lj14_2",itmp);
    lj14_3 = (double **) force->pair->extract("lj14_3",itmp);
    lj14_4 = (double **) force->pair->extract("lj14_4",itmp);
    int *ptr = (int *) force->pair->extract("implicit",itmp);
    if (!lj14_1 || !lj14_2 || !lj14_3 || !lj14_4 || !ptr)
      error->all("Dihedral charmm/omp is incompatible with Pair style");
    implicit = *ptr;

    int *omp = (int *) force->pair->extract("omp",itmp);
    if (!omp)
      error->all("Dihedral charmm/omp requires use of /omp Pair style variant");
    pair = static_cast<PairOMP *>(force->pair);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralCharmmOMP::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&multiplicity[1],sizeof(int),atom->ndihedraltypes,fp);
  fwrite(&shift[1],sizeof(int),atom->ndihedraltypes,fp);
  fwrite(&weight[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralCharmmOMP::read_restart(FILE *fp)
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

/* ---------------------------------------------------------------------- */

double DihedralCharmmOMP::memory_usage()
{
  const int n=atom->ntypes;

  double bytes = DihedralOMP::memory_usage();

  bytes +=  6*(n+1) * sizeof(double);
  bytes +=  1*(n+1) * sizeof(int);

  return bytes;
}

