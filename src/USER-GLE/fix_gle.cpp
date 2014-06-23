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
   Contributing author: Michele Ceriotti (Oxford), Joe Morrone (Columbia)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_gle.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL,ATOM};
//#define GLE_DEBUG 1

/* syntax for fix_gle:
 * fix_gle nfix id-group gle  ns  amatrix [temp|cmatrix] seed
 *
 *                                                                        */

/* ---------------------------------------------------------------------- */

namespace GLE {
  /*!!MC!! here begins the GLE stuff*/
#define midx(n,i,j) ((i)*(n)+(j))

//"stabilized" cholesky decomposition. does a LDL^t decomposition, then sets to zero the negative diagonal elements and gets MM^t
void StabCholesky(unsigned long n, const double* MMt, double* M)
{
    double L[n*n], D[n];

    unsigned long i,j,k;
    for(i=0; i<n; ++i) D[i]=0.;
    for(i=0; i<n*n; ++i) L[i]=0.;

    for(i=0; i<n; ++i)
    {
        L[midx(n,i,i)]=1.;
        for (j=0; j<i; j++)
        {
            printf("%d %d\n", i,j);
            L[midx(n,i,j)]=MMt[midx(n,i,j)];
            for (k=0; k<j; ++k) L[midx(n,i,j)]-=L[midx(n,i,k)]*L[midx(n,j,k)]*D[k];
            if (D[j]!=0.) L[midx(n,i,j)]/=D[j];
            else L[midx(n,i,j)]=0.0;
        }
        D[i]=MMt[midx(n,i,i)];
        for (k=0; k<i; ++k) D[i]-=L[midx(n,i,k)]*L[midx(n,i,k)]*D[k];
    }

    for(i=0; i<n; ++i) D[i]=(D[i]>0.?sqrt(D[i]):0.);

    for(i=0; i<n; ++i) for (j=0; j<n; j++) M[midx(n,i,j)]=L[midx(n,i,j)]*D[j];
}

void MyMult(unsigned long n, unsigned long m, unsigned long r, const double* A, const double* B, double* C, double cf=0.0)
{
   // !! TODO !! should probably call BLAS (or check if some efficient matrix-matrix multiply is implemented somewhere in LAMMPS)
   // (rows x cols)  :: A is n x r, B is r x m, C is n x m
   unsigned long i,j,k; double *cij;
   for (i=0; i<n; ++i)
      for (j=0; j<m; ++j)
      {
         cij=&C[midx(m,i,j)]; *cij *= cf;
         for (k=0; k<r; ++k) *cij+=A[midx(r,i,k)]*B[midx(m,k,j)];
      }
}

void MyTrans(unsigned long n, const double* A, double* AT)
{
   for (unsigned long i=0; i<n; ++i) for (unsigned long j=0; j<n; ++j) AT[j*n+i]=A[i*n+j];
}

void MyPrint(unsigned long n, const double* A)
{
   for (unsigned long k=0; k<n*n; ++k) { printf("%10.5e ", A[k]); if ((k+1)%n==0) printf("\n");}
}

//matrix exponential by scaling and squaring.
void MatrixExp(unsigned long n, const double* M, double* EM, unsigned long j=8, unsigned long k=8)
{
   double tc[j+1], SM[n*n], TMP[n*n];
   double onetotwok=pow(0.5,1.0*k);


   tc[0]=1;
   for (int i=0; i<j; ++i) tc[i+1]=tc[i]/(i+1.0);

   for (int i=0; i<n*n; ++i) { SM[i]=M[i]*onetotwok; EM[i]=0.0; TMP[i]=0.0; }

   for (int i=0; i<n; ++i) EM[midx(n,i,i)]=tc[j];

   //taylor exp of scaled matrix
   for (int p=j-1; p>=0; p--)
   {
      MyMult(n, n, n, SM, EM, TMP); for (int i=0; i<n*n; ++i) EM[i]=TMP[i];
      for (int i=0; i<n; ++i) EM[midx(n,i,i)]+=tc[p];
   }

   for (int p=0; p<k; p++)
   {  MyMult(n, n, n, EM, EM, TMP); for (int i=0; i<n*n; ++i) EM[i]=TMP[i]; }
}
}

FixGLE::FixGLE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix gle command");

  ns = atof(arg[3]);

  if (ns>atom->ns) error->all(FLERR, "Atom kind has not enough GLE slots");


  A = new double[(ns+1)*(ns+1)];
  C = new double[(ns+1)*(ns+1)];
  T=S=gle_tmp=gle_rnd=NULL; gle_buff=0;

  // LOADS A matrix
  printf("Reading A from %s\n", arg[4]);
  FILE* fgle = fopen(arg[4], "r");

  if (fgle==NULL) error->all(FLERR, "Cannot open A matrix GLE input");
  for (int i=0; i<(ns+1)*(ns+1); ++i)  //ASSUMES A IS IN THE CORRECT INVERSE TIME UNITS
    if (!fscanf(fgle,"%lf",&(A[i]))) error->all(FLERR, "Cannot read A matrix GLE input");

  fclose(fgle);

  temp=atof(arg[5]);
  if (temp==0.0)
  {
    printf("Reading C from %s\n", arg[5]);
    fgle = fopen(arg[5], "r");
    if (fgle==NULL) error->all(FLERR, "Cannot open C matrix GLE input");
    for (int i=0; i<(ns+1)*(ns+1); ++i)  //ASSUMES C IS IN THE CORRECT TEMPERATURE UNITS
      if (!fscanf(fgle,"%lf",&(C[i]))) error->all(FLERR, "Cannot read C matrix GLE input");
  }
  else
  {
    for (int i=0; i<(ns+1)*(ns+1); ++i) C[i]=0.0;
    for (int i=0; i<(ns+1)*(ns+1); i+=(ns+2)) C[i]=temp*force->boltz/force->mvv2e;
  }

#ifdef GLE_DEBUG
  printf("A Matrix\n");
  GLE::MyPrint(ns+1,A);
  printf("C Matrix\n");
  GLE::MyPrint(ns+1,C);
#endif

  tau=10;
  int seed = atoi(arg[6]);

  // initialize Marsaglia RNG with processor-unique seed
  // NB: this means runs will not be the same with different numbers of processors
  random = new RanMars(lmp,seed + comm->me);

  if (seed <= 0) error->all(FLERR,"Illegal fix gle command");

  // allocate per-type arrays for mass-scaling
  sqrt_m = new double[atom->ntypes+1];

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixGLE::~FixGLE()
{
  delete random;
  delete [] sqrt_m;
  delete [] A; delete [] C; delete [] S; delete [] T;
  delete [] gle_tmp; delete [] gle_rnd;
}

/* ---------------------------------------------------------------------- */

int FixGLE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  mask |= THERMO_ENERGY;


  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGLE::init()
{

  dogle = 1;
  dorattle = 1;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // compute Langevin terms
  T = new double[(ns+1)*(ns+1)]; S = new double[(ns+1)*(ns+1)];
  double *tmp1 = new double[(ns+1)*(ns+1)], *tmp2 = new double[(ns+1)*(ns+1)];
  for (int i=0; i<(ns+1)*(ns+1); ++i) { tmp1[i]=-A[i]*update->dt*0.5; tmp2[i]=S[i]=0.0; }
  GLE::MatrixExp(ns+1,tmp1,T);

  GLE::MyMult(ns+1,ns+1,ns+1,T,C,tmp1);
  GLE::MyTrans(ns+1,T,tmp2);
  GLE::MyMult(ns+1,ns+1,ns+1,tmp1,tmp2,S);
  for (int i=0; i<(ns+1)*(ns+1); ++i) tmp1[i]=C[i]-S[i];

  GLE::StabCholesky(ns+1, tmp1, S);   //!TODO use symmetric square root, which is more stable.

#ifdef GLE_DEBUG
  printf("T Matrix\n");
  GLE::MyPrint(ns+1,T);
  printf("S Matrix\n");
  GLE::MyPrint(ns+1,S);
#endif

  delete[] tmp1;  delete[] tmp2;

  // set force prefactors
  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      sqrt_m[i] = sqrt(atom->mass[i]);
    }
  }

  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }
}

/* ---------------------------------------------------------------------- */

void FixGLE::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

void FixGLE::gle_integrate()
{
  double **v = atom->v;
  double **s = atom->s;
  double *rmass = atom->rmass, smi, ismi;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;



  if (nlocal > gle_buff)
  {
    if (gle_tmp !=NULL) {
      delete [] gle_tmp; delete [] gle_rnd;
    }
    gle_tmp = new double[nlocal*(ns+1)];
    gle_rnd = new double[nlocal*(ns+1)];
    gle_buff = nlocal;
  }

#ifdef GLE_DEBUG
  printf("!MC! GLE THERMO STEP dt=%f\n",update->dt);
#endif

  // gets kinetic energy before doing anything to the velocities
  double deltae=0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=rmass[i];
      for (int j = 0; j<3; ++j) { deltae-=ismi*v[i][j]*v[i][j]; }
    } }
  } else {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=atom->mass[type[i]];
      for (int j = 0; j<3; ++j) { deltae-=ismi*v[i][j]*v[i][j]; }
    } }
  }

  // s is just taken to be the mass-scaled velocity!
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        smi=sqrt(rmass[i]);
        for (int j = 0; j<3; ++j) { s[i][j] = v[i][j]*smi; }
      } }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        smi=sqrt_m[type[i]];
        for (int j = 0; j<3; ++j) { s[i][j] = v[i][j]*smi; }
      } }
  }


  for (int k = 0; k<3; k++) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) for (int j = 0; j<(ns+1); ++j)  gle_rnd[j*(nlocal)+i] = s[i][3*j+k];
    }
    GLE::MyMult(ns+1,nlocal,ns+1,T,gle_rnd,gle_tmp);
    for (int i = 0; i < nlocal*(ns+1); i++) gle_rnd[i] = random->gaussian();
    GLE::MyMult(ns+1,nlocal,ns+1,S,gle_rnd,gle_tmp,1.0);
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) for (int j = 0; j<(ns+1); ++j)  s[i][3*j+k] = gle_tmp[j*(nlocal)+i];
    }
  }

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=1.0/sqrt(rmass[i]);
      for (int j = 0; j<3; ++j) { v[i][j] = s[i][j]*ismi; }
    } }
  } else {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=1.0/sqrt_m[type[i]];
      for (int j = 0; j<3; ++j) { v[i][j] = s[i][j]*ismi; }
    } }
  }

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=rmass[i];
      for (int j = 0; j<3; ++j) { deltae+=ismi*v[i][j]*v[i][j]; }
    } }
  } else {
    for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ismi=atom->mass[type[i]];
      for (int j = 0; j<3; ++j) { deltae+=ismi*v[i][j]*v[i][j]; }
    } }
  }

  energy += deltae*0.5*force->mvv2e;
}

void FixGLE::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


  if (dogle) gle_integrate();

#ifdef GLE_DEBUG
  printf("!MC! GLE P1 STEP dt=%f\n",dtv);
  printf("!MC! GLE Q STEP\n");
#endif

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

void FixGLE::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

#ifdef GLE_DEBUG
  printf("!MC! GLE P2 STEP dt=%f\n",dtv);
#endif
  if (dogle) gle_integrate();
}
/* ---------------------------------------------------------------------- */

void FixGLE::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel==nlevels_respa-1) gle_integrate();
  dogle=0;
  if (ilevel == 0) initial_integrate(vflag);
  else { dorattle=1; final_integrate();} // If RATTLE S2 should not be called here, dorattle=0
}

void FixGLE::final_integrate_respa(int ilevel, int iloop)
{
  dtv = step_respa[ilevel]; //!MC!
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dogle=0;
  if(0) {  // If true turn on stage 2 rattle for final integrate only
    if(ilevel==nlevels_respa-1)
      dorattle=1;
    else
      dorattle=0;
  } else { dorattle=1; }
  final_integrate();
  if (ilevel==nlevels_respa-1) gle_integrate();
}


double FixGLE::compute_scalar()
{
  double energy_me = energy;
  double energy_all;
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);

  return -energy_all;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixGLE::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"t_target") == 0) {
    return &temp;
  }
  return NULL;
}

