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

/* ----------------------------------------------------------------------
   Contributing authors: Michele Ceriotti (EPFL), Joe Morrone (Stony Brook),
                         Axel Kohylmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_gle.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL,ATOM};

//#define GLE_DEBUG 1

#define MAXLINE 1024

/* syntax for fix_gle:
 * fix nfix id-group gle ns Tstart Tstop seed amatrix [noneq cmatrix] [every nmts]
 *
 *                                                                        */

/* ---------------------------------------------------------------------- */

namespace GLE {
/* GLE Utility functions -- might perhaps use standard LAMMPS routines if available */
#define midx(n,i,j) ((i)*(n)+(j))

//"stabilized" cholesky decomposition. does a LDL^t decomposition, then sets to zero the negative diagonal elements and gets MM^t
void StabCholesky(int n, const double* MMt, double* M)
{
  double *L = new double[n*n];
  double *D = new double[n];

  int i,j,k;
  for (i=0; i<n; ++i) D[i]=0.0;
  for (i=0; i<n*n; ++i) L[i]=0.0;

  for (i=0; i<n; ++i) {
    L[midx(n,i,i)]=1.0;
    for (j=0; j<i; j++) {
      L[midx(n,i,j)]=MMt[midx(n,i,j)];
      for (k=0; k<j; ++k) L[midx(n,i,j)]-=L[midx(n,i,k)]*L[midx(n,j,k)]*D[k];
      if (D[j]!=0.) L[midx(n,i,j)]/=D[j];
      else L[midx(n,i,j)]=0.0;
    }
    D[i]=MMt[midx(n,i,i)];
    for (k=0; k<i; ++k) D[i]-=L[midx(n,i,k)]*L[midx(n,i,k)]*D[k];
  }

  for (i=0; i<n; ++i) {
#ifdef GLE_DEBUG
    if (D[i]<0) fprintf(stderr,"GLE Cholesky: Negative diagonal term %le, has been set to zero.\n", D[i]);
#endif
    D[i]=(D[i]>0.0) ? sqrt(D[i]):0.0;
  }

  for (i=0; i<n; ++i)
    for (j=0; j<n; j++) M[midx(n,i,j)]=L[midx(n,i,j)]*D[j];

  delete[] D;
  delete[] L;
}

void MyMult(int n, int m, int r, const double* A, const double* B, double* C, double cf=0.0)
{
   // !! TODO !! should probably call BLAS (or check if some efficient matrix-matrix multiply is implemented somewhere in LAMMPS)
   // (rows x cols)  :: A is n x r, B is r x m, C is n x m
   int i,j,k; double *cij;
   for (i=0; i<n; ++i)
      for (j=0; j<m; ++j)
      {
         cij=&C[midx(m,i,j)]; *cij *= cf;
         for (k=0; k<r; ++k) *cij+=A[midx(r,i,k)]*B[midx(m,k,j)];
      }
}

#define MIN(A,B) ((A) < (B) ? (A) : (B))
  // BLAS-like version of MyMult(). AK 2014-08-06

inline void AkMult(const int n, const int m, const int r,
            const double * const A, const double * const B,
            double * const C, const double cf=0.0)
{
  // block buffer
  const int blk=64;
  double buf[blk*blk];
  int i,j,k,ib,jb,kb;

  for (i = 0; i < (n*m); ++i) C[i] *= cf;

  for (kb = 0; kb < r; kb +=blk) {
    // fill buffer
    const int ke = MIN(kb+blk,r);
    for (ib = 0; ib < n; ib += blk) {
      const int ie = MIN(ib+blk,n);
      for (i = ib; i < ie; ++i) {
        for (k = kb; k < ke; ++k) {
          buf[midx(blk,k-kb,i-ib)] = A[midx(r,i,k)];
        }
      }

      for (jb = 0; jb < m; jb += blk) {
        double tmp;
        const int je = MIN(jb+blk,m);

        for (j = jb; j < je; ++j) {
          for (i = ib; i < ie; ++i) {
            tmp = 0.0;
            for (k = kb; k < ke; ++k)
              tmp += buf[midx(blk,k-kb,i-ib)] * B[midx(m,k,j)];

            C[midx(m,i,j)] += tmp;
          }
        }
      }
    }
  }
}

void MyTrans(int n, const double* A, double* AT)
{
   for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) AT[j*n+i]=A[i*n+j];
}

void MyPrint(int n, const double* A)
{
   for (int k=0; k<n*n; ++k) { printf("%10.5e ", A[k]); if ((k+1)%n==0) printf("\n");}
}

//matrix exponential by scaling and squaring.
void MatrixExp(int n, const double* M, double* EM, int j=8, int k=8)
{
  double *tc = new double[j+1];
  double *SM = new double[n*n];
  double *TMP = new double[n*n];
  double onetotwok=pow(0.5,1.0*k);


  tc[0]=1;
  for (int i=0; i<j; ++i) tc[i+1]=tc[i]/(i+1.0);

  for (int i=0; i<n*n; ++i) { SM[i]=M[i]*onetotwok; EM[i]=0.0; TMP[i]=0.0; }

  for (int i=0; i<n; ++i) EM[midx(n,i,i)]=tc[j];

  //taylor exp of scaled matrix
  for (int p=j-1; p>=0; p--) {
    MyMult(n, n, n, SM, EM, TMP); for (int i=0; i<n*n; ++i) EM[i]=TMP[i];
    for (int i=0; i<n; ++i) EM[midx(n,i,i)]+=tc[p];
  }

  for (int p=0; p<k; p++) {
     MyMult(n, n, n, EM, EM, TMP);
     for (int i=0; i<n*n; ++i) EM[i]=TMP[i];
  }
  delete[] tc;
  delete[] SM;
  delete[] TMP;
}
}

/* ---------------------------------------------------------------------- */

FixGLE::FixGLE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR,"Illegal fix gle command. Expecting: fix <fix-ID>"
               " <group-ID> gle <ns> <Tstart> <Tstop> <seed> <Amatrix>");

  ecouple_flag = 1;
  restart_peratom = 1;
  time_integrate = 1;

  // number of additional momenta
  ns = utils::inumeric(FLERR,arg[3],false,lmp);
  ns1sq = (ns+1)*(ns+1);

  // allocate GLE matrices
  A  = new double[ns1sq];
  C  = new double[ns1sq];
  T  = new double[ns1sq];
  S  = new double[ns1sq];
  TT = new double[ns1sq];
  ST = new double[ns1sq];

  // start temperature (t ramp)
  t_start = utils::numeric(FLERR,arg[4],false,lmp);

  // final temperature (t ramp)
  t_stop = utils::numeric(FLERR,arg[5],false,lmp);

  // PRNG seed
  int seed = utils::inumeric(FLERR,arg[6],false,lmp);

  // LOADING A matrix
  FILE *fgle = nullptr;
  char *fname = arg[7];
  if (comm->me == 0) {
    fgle = utils::open_potential(fname,lmp,nullptr);
    if (fgle == nullptr)
      error->one(FLERR,"Cannot open A-matrix file {}: {}",fname, utils::getsyserror());
    utils::logmesg(lmp,"Reading A-matrix from {}\n", fname);
  }

  // read each line of the file, skipping blank lines or leading '#'

  char line[MAXLINE],*ptr;
  int n,nwords,ndone=0,eof=0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fgle);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fgle);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';

    nwords = utils::count_words(line);
    if (nwords == 0) continue;

    ptr = strtok(line," \t\n\r\f");
    do {
      A[ndone] = atof(ptr);
      ptr = strtok(nullptr," \t\n\r\f");
      ndone++;
    } while ((ptr != nullptr) && (ndone < ns1sq));
  }

  fnoneq=0; gle_every=1; gle_step=0;
  for (int iarg=8; iarg<narg; iarg+=2) {
    if (strcmp(arg[iarg],"noneq") == 0) {
      fnoneq = 1;
      if (iarg+2>narg)
        error->all(FLERR,"Did not specify C matrix for non-equilibrium GLE");

      fname = arg[iarg+1];

    } else if (strcmp(arg[iarg],"every") == 0) {

      if (iarg+2>narg)
        error->all(FLERR, "Did not specify interval for applying the GLE");
      gle_every=utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    }
  }

  // set C matrix

  if (fnoneq == 0) {
    t_target=t_start;
    const double kT = t_target * force->boltz / force->mvv2e;
    memset(C,0,sizeof(double)*ns1sq);
    for (int i=0; i<ns1sq; i+=(ns+2))
      C[i]=kT;

  } else {
    if (comm->me == 0) {
      fgle = utils::open_potential(fname,lmp,nullptr);
      if (fgle == nullptr)
        error->one(FLERR,"Cannot open C-matrix file {}: {}",fname, utils::getsyserror());
      utils::logmesg(lmp,"Reading C-matrix from {}\n", fname);
    }

    // read each line of the file, skipping blank lines or leading '#'
    ndone = eof = 0;
    const double cfac = force->boltz / force->mvv2e;

    while (1) {
      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fgle);
        if (ptr == nullptr) {
          eof = 1;
          fclose(fgle);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);

      // strip comment, skip line if blank

      if ((ptr = strchr(line,'#'))) *ptr = '\0';

      nwords = utils::count_words(line);
      if (nwords == 0) continue;

      ptr = strtok(line," \t\n\r\f");
      do {
        C[ndone] = cfac*atof(ptr);
        ptr = strtok(nullptr," \t\n\r\f");
        ndone++;
      } while ((ptr != nullptr) && (ndone < ns1sq));
    }
  }

#ifdef GLE_DEBUG
  printf("A Matrix\n");
  GLE::MyPrint(ns+1,A);
  printf("C Matrix\n");
  GLE::MyPrint(ns+1,C);
#endif

  // initialize Marsaglia RNG with processor-unique seed
  // NB: this means runs will not be the same with different numbers of processors
  if (seed <= 0) error->all(FLERR,"Illegal fix gle command");
  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for mass-scaling
  sqrt_m=nullptr;
  memory->grow(sqrt_m, atom->ntypes+1,"gle:sqrt_m");

  // allocates space for additional degrees of freedom
  gle_s=nullptr;
  // allocates space for temporaries
  gle_tmp1=gle_tmp2=nullptr;

  grow_arrays(atom->nmax);
  init_gles();

  // add callbacks to enable restarts
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  energy = 0.0;
}


/* --- Frees up memory used by temporaries and buffers ------------------ */

FixGLE::~FixGLE()
{
  delete random;
  delete [] A;
  delete [] C;
  delete [] S;
  delete [] T;
  delete [] TT;
  delete [] ST;

  memory->destroy(sqrt_m);
  memory->destroy(gle_s);
  memory->destroy(gle_tmp1);
  memory->destroy(gle_tmp2);
}

/* ---------------------------------------------------------------------- */

int FixGLE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ------- Initializes one-time quantities for GLE ---------------------- */

void FixGLE::init()
{

  dogle = 1;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // set force prefactors
  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      sqrt_m[i] = sqrt(atom->mass[i]);
    }
  }

  if (utils::strmatch(update->integrate_style,"^respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }

  init_gle();
}

/* ------- Initializes integrator matrices (change with T and dt) ------- */

void FixGLE::init_gle()
{
  // compute Langevin terms

  double *tmp1 = new double[ns1sq];
  double *tmp2 = new double[ns1sq];

  for (int i=0; i<ns1sq; ++i) {
    tmp1[i]=-A[i]*update->dt*0.5*gle_every;
    tmp2[i]=S[i]=0.0;
  }
  GLE::MatrixExp(ns+1,tmp1,T);

  GLE::MyMult(ns+1,ns+1,ns+1,T,C,tmp1);
  GLE::MyTrans(ns+1,T,tmp2);
  GLE::MyMult(ns+1,ns+1,ns+1,tmp1,tmp2,S);

  for (int i=0; i<ns1sq; ++i) tmp1[i]=C[i]-S[i];

  GLE::StabCholesky(ns+1, tmp1, S);   //!TODO use symmetric square root, which is more stable.

#ifdef GLE_DEBUG
  printf("T Matrix\n");
  GLE::MyPrint(ns+1,T);
  printf("S Matrix\n");
  GLE::MyPrint(ns+1,S);
#endif

  // transposed evolution matrices to have fast index multiplication in gle_integrate
  GLE::MyTrans(ns+1,T,TT);
  GLE::MyTrans(ns+1,S,ST);
  delete[] tmp1;
  delete[] tmp2;
}

/* ------- Sets initial values of additional DOF for free particles ----- */
void FixGLE::init_gles()
{

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *rootC  = new double[ns1sq];
  double *rootCT = new double[ns1sq];
  double *newg   = new double[3*(ns+1)*nlocal];
  double *news   = new double[3*(ns+1)*nlocal];

  GLE::StabCholesky(ns+1, C, rootC);
  GLE::MyTrans(ns+1,rootC,rootCT);

  memset(news,0,sizeof(double)*3*(ns+1)*nlocal);
  for (int i = 0; i < nlocal*3*(ns+1); ++i)
    newg[i] = random->gaussian();

  GLE::AkMult(nlocal*3,ns+1,ns+1, newg, rootCT, news);

  int nk=0; // unpacks temporary into gle_s
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // first loads velocities
      for (int k = 0; k<3; k++) {
        for (int j=0; j<ns; ++j)
          gle_s[i][k*ns+j]=news[nk++];
      }
    }
  }
  delete[] rootC;
  delete[] rootCT;
  delete[] news;
  delete[] newg;
  return;
}

/* ---------------------------------------------------------------------- */

void FixGLE::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixGLE::gle_integrate()
{
  double **v = atom->v;
  double *rmass = atom->rmass, smi, ismi;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

#ifdef GLE_DEBUG
  printf("!MC! GLE THERMO STEP dt=%f\n",update->dt);
#endif

  // loads momentum data (mass-scaled) into the temporary vectors for the propagation
  int nk=0, ni=0; double deltae=0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ni++;
      if (rmass) {
        smi = sqrt(rmass[i]);
      } else {
        smi=sqrt_m[type[i]];
      }
      for (int k = 0; k<3; k++) {
        //also zeroes tmp1
        gle_tmp1[nk]=0.0;
        // first loads velocities and accumulates conserved quantity
        gle_tmp2[nk]=v[i][k]*smi; deltae+= gle_tmp2[nk]*gle_tmp2[nk]; ++nk;
        // then copies in the additional momenta
        for (int j=0; j<ns; ++j)
          gle_tmp2[nk++] = gle_s[i][k*ns+j];
      }
    }
  }

  // s(t+dt) = T.s ....
  GLE::AkMult(ni*3,ns+1,ns+1,gle_tmp2,TT,gle_tmp1);

  //fills up a vector of random numbers
  for (int i = 0; i < 3*ni*(ns+1); i++) gle_tmp2[i] = random->gaussian();

  // ... + S \xi
  GLE::AkMult(ni*3,ns+1,ns+1,gle_tmp2,ST,gle_tmp1,1.0);

  // unloads momentum data (mass-scaled) from the temporary vectors
  nk=0;
  for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) {
     if (rmass) ismi = 1.0/sqrt(rmass[i]); else ismi=1.0/sqrt_m[type[i]];
     for (int k = 0; k<3; k++)
     {
        // fetches new velocities and completes computation of the conserved quantity change
        v[i][k]=gle_tmp1[nk]*ismi; deltae-= gle_tmp1[nk]*gle_tmp1[nk]; ++nk;
        // stores the additional momenta in the gle_s buffer
        for (int j=0; j<ns; ++j)
          gle_s[i][k*ns+j]=gle_tmp1[nk++];
     }
  }

  energy += deltae*0.5*force->mvv2e;
}

void FixGLE::initial_integrate(int /*vflag*/)
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

  gle_step--;
  if (dogle && gle_step<1) gle_integrate();

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
  if (dogle && gle_step<1) { gle_integrate(); gle_step=gle_every; }

  // Change the temperature for the next step
  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
  if (t_stop != t_start && fnoneq == 0) {
    // only updates if it is really necessary
    const double kT = t_target * force->boltz / force->mvv2e;
    memset(C,0,sizeof(double)*ns1sq);
    for (int i=0; i<ns1sq; i+=(ns+2)) C[i]=kT;
    init_gle();
  }

}
/* ---------------------------------------------------------------------- */

void FixGLE::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel==nlevels_respa-1) gle_integrate();
  dogle=0;
  if (ilevel == 0) initial_integrate(vflag);
  else { final_integrate();}
}

void FixGLE::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dogle=0;
  final_integrate();
  if (ilevel==nlevels_respa-1) gle_integrate();
}


double FixGLE::compute_scalar()
{
  double energy_me = energy;
  double energy_all;
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);

  return energy_all;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixGLE::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  }
  return nullptr;
}


/* ----------------------------------------------------------------------
   Called when a change to the target temperature is requested mid-run
------------------------------------------------------------------------- */

void FixGLE::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
  if (fnoneq == 0)
  {
     // only updates if it is really necessary
     for (int i=0; i<ns1sq; ++i) C[i]=0.0;
     for (int i=0; i<ns1sq; i+=(ns+2)) C[i]=t_target*force->boltz/force->mvv2e;
     init_gle();
  }
  else
  {
     error->all(FLERR, "Cannot change temperature for a non-equilibrium GLE run");
  }
}

/* ----------------------------------------------------------------------
   Called when a change to the timestep is requested mid-run
------------------------------------------------------------------------- */

void FixGLE::reset_dt()
{
  // set the time integration constants
  dtv = update->dt;
  dtf = 0.5 * update->dt * (force->ftm2v);
  init_gle();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixGLE::memory_usage()
{
  double bytes = (double)atom->nmax*(3*ns+2*3*(ns+1))*sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixGLE::grow_arrays(int nmax)
{
  memory->grow(gle_s, nmax, 3*ns,"gle:gle_s");
  memory->grow(gle_tmp1, nmax*(ns+1)*3,"gle:tmp1");
  memory->grow(gle_tmp2, nmax*(ns+1)*3,"gle:tmp2");
  //zeroes out temporary buffers
  for (int i=0; i< nmax*(ns+1)*3; ++i) gle_tmp1[i]=0.0;
  for (int i=0; i< nmax*(ns+1)*3; ++i) gle_tmp2[i]=0.0;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixGLE::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int k = 0; k < 3*ns; k++) gle_s[j][k] = gle_s[i][k];
}

/* ----------------------------------------------------------------------
   Pack extended variables assoc. w/ atom i into buffer for exchange
   with another processor
------------------------------------------------------------------------- */

int FixGLE::pack_exchange(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < 3*ns; k++) buf[m++] = gle_s[i][k];
  return m;
}

/* ----------------------------------------------------------------------
   Unpack extended variables from exchange with another processor
------------------------------------------------------------------------- */

int FixGLE::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < 3*ns; k++) gle_s[nlocal][k] = buf[m++];
  return m;
}


/* ----------------------------------------------------------------------
   Pack extended variables assoc. w/ atom i into buffer for
   writing to a restart file
------------------------------------------------------------------------- */

int FixGLE::pack_restart(int i, double *buf)
{
  int m = 0;
  // pack buf[0] this way because other fixes unpack it
  buf[m++] = 3*ns + 1;
  for (int k = 0; k < 3*ns; k=k+3)
  {
    buf[m++] = gle_s[i][k];
    buf[m++] = gle_s[i][k+1];
    buf[m++] = gle_s[i][k+2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   Unpack extended variables to restart the fix from a restart file
------------------------------------------------------------------------- */

void FixGLE::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to the nth set of extended variables
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i< nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int k = 0; k < 3*ns; k=k+3)
  {
    gle_s[nlocal][k] = extra[nlocal][m++];
    gle_s[nlocal][k+1] = extra[nlocal][m++];
    gle_s[nlocal][k+2] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   Returns the number of items in atomic restart data associated with
   local atom nlocal.  Used in determining the total extra data stored by
   fixes on a given processor.
------------------------------------------------------------------------- */

int FixGLE::size_restart(int /*nlocal*/)
{
  return 3*ns+1;
}

/* ----------------------------------------------------------------------
   Returns the maximum number of items in atomic restart data
   Called in Modify::restart for peratom restart.
------------------------------------------------------------------------- */

int FixGLE::maxsize_restart()
{
  return 3*ns+1;
}
