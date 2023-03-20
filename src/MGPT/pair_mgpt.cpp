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

/* ----------------------------------------------------------------------
   Contributing authors: Tomas Oppelstrup, LLNL (oppelstrup2@llnl.gov)
   and John Moriarty, LLNL (moriarty2@llnl.gov)

   Fast MGPT algorithm developed by Tomas Oppelstrup (2015) based on the
   matrix MGPT v4.4 FORTRAN routine of John Moriarty (2006) as converted
   to C++ for LAMMPS application by Jamie Marian and Alexander Stukowski
   (2011).  See LLNL copyright notice at bottom of this file.
------------------------------------------------------------------------- */

#include "pair_mgpt.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cassert>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

//#define TIMING_ON

#ifdef TIMING_ON
#include <sys/time.h>
#include <time.h>
//#include "rdtsc.h"
#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#endif

static double gettime(int x = 0) {
  if (1) {
    /*
      struct timeval tv;
      gettimeofday(&tv,nullptr);
      return tv.tv_sec + 1e-6 * tv.tv_usec;
    */
    /*
      const double x = 1.0 / CLOCKS_PER_SEC;
      return clock() * x;
    */

    //const double invfreq = 1.0 / 2394.108e6;
    /*
    const double invfreq = 1.0 / 700e6;
      unsigned long long int x = rdtsc();
      return x*invfreq;
    */

    const double invfreq = 1.0 / 1.6e9;
    unsigned long long int x = GetTimeBase();
    return x*invfreq;


  } else
    return 0.0;
}
#else
static double gettime(int /*x*/ = 0) { return 0.0; }
#endif


/* ---------------------------------------------------------------------- */

PairMGPT::PairMGPT(LAMMPS *lmp) : Pair(lmp)
{
        single_enable = 0;
        one_coeff = 1;
        ghostneigh = 1;
}

PairMGPT::~PairMGPT()
{
        if (allocated) {
                memory->destroy(setflag);
                memory->destroy(cutsq);
                memory->destroy(cutghost);
        }
}

/* ---------------------------------------------------------------------- */


static double t_make_b2 = 0.0,n_make_b2 = 0.0;

template<typename intype,typename outtype,int ni,int nj> void fmatconv(intype *array) {
  outtype *cast = (outtype *) array;
  for (int i = 0; i<ni; i++)
    for (int j = 0; j<nj; j++)
      cast[i*nj+j] = array[i*nj+j];
}

void PairMGPT::make_bond(const double xx[][3],int i,int j,bond_data *bptr) {
  double rrij[3],rij;
  int p;

  double t0,t1;

  /* Check that alignment requirements for SIMD code are fulfilled */
  assert( (((unsigned long long int) (bptr->H.m )) & 31) == 0 );
  assert( (((unsigned long long int) (bptr->Hx.m)) & 31) == 0 );
  assert( (((unsigned long long int) (bptr->Hy.m)) & 31) == 0 );
  assert( (((unsigned long long int) (bptr->Hz.m)) & 31) == 0 );

  rij = 0.0;
  for (p = 0; p<3; p++) {
    rrij[p] = xx[i][p] - xx[j][p];
    rij = rij + rrij[p]*rrij[p];
  }

  /* Zero all matrix elements */
  for (i = 0; i<8; i++)
    for (j = 0; j<8; j++) {
      bptr->H.m[i][j] = 0.0;
      bptr->Hx.m[i][j] = 0.0;
      bptr->Hy.m[i][j] = 0.0;
      bptr->Hz.m[i][j] = 0.0;
      bptr->Hz.m[j][i] = 0.0;
    }

  if (rij <= rcrit*rcrit) {
    t0 = gettime();
    if (lang == 3) {
      hamltn_5_raw(rrij[0],rrij[1],rrij[2],
                   bptr->H.m ,bptr->Hx.m,
                   bptr->Hy.m,bptr->Hz.m,&bptr->fl_deriv_sum);
    } else {
      hamltn_7_raw(rrij[0],rrij[1],rrij[2],
                   bptr->H.m ,bptr->Hx.m,
                   bptr->Hy.m,bptr->Hz.m,&bptr->fl_deriv_sum);
    }

    t1 = gettime();
    t_make_b2 += t1-t0;
    n_make_b2++;
  } else {
    bptr->fl_deriv_sum = 0.0;
  }

  if (linalg.single) {
    fmatconv<double,float,7,8>(&(bptr->H.m[1][0]));
    fmatconv<double,float,7,8>(&(bptr->Hx.m[1][0]));
    fmatconv<double,float,7,8>(&(bptr->Hy.m[1][0]));
    fmatconv<double,float,7,8>(&(bptr->Hz.m[1][0]));
  }
}

static double t_trace = 0.0,n_trace = 0.0;
/*
static inline double mtrace(int n,double A[8][8],double B[8][8]) {
  double t0,t1;
  double s;

  t0 = gettime();
  if (n == 5) s = mtrace_5(A,B);
  else if (n == 7) s = mtrace_7(A,B);
  else {
    s = 0.0;
    for (int i = 1; i<=n; i++)
      for (int j = 1; j<=n; j++)
        s = s + A[i][j]*B[i][j];
  }
  t1 = gettime();
  t_trace += t1-t0;
  n_trace++;

  return s;
}
*/

void PairMGPT::make_triplet(bond_data *ij_bond,bond_data *ik_bond,
                             triplet_data *triptr) {
#if 1
  const trmul_fun tr_mul = linalg.tr_mul;
  tr_mul(&(ij_bond->H.m[1][0]), &(ik_bond->H.m[1][0]) ,&(triptr->H1H2.m[1][0]) );
  tr_mul(&(ij_bond->Hx.m[1][0]),&(ik_bond->H.m[1][0]) ,&(triptr->H1xH2.m[1][0]));
  tr_mul(&(ij_bond->Hy.m[1][0]),&(ik_bond->H.m[1][0]) ,&(triptr->H1yH2.m[1][0]));
  tr_mul(&(ij_bond->Hz.m[1][0]),&(ik_bond->H.m[1][0]) ,&(triptr->H1zH2.m[1][0]));
  tr_mul(&(ij_bond->H.m[1][0]) ,&(ik_bond->Hx.m[1][0]),&(triptr->H1H2x.m[1][0]));
  tr_mul(&(ij_bond->H.m[1][0]) ,&(ik_bond->Hy.m[1][0]),&(triptr->H1H2y.m[1][0]));
  tr_mul(&(ij_bond->H.m[1][0]) ,&(ik_bond->Hz.m[1][0]),&(triptr->H1H2z.m[1][0]));
#else
  transprod(ij_bond->H, ik_bond->H ,triptr->H1H2 );
  transprod(ij_bond->Hx,ik_bond->H ,triptr->H1xH2);
  transprod(ij_bond->Hy,ik_bond->H ,triptr->H1yH2);
  transprod(ij_bond->Hz,ik_bond->H ,triptr->H1zH2);
  transprod(ij_bond->H ,ik_bond->Hx,triptr->H1H2x);
  transprod(ij_bond->H ,ik_bond->Hy,triptr->H1H2y);
  transprod(ij_bond->H ,ik_bond->Hz,triptr->H1H2z);
#endif
}

static double t_make_t = 0.0,t_make_b = 0.0,n_make = 0.0;

PairMGPT::triplet_data *PairMGPT::get_triplet(const double xx[][3],int i,int j,int k,
                                                Hash<bond_data,Doublet> *bhash,
                                                triplet_data *twork,
                                                double *dvir_ij_p,double *dvir_ik_p) {
  const int recompute = 0;
  static bond_data bij_work,bik_work;

  double t0,t1;

  bond_data *bij = nullptr,*bik = nullptr;
  triplet_data *tptr = nullptr;

  t0 = gettime();
  if (recompute == 0) {
    bij = bhash->Lookup(Doublet(i,j));
    bik = bhash->Lookup(Doublet(i,k));
  }

  if (bij == nullptr) {
    if (recompute == 0)
      bij = bhash->Insert(Doublet(i,j));
    else
      bij = &bij_work;
    if (i < j)
      make_bond(xx,i,j,bij);
    else
      make_bond(xx,j,i,bij);
  }

  if (bik == nullptr) {
    if (recompute == 0)
      bik = bhash->Insert(Doublet(i,k));
    else
      bik = &bik_work;
    if (i < k)
      make_bond(xx,i,k,bik);
    else
      make_bond(xx,k,i,bik);
  }
  t1 = gettime();
  t_make_b += t1-t0;

  t0 = gettime();
  if (bij != nullptr && bik != nullptr) {
    tptr = twork;
    make_triplet(bij,bik,tptr);
    *dvir_ij_p = bij->fl_deriv_sum;
    *dvir_ik_p = bik->fl_deriv_sum;
  } else {
    *dvir_ij_p = 0.0;
    *dvir_ik_p = 0.0;
  }
  t1 = gettime();
  t_make_t += t1-t0;
  n_make++;
  return tptr;
}


double PairMGPT::numderiv3t(double xx[][3],int i,int j,int k,int p) {
  static bond_data Bij,Bjk,Bki;
  const double delta = 1e-5;

  const double xsave = xx[i][p];
  double e1,e2;

  const double vc = splinepot.vc;

  xx[i][p] = xsave + delta;
  make_bond(xx,i,j,&Bij);
  make_bond(xx,j,k,&Bjk);
  make_bond(xx,k,i,&Bki);
  e1 = trace(prodmat(Bij.H,Bjk.H),Bki.H) * (vc/anorm3);

  xx[i][p] = xsave - delta;
  make_bond(xx,i,j,&Bij);
  if (false) { /* This bond doesn't change when i is perturbed */
    make_bond(xx,j,k,&Bjk);
  }
  make_bond(xx,k,i,&Bki);
  e2 = trace(prodmat(Bij.H,Bjk.H),Bki.H) * (vc/anorm3);

  xx[i][p] = xsave;
  return (e1 - e2)/(2.0*delta);
}

double PairMGPT::numderiv3v(double xx[][3],int i,int j,int k,int p,int ipert) {
  static bond_data Bij,Bik;
  const double delta = 1e-5;

  const double xsave = xx[ipert][p];
  double e1,e2;

  const double vd = splinepot.vd;

  xx[ipert][p] = xsave + delta;
  make_bond(xx,i,j,&Bij);
  make_bond(xx,i,k,&Bik);
  e1 = trace(prodmat(Bij.H,Bij.H),prodmat(Bik.H,Bik.H)) * (vd/anorm4);

  xx[ipert][p] = xsave - delta;
  make_bond(xx,i,j,&Bij);
  make_bond(xx,i,k,&Bik);
  e2 = trace(prodmat(Bij.H,Bij.H),prodmat(Bik.H,Bik.H)) * (vd/anorm4);

  xx[ipert][p] = xsave;
  return (e1 - e2)/(2.0*delta);
}


double PairMGPT::numderiv4(double xx[][3],int i,int j,int k,int m,int p) {
  static bond_data Bij,Bjk,Bkm,Bmi;
  const double delta = 1e-5;

  const double xsave = xx[i][p];
  double e1,e2;

  const double ve = splinepot.ve;

  xx[i][p] = xsave + delta;
  make_bond(xx,i,j,&Bij);
  make_bond(xx,j,k,&Bjk);
  make_bond(xx,k,m,&Bkm);
  make_bond(xx,m,i,&Bmi);
  e1 = trace(prodmat(Bij.H,Bjk.H),prodmat(Bkm.H,Bmi.H)) * (ve/anorm4);

  xx[i][p] = xsave - delta;
  make_bond(xx,i,j,&Bij);
  if (false) { /* Only the i coordinates changed... */
    make_bond(xx,j,k,&Bjk);
    make_bond(xx,k,m,&Bkm);
  }
  make_bond(xx,m,i,&Bmi);
  e2 = trace(prodmat(Bij.H,Bjk.H),prodmat(Bkm.H,Bmi.H)) * (ve/anorm4);

  xx[i][p] = xsave;
  return (e1 - e2)/(2.0*delta);
}

static double dtol = 1e-6;
void PairMGPT::force_debug_3t(double xx[][3],
                    int i0,int j0,int k0,
                    int i ,int j ,int k ,
                    double dfix,double dfiy,double dfiz,
                    double dfjx,double dfjy,double dfjz,
                    double dfkx,double dfky,double dfkz) {
  double dfi[3],dfj[3],dfk[3];
  dfi[0] = dfix; dfi[1] = dfiy; dfi[2] = dfiz;
  dfj[0] = dfjx; dfj[1] = dfjy; dfj[2] = dfjz;
  dfk[0] = dfkx; dfk[1] = dfky; dfk[2] = dfkz;

  for (int p = 0; p<3; p++) {
    /* Compute numerical derivatives by displacing atoms i,j,k */
    double ndfi,ndfj,ndfk;
    ndfi = -numderiv3t(xx,i,j,k,p);
    ndfj = -numderiv3t(xx,j,k,i,p);
    ndfk = -numderiv3t(xx,k,i,j,p);

    if ((fabs(dfi[p] - ndfi) > dtol &&
        fabs(dfi[p] - ndfi) > dtol*fabs(ndfi)) ||
       (fabs(dfj[p] - ndfj) > dtol &&
        fabs(dfj[p] - ndfj) > dtol*fabs(ndfj)) ||
       (fabs(dfk[p] - ndfk) > dtol &&
        fabs(dfk[p] - ndfk) > dtol*fabs(ndfk))) {
      printf("Force error in T12 & T23 & T31 :: i,j,k = %d,%d,%d\n",i0,j0,k0);
      printf("    dE/d%c[i] = %20.10e    %20.10e\n", 'x'+p,ndfi, dfi[p]);
      printf("    dE/d%c[j] = %20.10e    %20.10e\n", 'x'+p,ndfj, dfj[p]);
      printf("    dE/d%c[k] = %20.10e    %20.10e\n", 'x'+p,ndfk, dfk[p]);
      printf("\n");
    }
  }
}

void PairMGPT::force_debug_3v(double xx[][3],
                    int i0,int j0,int k0,
                    int i ,int j ,int k ,
                    double dfix,double dfiy,double dfiz,
                    double dfjx,double dfjy,double dfjz,
                    double dfkx,double dfky,double dfkz) {
  double dfi[3],dfj[3],dfk[3];
  dfi[0] = dfix; dfi[1] = dfiy; dfi[2] = dfiz;
  dfj[0] = dfjx; dfj[1] = dfjy; dfj[2] = dfjz;
  dfk[0] = dfkx; dfk[1] = dfky; dfk[2] = dfkz;

  for (int p = 0; p<3; p++) {
    /* Compute numerical derivatives by displacing atoms i,j,k */
    double ndfi,ndfj,ndfk;
    ndfi = -numderiv3v(xx,i,j,k,p,i0);
    ndfj = -numderiv3v(xx,i,j,k,p,j0);
    ndfk = -numderiv3v(xx,i,j,k,p,k0);

    if ((fabs(dfi[p] - ndfi) > dtol &&
        fabs(dfi[p] - ndfi) > dtol*fabs(ndfi)) ||
       (fabs(dfj[p] - ndfj) > dtol &&
        fabs(dfj[p] - ndfj) > dtol*fabs(ndfj)) ||
       (fabs(dfk[p] - ndfk) > dtol &&
        fabs(dfk[p] - ndfk) > dtol*fabs(ndfk))) {
      printf("Force error in T12 :: i,j,k = %d,%d,%d\n",i0,j0,k0);
      printf("    dE/d%c[i] = %20.10e    %20.10e\n", 'x'+p,ndfi, dfi[p]);
      printf("    dE/d%c[j] = %20.10e    %20.10e\n", 'x'+p,ndfj, dfj[p]);
      printf("    dE/d%c[k] = %20.10e    %20.10e\n", 'x'+p,ndfk, dfk[p]);
      printf("\n");
    }
  }
}

void PairMGPT::force_debug_4(double xx[][3],
                   int i0,int j0,int k0,int m0,
                   int i ,int j ,int k ,int m ,
                   double dfix,double dfiy,double dfiz,
                   double dfjx,double dfjy,double dfjz,
                   double dfkx,double dfky,double dfkz,
                   double dfmx,double dfmy,double dfmz) {

  double dfi[3],dfj[3],dfk[3],dfm[3];
  dfi[0] = dfix; dfi[1] = dfiy; dfi[2] = dfiz;
  dfj[0] = dfjx; dfj[1] = dfjy; dfj[2] = dfjz;
  dfk[0] = dfkx; dfk[1] = dfky; dfk[2] = dfkz;
  dfm[0] = dfmx; dfm[1] = dfmy; dfm[2] = dfmz;

  const int ii0[] = {i0,j0,k0,m0},ii[] = {i,j,k,m,i,j,k};

  for (int p = 0; p<3; p++) {
    /* Compute numerical derivatives by displacing atoms i,j,k,m */
    double ndfi,ndfj,ndfk,ndfm;
    if (true) {
      double ndf[] = {0.0,0.0,0.0,0.0};
      for (int s = 0; s<4; s++)
        for (int t = 0; t<4; t++)
          if (ii[s] == ii0[t])
            ndf[t] = -numderiv4(xx,ii[s],ii[s+1],ii[s+2],ii[s+3],p);
      ndfi = ndf[0]; ndfj = ndf[1];
      ndfk = ndf[2]; ndfm = ndf[3];
    } else {
      ndfi = -numderiv4(xx,i,j,k,m,p);
      ndfj = -numderiv4(xx,j,k,m,i,p);
      ndfk = -numderiv4(xx,k,m,i,j,p);
      ndfm = -numderiv4(xx,m,i,j,k,p);
    }

    if ((fabs(dfi[p] - ndfi) > dtol &&
        fabs(dfi[p] - ndfi) > dtol*fabs(ndfi)) ||
       (fabs(dfj[p] - ndfj) > dtol &&
        fabs(dfj[p] - ndfj) > dtol*fabs(ndfj)) ||
       (fabs(dfk[p] - ndfk) > dtol &&
        fabs(dfk[p] - ndfk) > dtol*fabs(ndfk)) ||
       (fabs(dfm[p] - ndfm) > dtol &&
        fabs(dfm[p] - ndfm) > dtol*fabs(ndfm))) {
      printf("Force error in T31 & T64 :: i,j,k,m = %d,%d,%d,%d\n",i0,j0,k0,m0);
      printf("    dE/d%c[i] = %20.10e    %20.10e\n", 'x'+p,ndfi, dfi[p]);
      printf("    dE/d%c[j] = %20.10e    %20.10e\n", 'x'+p,ndfj, dfj[p]);
      printf("    dE/d%c[k] = %20.10e    %20.10e\n", 'x'+p,ndfk, dfk[p]);
      printf("    dE/d%c[m] = %20.10e    %20.10e\n", 'x'+p,ndfm, dfm[p]);
      printf("\n");
    }
  }
}



/*
#define trd_update_4(T12,T45,coord)                   \
  do {                                                \
    trd1 = transtrace(T12->H1##coord##H2,T45->H1H2 ); \
    trd2 = transtrace(T12->H1H2##coord,T45->H1H2   ); \
    trd3 = transtrace(T12->H1H2 ,T45->H1##coord##H2); \
    trd4 = transtrace(T12->H1H2 ,T45->H1H2##coord  ); \
  } while (0)
*/
#define trd_update_4(T12,T45) \
  do {                                         \
    tr_trace3(&(T45->H1H2.m[1][0]),            \
              &(T12->H1xH2.m[1][0]),&utr1x.d,  \
              &(T12->H1yH2.m[1][0]),&utr1y.d,  \
              &(T12->H1zH2.m[1][0]),&utr1z.d); \
    tr_trace3(&(T45->H1H2.m[1][0]),            \
              &(T12->H1H2x.m[1][0]),&utr2x.d,  \
              &(T12->H1H2y.m[1][0]),&utr2y.d,  \
              &(T12->H1H2z.m[1][0]),&utr2z.d); \
    tr_trace3(&(T12->H1H2.m[1][0]),            \
              &(T45->H1xH2.m[1][0]),&utr3x.d,  \
              &(T45->H1yH2.m[1][0]),&utr3y.d,  \
              &(T45->H1zH2.m[1][0]),&utr3z.d); \
    tr_trace3(&(T12->H1H2.m[1][0]),            \
              &(T45->H1H2x.m[1][0]),&utr4x.d,  \
              &(T45->H1H2y.m[1][0]),&utr4y.d,  \
              &(T45->H1H2z.m[1][0]),&utr4z.d); \
    if (linalg.single) {                        \
      trd1x = utr1x.f; trd2x = utr2x.f; trd3x = utr3x.f; trd4x = utr4x.f; \
      trd1y = utr1y.f; trd2y = utr2y.f; trd3y = utr3y.f; trd4y = utr4y.f; \
      trd1z = utr1z.f; trd2z = utr2z.f; trd3z = utr3z.f; trd4z = utr4z.f; \
    } else {                                                              \
      trd1x = utr1x.d; trd2x = utr2x.d; trd3x = utr3x.d; trd4x = utr4x.d; \
      trd1y = utr1y.d; trd2y = utr2y.d; trd3y = utr3y.d; trd4y = utr4y.d; \
      trd1z = utr1z.d; trd2z = utr2z.d; trd3z = utr3z.d; trd4z = utr4z.d; \
    }                                          \
  } while (0)

#define dfix_update_4a(coord) \
  do { \
    dfi##coord = ( (-sij)*trd1##coord + (-sim)*trd3##coord ) * (ve / anorm4); \
    dfj##coord = ( ( sij)*trd1##coord + (-sjk)*trd2##coord ) * (ve / anorm4); \
    dfk##coord = ( ( sjk)*trd2##coord + (-skm)*trd4##coord ) * (ve / anorm4); \
    dfm##coord = ( ( sim)*trd3##coord + ( skm)*trd4##coord ) * (ve / anorm4); \
  } while (0)


#define dfix_update_4b(coord) \
  do { \
    dfi##coord = ( ( ski)*trd1##coord + (-sim)*trd3##coord ) * (ve / anorm4); \
    dfj##coord = ( (-sjk)*trd2##coord + (-sjm)*trd4##coord ) * (ve / anorm4); \
    dfk##coord = ( (-ski)*trd1##coord + ( sjk)*trd2##coord ) * (ve / anorm4); \
    dfm##coord = ( ( sim)*trd3##coord + ( sjm)*trd4##coord ) * (ve / anorm4); \
  } while (0);

#define dfix_update_4c(coord) \
  do { \
    dfi##coord = ( (-sij)*trd1##coord + ( ski)*trd2##coord ) * (ve / anorm4); \
    dfj##coord = ( ( sij)*trd1##coord + (-sjm)*trd3##coord ) * (ve / anorm4); \
    dfk##coord = ( (-ski)*trd2##coord + (-skm)*trd4##coord ) * (ve / anorm4); \
    dfm##coord = ( ( sjm)*trd3##coord + ( skm)*trd4##coord ) * (ve / anorm4); \
  } while (0);

#define accumulate_forces_2(w) \
  do { \
    fix = fix + dfix*(w);    \
    fiy = fiy + dfiy*(w);    \
    fiz = fiz + dfiz*(w);    \
                             \
    fjx = fjx + dfjx*(w);    \
    fjy = fjy + dfjy*(w);    \
    fjz = fjz + dfjz*(w);    \
  } while (0)

#define accumulate_forces_3(w) \
  do { \
    accumulate_forces_2(w); \
    fkx = fkx + dfkx*(w);    \
    fky = fky + dfky*(w);    \
    fkz = fkz + dfkz*(w);    \
  } while (0)

#define accumulate_forces_4(w) \
  do { \
    accumulate_forces_3(w); \
    fmx = fmx + dfmx*(w);    \
    fmy = fmy + dfmy*(w);    \
    fmz = fmz + dfmz*(w);    \
  } while (0)



#define restrict __restrict__
#ifdef __bg__
#define const
#endif
static int ntr_calls = 0;
static trtrace3_fun tr_internal;
static void tr_count(const double * restrict A,
                     const double * restrict B1,double * restrict t1,
                     const double * restrict B2,double * restrict t2,
                     const double * restrict B3,double * restrict t3) {
  tr_internal(A,B1,t1,B2,t2,B3,t3);
  ntr_calls++;
}
#ifdef __bg__
#undef const
#endif
#undef restrict


int PairMGPT::Matrix::sz;
void PairMGPT::compute_x(const int *nnei,const int * const *nlist,
                          double *e_s,double *e_p,double *e_t,double *e_q,
                          int evflag,int newton_pair) {
  Hash<bond_data,Doublet> bond_hash(100000);
  int i,j,k,m,ix,jx,kx,mx,itag,jtag,p;

  double e_single,e_pair,e_triplet,e_triplet_c,e_quad;
  double volvir2;

  double nbc = 0.0,tbl = 0.0,tbm = 0.0;
  const int lmax_local = lmax;

  //if(evflag) printf("##### ev flag is set... wasting cycles...\n");

  *e_s = -99.0;
  *e_p = -99.0;
  *e_t = -99.0;
  *e_q = -99.0;

  double t0,t1;

  t0 = gettime(1);
  e_single = e_pair = e_triplet = e_triplet_c = e_quad = 0.0;
  volvir2 = 0.0;

  t_make_t = t_make_b = t_make_b2 = t_trace = 0.0;
  n_make = n_make_b2 = n_trace =  0.0;

  double tx0,tx1,tsort = 0.0,tpair = 0.0,tlookup = 0.0;
  double ttriplet = 0.0,tquad = 0.0,tmem = 0.0;
  double ntsort = 0.0,ntpair = 0.0,ntlookup = 0.0;
  double nttriplet = 0.0,ntquad = 0.0,ntmem = 0.0,ntquaditer = 0.0;
  double mcount = 0.0,mcount2 = 0.0, qcount = 0.0;

  double fix,fjx,fkx,fmx,dfix,dfjx,dfkx,dfmx;
  double fiy,fjy,fky,fmy,dfiy,dfjy,dfky,dfmy;
  double fiz,fjz,fkz,fmz,dfiz,dfjz,dfkz,dfmz;

  double fsave[4][3] = { {0.0} } /* {{0.0}} is to get rid of uninitialized use warning */;

  //const int numerical_pair_forces = (nbody_flag/16)%2;
  const int pair_forces   = (nbody_flag/2)%2,three_body_forces   = (nbody_flag/4)%2,four_body_forces   = (nbody_flag/8)%2;
  const int pair_energies = (nbody_flag/2)%2,three_body_energies = (nbody_flag/4)%2,four_body_energies = (nbody_flag/8)%2;
  const int single_energies = nbody_flag%2;

  const int triplet_debug = 0,quad_debug = 0;

  /* Energy and force scale factor for unit conversion. */
  const double e_scale = 0.5;

#ifdef NEIGHMASK
#define NIDX(x) (x)
#else
#define NIDX(x) ((x) & NEIGHMASK)
#endif

  int nneitot,*first,*nlist_short;
  double w2,w3,w4;
  triplet_data T12work,T23work,T31work,T45work,T56work,T64work;
  triplet_data *T12,*T23,*T31,*T45,*T56,*T64;
  int c_ij,c_jk,c_ki,c_im,c_jm,c_km;
  int mi,mj,mk;
  double tr0,tr1,tr2,tr3;
  double v33,v43;
  double rcut2_pair = rmax*rmax,rcut2_bond = rcrit*rcrit,rij2;
  int ntot,nloc;

  double dvir_ij,dvir_jk,dvir_ki,dvir_im,dvir_jm,dvir_km;
  double vir3t = 0.0,vir3v = 0.0,vir4 = 0.0;

  double (*xx)[3],(*ff)[3],(*ss)[3];

#ifdef TIMING_ON
  tr_internal = linalg.tr_trace; ntr_calls = 0;
  const trtrace3_fun tr_trace3 = tr_count;
#else
  const trtrace3_fun tr_trace3 = linalg.tr_trace;
#endif
  union {
    double d;
    float f;
  } utr1x,utr2x,utr3x,utr4x,utr1y,utr2y,utr3y,utr4y,utr1z,utr2z,utr3z,utr4z;
  double trd1x,trd2x,trd3x,trd4x;
  double trd1y,trd2y,trd3y,trd4y;
  double trd1z,trd2z,trd3z,trd4z;


  tx0 = gettime();

  double rhoinv;
  {
    double vtot = 1.0;
    double ntot = atom->natoms;
    for (i = 0; i<3; i++)
      vtot = vtot * (domain->boxhi[i] - domain->boxlo[i]);
    rhoinv = vtot / ntot;
  }

  /* Make sure triplet data work area is aligned and zeroed out. */ {
    assert(T12work.align_check() == 0);
    assert(T23work.align_check() == 0);
    assert(T31work.align_check() == 0);
    assert(T45work.align_check() == 0);
    assert(T56work.align_check() == 0);
    assert(T64work.align_check() == 0);

    T12work.zero();  T23work.zero();  T31work.zero();
    T45work.zero();  T56work.zero();  T64work.zero();
  }

  ntot = atom->nlocal + atom->nghost;
  nloc = atom->nlocal;
  //printf("[%3d] Allocating local array, size is %d atoms...\n",comm->me,j);
  xx = (double (*)[3]) memory->smalloc(sizeof(double [3]) * ntot,"mgpt: local position vector.");
  ff = (double (*)[3]) memory->smalloc(sizeof(double [3]) * ntot,"mgpt: local force vector.");
  //printf("[%3d] Initializing arrays...\n",comm->me);

  const int triclinic = domain->triclinic;
  double alpha[3] = {0.0,0.0,0.0};
  if (triclinic) {
    double E[3][3],EX[3][3];
    int cyc[] = {0,1,2,0,1};

    ss = (double (*)[3]) memory->smalloc(sizeof(double [3]) * ntot,
                                         "mgpt: local reduced coordinate vector.");

    for (i = 0; i<3; i++) {
      for (j = 0; j<3; j++)
        E[i][j] = 0.0;
      E[i][i] = domain->subhi_lamda[i] - domain->sublo_lamda[i];
      domain->lamda2x(E[i],EX[i]);
    }
    for (i = 0; i<3; i++) {
      int i1 = cyc[i+1],i2 = cyc[i+2];
      double dot = 0.0,ns2 = 0.0;
      for (j = 0; j<3; j++) {
        int j1 = cyc[j+1],j2 = cyc[j+2];
        double cj = EX[i1][j1]*EX[i2][j2] - EX[i1][j2]*EX[i2][j1];
        ns2 = ns2 + cj*cj;
        dot = dot + EX[i][j]*cj;
      }
      alpha[i] = E[i][i] / (dot/sqrt(ns2));
      if (comm->me == 0) {
        static int count = 0;
        if (count < 3)
          printf("@@@ alpha(%d) = %15.5e\n",i+1,alpha[i]);
        count++;
      }
      if (alpha[i] < 0.0) alpha[i] = -alpha[i];
    }
  } else
    ss = xx;

  nneitot = 0;
  for (ix = 0; ix<ntot; ix++) {
    for (p = 0; p<3; p++) {
      xx[ix][p] = atom->x[ix][p];
      ff[ix][p] = 0.0;
    }
    if (triclinic)
      domain->x2lamda(xx[ix],ss[ix]);
    nneitot = nneitot + nnei[ix];
  }

  first = (int *) memory->smalloc(sizeof(int) * (ntot+1),"mgpt: first");
  nlist_short = (int *) memory->smalloc(sizeof(int) * nneitot,"mgpt: nlist_short");

  tx1 = gettime();
  tmem += tx1-tx0;
  ntmem++;

  //printf("[%3d] Starting calculation...\n",comm->me);


  fix = fjx = fkx = fmx = 0.0;
  fiy = fjy = fky = fmy = 0.0;
  fiz = fjz = fkz = fmz = 0.0;

  int c_p = 0, c_t = 0, c_q = 0;

  if (false)
    if (domain->triclinic) {
      if (comm->me == 0)
        printf("Can not handle triclinic box yet\n");
      error->all(FLERR,"Can not handle triclinic cell with mgpt yet.");
    }

  /*
  for (i = 0; i<nloc; i++) {
    printf("Atom %3d:: %10.3f  %10.3f  %10.3f\n",
           i,xx[i][0],xx[i][1],xx[i][2]);
  }
  */

  first[0] = 0;
  for (i = 0; i<ntot; i++) {
    fix = fiy = fiz = 0.0;

    first[i+1] = first[i];

    const int c1 = c1_outside(ss[i],triclinic,alpha);

    tx0 = gettime();
    for (jx = 0; jx<nnei[i]; jx++) {
      fjx = fjy = fjz = 0.0;

      j = NIDX( nlist[i][jx] );

      rij2 = 0.0;
      for (p = 0; p<3; p++) {
        double t = xx[i][p] - xx[j][p];
        rij2 = rij2 + t*t;
      }

      if (c1 == 0 && rij2 < rcut2_pair) {
        if (j < i) {
          w2 = get_weight(triclinic,ss[i],ss[j]);

          if (w2 > 0.0) {
            /*
              Compute pair energy/force
            */
            double de_pair,df,rij = sqrt(rij2);
            splinepot.eval_pot(rij,&de_pair,&df);
            de_pair = de_pair * e_scale * w2;
            df = df / rij * w2;


            if (pair_energies == 0) de_pair = 0.0;
            e_pair = e_pair + de_pair;
            c_p++;

            if (pair_forces == 0) df = 0.0;

            if (volpres_flag && pair_energies) {
              double dvir;
              splinepot.eval_vir(rij,&dvir);
              volvir2 = volvir2 - dvir * w2;

              /* Per-atom virial contribution of volumetric energy term */
              if (vflag_atom)
                for (int pp = 0; pp<3; pp++) {
                  //virial[i] = virial[i] + rhoinv*e_scale*volvir2;
                  vatom[i][pp] -= 0.5 * rhoinv*e_scale*dvir*w2;
                  vatom[j][pp] -= 0.5 * rhoinv*e_scale*dvir*w2;
                }
            }

            double drijx = xx[j][0] - xx[i][0];
            double drijy = xx[j][1] - xx[i][1];
            double drijz = xx[j][2] - xx[i][2];

            fix = fix + df*drijx;
            fjx = fjx - df*drijx;

            fiy = fiy + df*drijy;
            fjy = fjy - df*drijy;

            fiz = fiz + df*drijz;
            fjz = fjz - df*drijz;

            if (evflag) {
              //ev_tally(i,j,nloc,newton_pair,de_pair,0.0,df,-drijx,-drijy,-drijz);
              /* To fix stress-per-atom scaling, and sign */
              ev_tally(i,j,nloc,newton_pair,de_pair,0.0,-df * e_scale,-drijx,-drijy,-drijz);
            }

            ff[j][0] += fjx * e_scale;
            ff[j][1] += fjy * e_scale;
            ff[j][2] += fjz * e_scale;

          }
        }
      }

      if (rij2 < rcut2_bond && c2_outside(ss[i],ss[j],triclinic,alpha) == 0) {
        /*
          Add j to short neighbor list for i.
          Insert j to keep list sorted.
        */

        p = first[i+1]-1;
        while (p >= first[i] && nlist_short[p] > j) {
          nlist_short[p+1] = nlist_short[p];
          p = p - 1;
        }
        nlist_short[p+1] = j;
        first[i+1] = first[i+1] + 1;
        if (first[i+1] > nneitot) {
          printf("nneitot = %d, short list full. i=%d\n",
                 nneitot,i);
          error->one(FLERR,"Shit! Short list full\n");
        }

      }
    }

    ff[i][0] += fix * e_scale;
    ff[i][1] += fiy * e_scale;
    ff[i][2] += fiz * e_scale;

    tx1 = gettime();
    tpair += tx1-tx0;
    ntpair += nnei[i];
  }

  for (i = 0; i<ntot; i++) {
    fix = fiy = fiz = 0.0;

    /*
      Use short lists for triplets and quadruplets.
      For open (2-bonded) triplets, can only use k<j, but not k<i.
      For closed (3-bonded) triplets, we can assume k<j<i.

      Quadruplets:
      Always use k<j<i, and require m<i.
      If 5-bonded with im bond, ignore the quadruplet.
      If 6-bonded, require m<k.

      For 4-bonded quadruplets, we can still use k<j, but also
      assume max(m,j)<i
    */

    if (three_body_energies || three_body_forces ||
       four_body_energies || four_body_forces)
      for (jx = first[i]; jx<first[i+1]; jx++) {
        fjx = fjy = fjz = 0.0;

        j = nlist_short[jx];

        for (kx = first[i]; kx<jx; kx++) {
          fkx = fky = fkz = 0.0;

          k = nlist_short[kx];

          /*
            Search lists of j and k, and see if
            1) j is in k-list (closed triplet)
            2) j and k have a common neighbor (closed quadruplet)
          */

          c_ij = c_ki = 1;

          const int sij = (i < j) ? 1 : -1;
          const int sjk = (j < k) ? 1 : -1;
          const int ski = (k < i) ? 1 : -1;


          T12 = T23 = T31 = nullptr;

          mj = first[j];
          /*
            Since i is in the j-list, and i > k and the list
            is sorted, the loop below terminates:-)
          */
          while (mj < first[j+1] && nlist_short[mj] < k) mj = mj + 1;
          if (mj < first[j+1] && nlist_short[mj] == k) {
            /* Closed triplet */
            c_jk = 1;

            if (j > i) continue; /* Require k<j<i for closed triplets */
          } else {
            /* Open triplet */
            c_jk = 0;
          }

          tx0 = gettime();

          w3 = get_weight(triclinic,ss[i],ss[j],ss[k]);

          int triplet_defer;
          if (w3 > 0.0) {
            triplet_defer = 0;

            dvir_ij = dvir_jk = dvir_ki = 0.0;
            if (c_ij && c_jk)
              T12 = get_triplet(xx,j,i,k,&bond_hash,&T12work,&dvir_ij,&dvir_jk);
            if (c_ki && c_jk)
              T23 = get_triplet(xx,k,i,j,&bond_hash,&T23work,&dvir_ki,&dvir_jk);
            if (c_ij && c_ki)
              T31 = get_triplet(xx,i,j,k,&bond_hash,&T31work,&dvir_ij,&dvir_ki);

            if (evflag) {
              fsave[0][0] = fix; fsave[0][1] = fiy; fsave[0][2] = fiz;
              fsave[1][0] = fjx; fsave[1][1] = fjy; fsave[1][2] = fjz;
              fsave[2][0] = fkx; fsave[2][1] = fky; fsave[2][2] = fkz;
              fix = fiy = fiz = 0.0;
              fjx = fjy = fjz = 0.0;
              fkx = fky = fkz = 0.0;
            }

            tr0 = tr1 = tr2 = tr3 = 0.0;
            double xvir3t,xvir3v;
            xvir3t = xvir3v = 0.0;

            if (T12 && T23) {
              bond_data *bki = bond_hash.Lookup(Doublet(k,i));

              if (three_body_energies && evflag) {
                tr0 = transtrace(T12->H1H2,bki->H);
                double dvir = ((dvir_ij + dvir_jk + bki->fl_deriv_sum)*splinepot.vc +
                                 splinepot.dvc)*tr0*w3/anorm3;
                vir3t = vir3t + dvir;
                xvir3t = xvir3t + dvir;
              }
              mcount2++;

              {
                const double vc = splinepot.vc;
                tr_trace3(&(bki->H.m[1][0]),
                            &(T12->H1xH2.m[1][0]),&utr1x.d,
                            &(T12->H1yH2.m[1][0]),&utr1y.d,
                            &(T12->H1zH2.m[1][0]),&utr1z.d);

                tr_trace3(&(bki->H.m[1][0]),
                            &(T12->H1H2x.m[1][0]),&utr2x.d,
                            &(T12->H1H2y.m[1][0]),&utr2y.d,
                            &(T12->H1H2z.m[1][0]),&utr2z.d);

                tr_trace3(&(T12->H1H2.m[1][0]),
                          &(bki->Hx.m[1][0]),&utr3x.d,
                          &(bki->Hy.m[1][0]),&utr3y.d,
                          &(bki->Hz.m[1][0]),&utr3z.d);

                if (linalg.single) {
                  trd1x = utr1x.f; trd2x = utr2x.f; trd3x = utr3x.f;
                  trd1y = utr1y.f; trd2y = utr2y.f; trd3y = utr3y.f;
                  trd1z = utr1z.f; trd2z = utr2z.f; trd3z = utr3z.f;
                } else {
                  trd1x = utr1x.d; trd2x = utr2x.d; trd3x = utr3x.d;
                  trd1y = utr1y.d; trd2y = utr2y.d; trd3y = utr3y.d;
                  trd1z = utr1z.d; trd2z = utr2z.d; trd3z = utr3z.d;
                }

                dfix = ( (-sij)*trd1x + ( ski)*trd3x ) * (vc / anorm3);
                dfjx = ( ( sij)*trd1x + (-sjk)*trd2x ) * (vc / anorm3);
                dfkx = ( ( sjk)*trd2x + (-ski)*trd3x ) * (vc / anorm3);

                dfiy = ( (-sij)*trd1y + ( ski)*trd3y ) * (vc / anorm3);
                dfjy = ( ( sij)*trd1y + (-sjk)*trd2y ) * (vc / anorm3);
                dfky = ( ( sjk)*trd2y + (-ski)*trd3y ) * (vc / anorm3);

                dfiz = ( (-sij)*trd1z + ( ski)*trd3z ) * (vc / anorm3);
                dfjz = ( ( sij)*trd1z + (-sjk)*trd2z ) * (vc / anorm3);
                dfkz = ( ( sjk)*trd2z + (-ski)*trd3z ) * (vc / anorm3);
              }

              if (triplet_debug)
                force_debug_3t(xx,i,j,k, i,j,k,
                               dfix,dfiy,dfiz,
                               dfjx,dfjy,dfjz,
                               dfkx,dfky,dfkz);

              if (three_body_forces)
                accumulate_forces_3(w3);
            }

            if (T12 != nullptr) {
              //printf("T12 i,j,k = %d,%d,%d\n",i,j,k);
              mcount++;
              if (three_body_energies && evflag) {
                tr1 = transtrace(T12->H1H2,T12->H1H2);
                double dvir = (2.0*(dvir_ij + dvir_jk)*splinepot.vd +
                                 splinepot.dvd)*tr1*w3/anorm4;
                vir3v = vir3v + dvir;
                xvir3v = xvir3v + dvir;
              }

              {
                const double vd = splinepot.vd;

                tr_trace3(&(T12->H1H2.m[1][0]),
                          &(T12->H1xH2.m[1][0]),&utr1x.d,
                          &(T12->H1yH2.m[1][0]),&utr1y.d,
                          &(T12->H1zH2.m[1][0]),&utr1z.d);
                tr_trace3(&(T12->H1H2.m[1][0]),
                          &(T12->H1H2x.m[1][0]),&utr2x.d,
                          &(T12->H1H2y.m[1][0]),&utr2y.d,
                          &(T12->H1H2z.m[1][0]),&utr2z.d);
                if (linalg.single) {
                  trd1x = utr1x.f; trd2x = utr2x.f;
                  trd1y = utr1y.f; trd2y = utr2y.f;
                  trd1z = utr1z.f; trd2z = utr2z.f;
                } else {
                  trd1x = utr1x.d; trd2x = utr2x.d;
                  trd1y = utr1y.d; trd2y = utr2y.d;
                  trd1z = utr1z.d; trd2z = utr2z.d;
                }

                dfix = 2.0*(-sij)*trd1x * (vd / anorm4);
                dfkx = 2.0*( sjk)*trd2x * (vd / anorm4);
                dfjx = -(dfix + dfkx);

                dfiy = 2.0*(-sij)*trd1y * (vd / anorm4);
                dfky = 2.0*( sjk)*trd2y * (vd / anorm4);
                dfjy = -(dfiy + dfky);

                dfiz = 2.0*(-sij)*trd1z * (vd / anorm4);
                dfkz = 2.0*( sjk)*trd2z * (vd / anorm4);
                dfjz = -(dfiz + dfkz);
              }

              if (triplet_debug) /* Compare forces to numerical derivatives */
                force_debug_3v(xx,i,j,k, j,i,k,
                               dfix,dfiy,dfiz,
                               dfjx,dfjy,dfjz,
                               dfkx,dfky,dfkz);


              if (three_body_forces)
                accumulate_forces_3(w3);
            }

            if (T23 != nullptr) {
              //printf("T23 i,j,k = %d,%d,%d\n",i,j,k);
              mcount++;
              if (three_body_energies && evflag) {
                tr2 = transtrace(T23->H1H2,T23->H1H2);
                double dvir = (2.0*(dvir_jk + dvir_ki)*splinepot.vd +
                                 splinepot.dvd)*tr2*w3/anorm4;
                vir3v = vir3v + dvir;
                xvir3v = xvir3v + dvir;
              }

              {
                const double vd = splinepot.vd;

                tr_trace3(&(T23->H1H2.m[1][0]),
                          &(T23->H1xH2.m[1][0]),&utr1x.d,
                          &(T23->H1yH2.m[1][0]),&utr1y.d,
                          &(T23->H1zH2.m[1][0]),&utr1z.d);
                tr_trace3(&(T23->H1H2.m[1][0]),
                          &(T23->H1H2x.m[1][0]),&utr2x.d,
                          &(T23->H1H2y.m[1][0]),&utr2y.d,
                          &(T23->H1H2z.m[1][0]),&utr2z.d);
                if (linalg.single) {
                  trd1x = utr1x.f; trd2x = utr2x.f;
                  trd1y = utr1y.f; trd2y = utr2y.f;
                  trd1z = utr1z.f; trd2z = utr2z.f;
                } else {
                  trd1x = utr1x.d; trd2x = utr2x.d;
                  trd1y = utr1y.d; trd2y = utr2y.d;
                  trd1z = utr1z.d; trd2z = utr2z.d;
                }

                dfix = 2.0*( ski)*trd1x * (vd / anorm4);
                dfjx = 2.0*(-sjk)*trd2x * (vd / anorm4);
                dfkx = -(dfix + dfjx);

                dfiy = 2.0*( ski)*trd1y * (vd / anorm4);
                dfjy = 2.0*(-sjk)*trd2y * (vd / anorm4);
                dfky = -(dfiy + dfjy);

                dfiz = 2.0*( ski)*trd1z * (vd / anorm4);
                dfjz = 2.0*(-sjk)*trd2z * (vd / anorm4);
                dfkz = -(dfiz + dfjz);
              }

              if (triplet_debug) /* Compare forces to numerical derivatives */
                force_debug_3v(xx,i,j,k, k,i,j,
                               dfix,dfiy,dfiz,
                               dfjx,dfjy,dfjz,
                               dfkx,dfky,dfkz);

              if (three_body_forces)
                accumulate_forces_3(w3);

            }

            if (T31 != nullptr) {
              //printf("T31 i,j,k = %d,%d,%d\n",i,j,k);
              mcount++;
              if (three_body_energies && evflag) {
                tr3 = transtrace(T31->H1H2,T31->H1H2);
                double dvir = (2.0*(dvir_ki + dvir_ij)*splinepot.vd +
                                 splinepot.dvd)*tr3*w3/anorm4;
                vir3v = vir3v + dvir;
                xvir3v = xvir3v + dvir;
              }

              {
                const double vd = splinepot.vd;

                tr_trace3(&(T31->H1H2.m[1][0]),
                          &(T31->H1xH2.m[1][0]),&utr1x.d,
                          &(T31->H1yH2.m[1][0]),&utr1y.d,
                          &(T31->H1zH2.m[1][0]),&utr1z.d);
                tr_trace3(&(T31->H1H2.m[1][0]),
                          &(T31->H1H2x.m[1][0]),&utr2x.d,
                          &(T31->H1H2y.m[1][0]),&utr2y.d,
                          &(T31->H1H2z.m[1][0]),&utr2z.d);
                if (linalg.single) {
                  trd1x = utr1x.f; trd2x = utr2x.f;
                  trd1y = utr1y.f; trd2y = utr2y.f;
                  trd1z = utr1z.f; trd2z = utr2z.f;
                } else {
                  trd1x = utr1x.d; trd2x = utr2x.d;
                  trd1y = utr1y.d; trd2y = utr2y.d;
                  trd1z = utr1z.d; trd2z = utr2z.d;
                }

                dfjx = 2.0*( sij)*trd1x * (vd / anorm4);
                dfkx = 2.0*(-ski)*trd2x * (vd / anorm4);
                dfix = -(dfjx + dfkx);

                dfjy = 2.0*( sij)*trd1y * (vd / anorm4);
                dfky = 2.0*(-ski)*trd2y * (vd / anorm4);
                dfiy = -(dfjy + dfky);

                dfjz = 2.0*( sij)*trd1z * (vd / anorm4);
                dfkz = 2.0*(-ski)*trd2z * (vd / anorm4);
                dfiz = -(dfjz + dfkz);

              }

              if (triplet_debug) /* Compare forces to numerical derivatives */
                force_debug_3v(xx,i,j,k, i,j,k,
                               dfix,dfiy,dfiz,
                               dfjx,dfjy,dfjz,
                               dfkx,dfky,dfkz);

              if (three_body_forces)
                accumulate_forces_3(w3);
            }

            v33 = tr0 / anorm3;
            v43 = (tr1 + tr2 + tr3) / anorm4;
            double de_triplet = (splinepot.vc*v33 + splinepot.vd*v43) * e_scale * w3;
            e_triplet = e_triplet + de_triplet;
            e_triplet_c = e_triplet_c + splinepot.vc*v33 * e_scale * w3;
            c_t++;

            //printf("xxxx %6d %6d %6d :: %20.10e\n",1,2,3,de_triplet);

            if (evflag) {
              double drji[3],drki[3];
              double fj[3] = {fjx,fjy,fjz},fk[3] = {fkx,fky,fkz};
              for (int p = 0; p<3; p++) {
                drji[p] = xx[j][p] - xx[i][p];
                drki[p] = xx[k][p] - xx[i][p];
                /* To fix stress-per-atom scaling. */
                fj[p] *= e_scale;
                fk[p] *= e_scale;
              }

              ev_tally3(i,j,k,de_triplet,0.0,fj,fk,drji,drki);

              if (volpres_flag && vflag_atom) {
                //virial[i] = virial[i] - (vir3v + vir3t) * rhoinv*e_scale;
                double dvir = -(xvir3v + xvir3t) * rhoinv*e_scale * (1.0/3.0);
                for (int pp = 0; pp<3; pp++) {
                  vatom[i][pp] += dvir;
                  vatom[j][pp] += dvir;
                  vatom[k][pp] += dvir;
                }
              }

              fix = fix+fsave[0][0]; fiy = fiy+fsave[0][1]; fiz = fiz+fsave[0][2];
              fjx = fjx+fsave[1][0]; fjy = fjy+fsave[1][1]; fjz = fjz+fsave[1][2];
              fkx = fkx+fsave[2][0]; fky = fky+fsave[2][1]; fkz = fkz+fsave[2][2];
            }

            tx1 = gettime();
            ttriplet += tx1 - tx0;
            nttriplet++;
          } else {
            triplet_defer = 1;
          }

          if (four_body_energies || four_body_forces)
            if (j < i) { /* Search for quadruplet */
              tx0 = gettime();

              mj = first[j];
              mk = first[k];
              /*
                i is in both the j-list and the k-list, and i > k,
                and lists are sorted, so the loop terminates.
              */
              while (nlist_short[mj] < i && nlist_short[mk] < i) {

                if (mj >= first[j+1] || mk >= first[k+1]) {
                  printf("Illegal quad...\n"
                         "  j=%d  first[j]=%d  first[j+1]=%d  mj=%d\n"
                         "  k=%d  first[k]=%d  first[k+1]=%d  mk=%d\n",
                         j,first[j],first[j+1],mj,
                         k,first[k],first[k+1],mk);
                  error->one(FLERR,"Shit, brkoen quad loop");
                }

                if (nlist_short[mj] == nlist_short[mk]) {
                  /* Closed quadruplet */
                  m = nlist_short[mj];
                  c_jm = c_km = 1;

                  const int sim = (i < m) ? 1 : -1;
                  const int sjm = (j < m) ? 1 : -1;
                  const int skm = (k < m) ? 1 : -1;

                  w4 = get_weight(triclinic,ss[i],ss[j],ss[k],ss[m]);

                  if (w4 > 0.0) {

                    /* Alrady know ij,jk,ki,jm,km bonds. Look for im bond. */
                    mi = first[i];
                    while (mi < first[i+1] && nlist_short[mi] < m) mi = mi + 1;
                    if (mi < first[i+1] && nlist_short[mi] == m)
                      c_im = 1;
                    else
                      c_im = 0;

                    if (c_im == 0 || c_jk == 0 || (c_jk && c_im && m < k)) {

                      if (triplet_defer) {
                        dvir_ij = dvir_jk = dvir_ki = 0.0;
                        if (c_ij && c_jk)
                          T12 = get_triplet(xx,j,i,k,&bond_hash,&T12work,&dvir_ij,&dvir_jk);
                        if (c_ki && c_jk)
                          T23 = get_triplet(xx,k,i,j,&bond_hash,&T23work,&dvir_ki,&dvir_jk);
                        if (c_ij && c_ki)
                          T31 = get_triplet(xx,i,j,k,&bond_hash,&T31work,&dvir_ij,&dvir_ki);
                        triplet_defer = 0;
                      }


                      fmx = fmy = fmz = 0.0;
                      double xvir4 = 0.0;

                      if (evflag) {
                        fsave[0][0] = fix; fsave[0][1] = fiy; fsave[0][2] = fiz;
                        fsave[1][0] = fjx; fsave[1][1] = fjy; fsave[1][2] = fjz;
                        fsave[2][0] = fkx; fsave[2][1] = fky; fsave[2][2] = fkz;
                        fsave[3][0] = fmx; fsave[3][1] = fmy; fsave[3][2] = fmz;
                        fix = fiy = fiz = 0.0;
                        fjx = fjy = fjz = 0.0;
                        fkx = fky = fkz = 0.0;
                        fmx = fmy = fmz = 0.0;
                      }

                      tr1 = tr2 = tr3 = 0.0;

                      dvir_im = dvir_jm = dvir_km = 0.0;
                      T45 = T56 = T64 = nullptr;
                      if (T12 != nullptr && c_km && c_im)
                        T45 = get_triplet(xx,m,i,k,&bond_hash,&T45work,&dvir_im,&dvir_km);
                      if (T23 != nullptr && c_im && c_jm)
                        T56 = get_triplet(xx,m,i,j,&bond_hash,&T56work,&dvir_im,&dvir_jm);
                      if (T31 != nullptr && c_jm && c_km)
                        T64 = get_triplet(xx,m,j,k,&bond_hash,&T64work,&dvir_jm,&dvir_km);

                      if (T12 != nullptr && T45 != nullptr) {
                        if (four_body_energies && evflag) {
                          tr1 = transtrace(T12->H1H2,T45->H1H2);
                          double dvir = ( (dvir_ij + dvir_jk + dvir_im + dvir_km)*splinepot.ve +
                                          splinepot.dve )*tr1*w4/anorm4;
                          vir4 = vir4 + dvir;
                          xvir4 = xvir4 + dvir;
                        }
                        qcount++;

                        {
                          const double ve = splinepot.ve;

                          trd_update_4(T12,T45);

                          dfix_update_4a(x);
                          dfix_update_4a(y);
                          dfix_update_4a(z);
                        }

                        if (quad_debug) /* Compare forces to numerical derivatives */
                          force_debug_4(xx,i,j,k,m, i,j,k,m,
                                        dfix,dfiy,dfiz , dfjx,dfjy,dfjz,
                                        dfkx,dfky,dfkz , dfmx,dfmy,dfmz);

                        if (four_body_forces)
                          accumulate_forces_4(w4);
                      }

                      if (T23 != nullptr && T56 != nullptr) {
                        if (four_body_energies && evflag) {
                          tr2 = transtrace(T23->H1H2,T56->H1H2);
                          double dvir = ( (dvir_ki + dvir_jk + dvir_im + dvir_jm)*splinepot.ve +
                                          splinepot.dve )*tr2*w4/anorm4;
                          vir4 = vir4 + dvir;
                          xvir4 = xvir4 + dvir;
                        }
                        qcount++;

                        {
                          const double ve = splinepot.ve;

                          trd_update_4(T23,T56);

                          dfix_update_4b(x);
                          dfix_update_4b(y);
                          dfix_update_4b(z);
                        }

                        if (quad_debug) /* Compare forces to numerical derivatives */
                          force_debug_4(xx,i,j,k,m, i,m,j,k,
                                        dfix,dfiy,dfiz , dfjx,dfjy,dfjz,
                                        dfkx,dfky,dfkz , dfmx,dfmy,dfmz);

                        if (four_body_forces)
                          accumulate_forces_4(w4);

                      }

                      if (T31 != nullptr && T64 != nullptr) {
                        if (four_body_energies && evflag) {
                          tr3 = transtrace(T31->H1H2,T64->H1H2);
                          double dvir = ( (dvir_ki + dvir_ij + dvir_jm + dvir_km)*splinepot.ve +
                                          splinepot.dve )*tr3*w4/anorm4;
                          vir4 = vir4 + dvir;
                          xvir4 = xvir4 + dvir;
                        }
                        qcount++;

                        {
                          const double ve = splinepot.ve;

                          /* X */
                          trd_update_4(T31,T64);

                          dfix_update_4c(x);
                          dfix_update_4c(y);
                          dfix_update_4c(z);
                        }

                        if (quad_debug) /* Compare forces to numerical derivatives */
                          force_debug_4(xx,i,j,k,m, i,j,m,k,
                                        dfix,dfiy,dfiz , dfjx,dfjy,dfjz,
                                        dfkx,dfky,dfkz , dfmx,dfmy,dfmz);

                        if (four_body_forces)
                          accumulate_forces_4(w4);
                      }

                      double de_quad = splinepot.ve*(tr1 + tr2 + tr3)/anorm4 * e_scale * w4;
                      e_quad = e_quad + de_quad;
                      if ((T12 && T45) ||
                         (T23 && T56) ||
                         (T31 && T64)) {
                        c_q++;
                      }

                      if (evflag) {
                        double drim[3],drjm[3],drkm[3];
                        double fi[3] = {fix,fiy,fiz};
                        double fj[3] = {fjx,fjy,fjz};
                        double fk[3] = {fkx,fky,fkz};
                        for (int p = 0; p<3; p++) {
                          drim[p] = xx[i][p] - xx[m][p];
                          drjm[p] = xx[j][p] - xx[m][p];
                          drkm[p] = xx[k][p] - xx[m][p];
                          fi[p] *= e_scale;
                          fj[p] *= e_scale;
                          fk[p] *= e_scale;
                        }

                        ev_tally4(i,j,k,m,de_quad,fi,fj,fk,drim,drjm,drkm);

                        if (volpres_flag && vflag_atom) {
                          //virial[i] = virial[i] - vir4 * rhoinv*e_scale;
                          double dvir = -xvir4 * rhoinv*e_scale * (1.0/4.0);
                          for (int pp = 0; pp<3; pp++) {
                            vatom[i][pp] += dvir;
                            vatom[j][pp] += dvir;
                            vatom[k][pp] += dvir;
                            vatom[m][pp] += dvir;
                          }
                        }

                        fix = fix+fsave[0][0]; fiy = fiy+fsave[0][1]; fiz = fiz+fsave[0][2];
                        fjx = fjx+fsave[1][0]; fjy = fjy+fsave[1][1]; fjz = fjz+fsave[1][2];
                        fkx = fkx+fsave[2][0]; fky = fky+fsave[2][1]; fkz = fkz+fsave[2][2];
                        fmx = fmx+fsave[3][0]; fmy = fmy+fsave[3][1]; fmz = fmz+fsave[3][2];
                      }

                      ff[m][0] += fmx * e_scale;
                      ff[m][1] += fmy * e_scale;
                      ff[m][2] += fmz * e_scale;

                    }
                  }
                  mj = mj + 1;
                  mk = mk + 1;
                } else if (nlist_short[mj] < nlist_short[mk]) {
                  mj = mj + 1;
                } else {
                  mk = mk + 1;
                }

              }
              tx1 = gettime();
              tquad += tx1 - tx0;
              ntquad++;
              ntquaditer++;
            }


          ff[k][0] += fkx * e_scale;
          ff[k][1] += fky * e_scale;
          ff[k][2] += fkz * e_scale;

        }
#undef transtrace

        ff[j][0] += fjx * e_scale;
        ff[j][1] += fjy * e_scale;
        ff[j][2] += fjz * e_scale;

      }

    ff[i][0] += fix * e_scale;
    ff[i][1] += fiy * e_scale;
    ff[i][2] += fiz * e_scale;

    if (single_energies == 1 && i < nloc) {
      const double evol0 = splinepot.evol0;
      if (eflag_global) {
        e_single = e_single + evol0 * e_scale;
        eng_vdwl = eng_vdwl + evol0 * e_scale;
      }
      if (eflag_atom) eatom[i] = eatom[i] + evol0 * e_scale;
      if (volpres_flag && vflag_atom) {
        for (int pp = 0; pp<3; pp++)
          vatom[i][pp] = vatom[i][pp] - rhoinv*splinepot.devol0*e_scale;
      }

    }

  }

  tx0 = gettime();
  for (i = 0; i<ntot; i++)
    for (p = 0; p<3; p++)
      atom->f[i][p] = atom->f[i][p] + ff[i][p];

  memory->sfree(nlist_short);
  memory->sfree(first);
  if (ss != xx) memory->sfree(ss);
  memory->sfree(ff);
  memory->sfree(xx);
  tx1 = gettime();
  tmem += tx1-tx0;
  ntmem++;

  t1 = gettime(1);

  //printf("compute_x: c_p = %d    c_t = %d    c_q = %d\n",c_p,c_t,c_q);


#ifdef TIMING_ON
  if (comm->me == 0) {
    double tsum = (tmem+tsort+tpair+tlookup+ttriplet+tquad);
    double nsum = (ntmem+ntsort+ntpair+ntlookup+nttriplet+ntquad);
    //double adj = ((t1-t0)-tsum)/nsum;
    /* Use adj = 6ns for RDTSC, and 58ns for gettimeofday,
       on monkfish.llnl.gov, 2.4GHz Intel

       Use adj = 35.945ns for RDTSC on uBGL (assumed rate set to 700MHz)
    */
    double adj = 35.945e-9;

    double
      memadj     = tmem     - adj*ntmem    ,
      sortadj    = tsort    - adj*ntsort   ,
      pairadj    = tpair    - adj*ntpair   ,
      lookupadj  = tlookup  - adj*ntlookup ,
      tripletadj = ttriplet - adj*nttriplet,
      quadadj    = tquad    - adj*ntquad   ,

      make_b_adj = t_make_b - adj*n_make,
      make_t_adj = t_make_t - adj*n_make,
      make_b2_adj = t_make_b2 - adj*n_make_b2,
      trace_adj = t_trace - adj*n_trace;

    printf("mgpt engy = %10.3fms\n",(t1-t0)*1e3);
    printf("       mem = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           tmem*1e3,ntmem,memadj*1e3,memadj/ntmem*1e9);
    printf("      sort = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           tsort*1e3,ntsort,sortadj*1e3,sortadj/ntsort*1e9);
    printf("      pair = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           tpair*1e3,ntpair,pairadj*1e3,pairadj/ntpair*1e9);
    printf("    lookup = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           tlookup*1e3,ntlookup,lookupadj*1e3,lookupadj/ntlookup*1e9);
    printf("   triplet = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           ttriplet*1e3,nttriplet,tripletadj*1e3,tripletadj/nttriplet*1e9);
    printf("      quad = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           tquad*1e3,ntquaditer,quadadj*1e3,quadadj/ntquaditer*1e9);
    printf("       sum = %10.3fms                adj = %10.3fms\n",
           tsum*1e3,(tsum - adj*nsum)*1e3);
    printf("\n    make_b = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           t_make_b*1e3,n_make,make_b_adj*1e3,make_b_adj/n_make*1e9);
    printf("   make_b2 = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n",
           t_make_b2*1e3,n_make_b2,make_b2_adj*1e3,make_b2_adj/n_make_b2*1e9);
    printf("    make_t = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n\n",
           t_make_t*1e3,n_make,make_t_adj*1e3,make_t_adj/n_make*1e9);
    printf("     trace = %10.3fms  n = %8.0f  adj = %10.3fms  one = %10.3fns\n\n",
           t_trace*1e3,n_trace,trace_adj*1e3,trace_adj/n_trace*1e9);

    printf("mcount (transpose + trace for triplet) = %.0f , %.0f  qcount = %.0f  lmax = %d\n",
           mcount,mcount2,qcount,lmax);

    printf("nbc=%.0f  tbl=%.3fms  tbm=%.3fms  one tbl=%.3fns  one tbm=%.3fns\n",
           nbc,(tbl-adj*nbc)*1e3,(tbm-adj*nbc)*1e3,(tbl/nbc-adj)*1e9,
           (tbm/nbc-adj)*1e9);
    printf("\n\nForces:\n");
    printf("fix = %.3f  fiy=%.3f  fiz=%.3f\n",fix,fiy,fiz);
    printf("fjx = %.3f  fjy=%.3f  fjz=%.3f\n",fjx,fjy,fjz);
    printf("fkx = %.3f  fky=%.3f  fkz=%.3f\n",fkx,fky,fkz);
    printf("\n");

    printf("Bonds   : nsearch=%d maxlen=%d avg.len=%.3f\n",
           bond_hash.NSearch(),bond_hash.MaxLength(),
           bond_hash.NStep()/(double) bond_hash.NSearch());

    printf("compute_x: c_p = %d    c_t = %d    c_q = %d\n",c_p,c_t,c_q);

    printf("@@ Total number of trace3 calls is %d, total number of make_triplet is %.1f\n",
           ntr_calls,n_make);

    {
      Hash<bond_data,Doublet>::Iterator iter = bond_hash.begin();
      int nitem = 0,nhit = 0;
      while (iter != bond_hash.end()) {
        nitem++;
        nhit += iter.link()->hits;
        iter.next();
      }
      printf("bond_hash    hits: nitems=%d   nhits=%d  hits/item = %.3f\n",
             nitem,nhit,nhit/(double) nitem);
    }
  }
#endif

  if (volpres_flag) {
    /*
      Include contributions to the pressure due to derivatines
      of the energy with respect to the potential input volume.
    */

    /* The following lines have moved to beginning of functions,
       since they are used in calculating per-atom virial contributions */
    /*
      double vtot = 1.0;
      double ntot = atom->natoms;
      for (i = 0; i<3; i++)
      vtot = vtot * (domain->boxhi[i] - domain->boxlo[i]);
      double rhoinv = vtot / ntot;
    */

    if (single_energies) // Virial correction for self energy
      for (i = 0; i<3; i++) {
        //virial[i] = virial[i] + nloc*pot_input_vol*pvol0*e_scale;
        virial[i] = virial[i] - nloc*rhoinv*splinepot.devol0*e_scale;
      }

    if (pair_energies) // Virial correction for pair energy
      for (i = 0; i<3; i++)
        virial[i] = virial[i] + rhoinv*e_scale*volvir2;

    if (three_body_energies) // Virial correction for three body enegries
      for (i = 0; i<3; i++) {
        //virial[i] = virial[i] - pot_input_vol*(e_triplet_c*pc + (e_triplet-e_triplet_c)*pd);
        virial[i] = virial[i] - (vir3v + vir3t) * rhoinv*e_scale;
      }

    if (four_body_energies) // Virial correction for four body enegries
      for (i = 0; i<3; i++) {
        //virial[i] = virial[i] - pot_input_vol*e_quad*pe;
        virial[i] = virial[i] - vir4 * rhoinv*e_scale;
      }

  }

  *e_s = e_single;
  *e_p = e_pair;
  *e_t = e_triplet;
  *e_q = e_quad;
}

void PairMGPT::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int newton_pair = force->newton_pair;
  double e_s,e_p,e_t,e_q;

  //printf("newton_pair = %d, newton = %d, tag_enable = %d\n",force->newton_pair,force->newton,atom->tag_enable);

  if (newton_pair == 0) {
    printf("This is a problem. MGPT requires newton_pair flag to be on. Exiting...\n");
    exit(1);
  }

  if (atom->tag_enable == 0) {
    printf("This is a problem. MGPT requires tag_enable flag to be on. Exiting...\n");
    exit(1);
  }

  compute_x(listfull->numneigh,listfull->firstneigh,&e_s,&e_p,&e_t,&e_q,evflag,newton_pair);

  if (false) { // Stupid force calculation / verification
    int ii,nmax=-1;
    for (ii = 0; ii<listfull->inum + listfull->gnum; ii++) {
      int i = listfull->ilist[ii];
      if (i > nmax) nmax = i;
    }
    nmax++;
    double *ffwork = new double[3*nmax];
    double *ffloc = new double[3*listfull->inum];
    double *ffloc2 = new double[3*listfull->inum];
    double **ffptr = new double *[nmax];
    for (ii = 0; ii<listfull->inum + listfull->gnum; ii++)
      ffptr[ii] = &ffwork[3*ii];

    printf("Computing boundary forces\n");
    for (ii = 0; ii<listfull->inum; ii++) {
      ffloc2[3*ii] = 0.0;
      ffloc2[3*ii+1] = 0.0;
      ffloc2[3*ii+2] = 0.0;
      int i = listfull->ilist[ii];
      for (int jj = 0; jj<listfull->inum+listfull->gnum; jj++) {
        int j = listfull->ilist[jj];
        if (atom->tag[i] == atom->tag[j])
          for (int p = 0; p<3; p++)
            ffloc2[3*ii+p] += atom->f[j][p];
      }
    }

    printf("Starting main displacement force calculation\n");
    for (ii = 0; ii<listfull->inum; ii++) {
      int i = listfull->ilist[ii];

      double **atom_f_save = atom->f;
      atom->f = ffptr;

      for (int p = 0; p<3; p++) {
        double xsave = atom->x[i][p];
        const double delta = 1e-3;

        atom->x[i][p] = xsave + delta;
        for (int jj = 0; jj<3*nmax; jj++) ffwork[jj] = 0.0;
        compute_x(listfull->numneigh,
                  listfull->firstneigh,
                  &e_s,&e_p,&e_t,&e_q,evflag,newton_pair);
        double e1 = e_s + e_p + e_t + e_q;

        atom->x[i][p] = xsave - delta;
        for (int jj = 0; jj<3*nmax; jj++) ffwork[jj] = 0.0;
        compute_x(listfull->numneigh,
                  listfull->firstneigh,
                  &e_s,&e_p,&e_t,&e_q,evflag,newton_pair);
        double e2 = e_s + e_p + e_t + e_q;

        ffloc[3*ii+p] = -(e1-e2)/(2*delta);

        atom->x[i][p] = xsave;
      }

      atom->f = atom_f_save;
      printf("Force on i=%4d:\n",i);
      printf("  Position   %20.10e  %20.10e  %20.10e\n",
             atom->x[i][0],atom->x[i][1],atom->x[i][2]);
      printf("  Exact      %20.10e  %20.10e  %20.10e\n",
             atom->f[i][0],atom->f[i][1],atom->f[i][2]);
      printf("  Numerical  %20.10e  %20.10e  %20.10e\n",
             ffloc[3*ii+0],ffloc[3*ii+1],ffloc[3*ii+2]);
      printf("  Boundary   %20.10e  %20.10e  %20.10e\n",
             ffloc2[3*ii+0],ffloc2[3*ii+1],ffloc2[3*ii+2]);
    }


    delete[] ffloc2;
    delete[] ffloc;
    delete[] ffptr;
    delete[] ffwork;
  }


  if (false) {
    printf("\nForces MGPT:\n");
    const int iimax = (listfull->inum < 10) ? listfull->inum : 10;
    for (int ii = 0; ii<iimax; ii++) {
      int i = listfull->ilist[ii];
      printf("%4d  =  %20.10e  %20.10e  %20.10e\n",
             i,atom->f[i][0],atom->f[i][1],atom->f[i][2]);
    }
    printf("\n\n");
  }

  if (vflag_fdotr) {
    //printf("##### Using virial_compute!!!\n");
    virial_fdotr_compute();
  }
}

void PairMGPT::allocate()
{
        allocated = 1;
        int n = atom->ntypes;

        memory->create(setflag,n+1,n+1,"pair:setflag");
        for (int i = 0; i <= n; i++)
                for (int j = 0; j <= n; j++)
                        setflag[i][j] = 0;

        memory->create(cutsq,n+1,n+1,"pair:cutsq");
        memory->create(cutghost,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairMGPT::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMGPT::coeff(int narg, char **arg)
{
  int single_precision = 0;

  if (narg < 5)
    error->all(FLERR,
               "Not enough arguments for mgpt (MGPT) pair coefficients.");

  if (!allocated) allocate();

  // Make sure I,J args are * *
  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  double vol;
  if (sscanf(arg[4], "%lg", &vol) != 1 || vol <= 0.0)
    error->all(FLERR,"Invalid volume in mgpt (MGPT) pair coefficients.");

  volpres_flag = 1;
  single_precision = 0;

  /* Parse arguments */ {
    int volpres_tag = 0,precision_tag = 0,nbody_tag = 0;

    int iarg = 5;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"volpress") == 0) { /* Volumetric pressure flag */
        if (iarg+2 > narg)
          error->all(FLERR,"Incorrect args for pair coefficients");
        if (strcmp(arg[iarg+1],"yes") == 0) volpres_flag = 1;
        else if (strcmp(arg[iarg+1],"no") == 0) volpres_flag = 0;
        else {
          char line[1024];
          sprintf(line,"(In %s:%d) Invalid value for volumetric pressure argument.\n"
                  "It should be \"volpress yes\" or \"volpress no\".\n"
                  "The value is \"%s\".\n",FLERR,arg[iarg+1]);
          error->all(FLERR,line);
        }
        volpres_tag = 1;
        iarg += 2;
        if (comm->me == 0) printf("* volpress: volpres_flag = %d [%s %s]\n",volpres_flag,arg[iarg-2],arg[iarg-1]);
      } else if (strcmp(arg[iarg],"nbody") == 0) {
        if (iarg+2 > narg)
          error->all(FLERR,"Incorrect args for pair coefficients");
        if (strspn(arg[iarg+1],"1234") == strlen(arg[iarg+1])) {
          nbody_flag = 0;
          for (int i = 0; i<4; i++)
            if (strchr(arg[iarg+1],'1'+i) != nullptr) {
              nbody_flag = nbody_flag + (1<<i);
              if (comm->me == 0) printf("Explicitly adding %d-tuple forces.\n",i+1);
            }
        } else {
          char line[1024];
          sprintf(line,"(In %s:%d) Invalid value for nbody flag.\n"
                  "It should be e.g. \"nbody=1234\" (for single, pair, triple, and quad forces/energiers)\n"
                  "For e.g. only pair and triple forces/energies, use \"nbody=23\".\n"
                  "The default is \"nbody=1234\".\n"
                  "The current value is \"%s\".\n",FLERR,arg[iarg+1]);
          error->all(FLERR,line);
        }
        nbody_tag = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg],"precision") == 0) {
        if (iarg+2 > narg)
          error->all(FLERR,"Incorrect args for pair coefficients");
        if (strcmp(arg[iarg+1],"single") == 0) single_precision = 1;
        else if (strcmp(arg[iarg+1],"double") == 0) single_precision = 0;
        else {
          char line[1024];
          sprintf(line,"(In %s:%d) Invalid value for precision argument.\n"
                  "It should be \"precision single\" or \"precision double\".\n"
                  "The value is \"%s\".\n",FLERR,arg[iarg+1]);
          error->all(FLERR,line);
        }
        precision_tag = 1;
        iarg += 2;
        if (comm->me == 0) printf("* precision: single_flag = %d [%s %s]\n",single_precision,arg[iarg-2],arg[iarg-1]);
      } else {
        char line[1024];
        sprintf(line,"(In %s:%d) Invalid argument. Allowed arguments are:\n"
                "    volpress {yes|no} , default = yes\n"
                "    precision {single|double} , default = double\n"
                "    nbody {[1234,]*} , default = whichever terms potential require\n"
                "The invalid argument is \"%s\".\n",FLERR,arg[iarg]);
        error->all(FLERR,line);
      }
    }

    if (comm->me == 0)
      printf("Volumetric pressure is %s.\n",volpres_flag ? "on" : "off");

    if (comm->me == 0) {
      FILE *parmin_fp = utils::open_potential(arg[2],lmp,nullptr);
      FILE *potin_fp = utils::open_potential(arg[3],lmp,nullptr);
      if (parmin_fp == nullptr || potin_fp == nullptr) {
        char str[128];
        sprintf(str,"Cannot open MGPT potential files %s %s",arg[2],arg[3]);
        error->one(FLERR,str);
      }
      fclose(parmin_fp);
      fclose(potin_fp);

      splinepot.readpot(arg[2],arg[3],vol);
      printf("evol0 = %.10e\n",splinepot.evol0);

      /* Set up default and requested nbody forces to include */ {
        int nbody_default = (1<<0) + (1<<1) + (1<<2) + (1<<3);

        if (splinepot.vd == 0.0 && splinepot.dvd == 0.0)
          nbody_default -= (1<<2); // No 3-body contributions
        if (splinepot.ve == 0.0 && splinepot.dve == 0.0)
          nbody_default -= (1<<3); // No 4-body contributions

        if (nbody_tag == 0) nbody_flag = nbody_default;

        if (nbody_flag != nbody_default) {
          printf("Warning: nbody=%d (suggested=%d) set to disregard multibody-forces in potential.\n",
                 nbody_flag,nbody_default);
        }
      }
    }
  }

  MPI_Bcast(&nbody_flag,sizeof(nbody_flag),MPI_BYTE,0,world);

  /*
    Broadcast structure to all processes. In receiving
    processes, pointes will be screwed up. We allocate
    memory, and then broadcast contents of arrays.
  */
  MPI_Bcast(&splinepot,sizeof(splinepot),MPI_BYTE,0,world);
  if (comm->me != 0) {
    splinepot.vpair_spline = new double[splinepot.nr-1][4];
    splinepot.dvpair_spline = new double[splinepot.nr-1][4];
  }
  MPI_Bcast(splinepot.vpair_spline,4*(splinepot.nr-1),MPI_DOUBLE,0,world);
  MPI_Bcast(splinepot.dvpair_spline,4*(splinepot.nr-1),MPI_DOUBLE,0,world);
  anorm3 = splinepot.anorm3;
  anorm4 = splinepot.anorm4;
  lmax = splinepot.lmax;
  lang = splinepot.lang;
  //ipot = splinepot.ipot;
  for (int i = 0; i<(int) (sizeof(ddl)/sizeof(double)); i++)
    ddl[i] = splinepot.ddl[i];
  for (int i = 0; i<lmax; i++) {
    for (int j = 0; j<lmax; j++)
      del0.m[i+1][j+1] = 0.0;
    del0.m[i+1][i+1] = 1.0;
  }

  /* Set matrix param, cutoff, LAMMPS param */
  Matrix::sz = lmax;

  rcrit = splinepot.rcrit;
  rmax = splinepot.rmax;
  cutoff = rmax;
  if (rcrit > rmax) cutoff = rcrit;

  // Set LAMMPS pair interaction flags.
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 1; j <= atom->ntypes; j++) {
      setflag[i][j] = 1;
      cutsq[i][j] = cutoff;
      cutghost[i][j] = cutoff;
    }
  }

  // Set atomic mass.
  for (int i = 1; i <= atom->ntypes; i++)
    atom->set_mass(FLERR,i, splinepot.mass);

  // Initialize linear algebra routines.
  linalg = mgpt_linalg(lmax,single_precision);
  if (comm->me == 0)
    printf("%s",linalg.msg);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairMGPT::init_style()
{
        if (force->newton_pair == 0)
          error->all(FLERR,"Pair style mgpt requires newton pair on.");

        // Need a half list and a full neighbor list with neighbors of ghosts
        neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST)->set_id(1);
        neighbor->add_request(this)->set_id(2);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */
void PairMGPT::init_list(int id, NeighList *ptr)
{
        if (id == 1) listfull = ptr;
        else if (id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairMGPT::init_one(int /*i*/, int /*j*/)
{
        return cutoff;
}

/************************************************************************
 ****  REIMPLEMENTATION OF FL AND HAMLTN WITH ANALYTICAL DERIVATIVES ****
 ************************************************************************/
/*
  Reimplementation of bond length potential, including
  derivatives with respect to x,y, and z.
*/
void PairMGPT::fl_deriv_new(double r,double ri,double xhat,double yhat,double zhat,
                            double &fl_0,double &fl_x,double &fl_y,double &fl_z,
                            double &fl_rp,double &fl_p1,double &fl_r0,double &fl_al) {
  const double rp = splinepot.rp,p1 = splinepot.p1,r0 = splinepot.r00,al = splinepot.al;
  const int mode = splinepot.mode;
  const double pn = splinepot.pn;

  double t,tx,ty,tz,t_rp_ti,t_p1_ti;
  double s;

  /*
  // Original code
    double term;
    double pn=1.0;
    if (mode <= 4)
      term = pow(rp/r, p1);
    else
      term = exp(-p1*(pow(r/rp, pn) - 1.0)/pn);
  */

  double rpi = 1.0/rp;
  if (mode <= 4) {
    t = pow(rp*ri,p1);
    s = -p1 * t * ri;
    t_rp_ti = p1*rpi;
    t_p1_ti = log(rp*ri);
  } else {
    if (pn == 1.0) {
      double p1_rpi = -p1*rpi;
      t = exp(p1 + r*p1_rpi);
      s = p1_rpi * t;
      t_rp_ti = -r*p1_rpi*rpi;
      t_p1_ti = 1.0 - r*rpi;
    } else {
      double pni = 1.0/pn;
      double rprpn = pow(r*rpi,pn);
      t = exp(-p1*pni*(rprpn - 1.0));
      s = -p1*rprpn*ri * t;
      t_rp_ti = p1*rprpn*rpi;
      t_p1_ti = pni - pni*rprpn;// -pni*(rprpn - 1.0);
    }
  }
  tx = s * xhat;
  ty = s * yhat;
  tz = s * zhat;

  fl_rp = t_rp_ti;
  fl_p1 = t_p1_ti;

  if (r <= r0) {
    fl_0 = t;
    fl_x = tx;
    fl_y = ty;
    fl_z = tz;

    fl_r0 = 0.0;
    fl_al = 0.0;
  } else {
    double q,qx,qy,qz,exp_q,q_r0,q_al;
    double r0i,u;

    r0i = 1.0/r0;
    u = r*r0i - 1.0;
    q = al*u*u;
    s = 2*al*u*r0i;
    qx = s * xhat;
    qy = s * yhat;
    qz = s * zhat;
    q_r0 = -2.0*al*u*r*r0i*r0i;
    q_al = u*u;

    exp_q = exp(-q);
    if (mode <= 2) {
      fl_0 = exp_q * t;
      fl_x = exp_q*(tx - t*qx);
      fl_y = exp_q*(ty - t*qy);
      fl_z = exp_q*(tz - t*qz);

      fl_r0 = -q_r0;
      fl_al = -q_al;
    } else {
      fl_0 = exp_q * (1.0 + q) * t;
      fl_x = exp_q * (tx + q*(tx - t*qx));
      fl_y = exp_q * (ty + q*(ty - t*qy));
      fl_z = exp_q * (tz + q*(tz - t*qz));

      fl_r0 = -q_r0 * q/(1.0 + q);
      fl_al = -q_al * q/(1.0 + q);
    }
  }
}


/*
  Macros to build elements of the bond matrix, and also
  its derivatives with repsect to x,y, and z.
*/
#define MAKE_ELEMENT_5(i,j)                                          \
  do {                                                               \
    const double dl0  = del0.m[i][j];                                \
    const double dl4  = gsl_##i    * gsl_##j;                        \
    const double dl4x = gsl_##i##x * gsl_##j + gsl_##i * gsl_##j##x; \
    const double dl4y = gsl_##i##y * gsl_##j + gsl_##i * gsl_##j##y; \
    const double dl4z = gsl_##i##z * gsl_##j + gsl_##i * gsl_##j##z; \
                                                                     \
    const double tmp  = w4*dl4  + w2*dl2 + w0*dl0;                   \
    const double tmpx = w4*dl4x + w2*dl2x;                           \
    const double tmpy = w4*dl4y + w2*dl2y;                           \
    const double tmpz = w4*dl4z + w2*dl2z;                           \
    const double tmpsum = tmpx*x + tmpy*y + tmpz*z;                  \
    M [j][i] = M[i][j] =  fl  *tmp;                                  \
    Mx[j][i] = Mx[i][j] = fl_x*tmp + fl_ri*(tmpx - x*tmpsum);        \
    My[j][i] = My[i][j] = fl_y*tmp + fl_ri*(tmpy - y*tmpsum);        \
    Mz[j][i] = Mz[i][j] = fl_z*tmp + fl_ri*(tmpz - z*tmpsum);        \
  } while (0)

#define MAKE_ELEMENT_7(i,j)                                          \
  do {                                                               \
    const double dl0  = del0.m[i][j];                                \
    const double dl6  = gsl_##i    * gsl_##j;                        \
    const double dl6x = gsl_##i##x * gsl_##j + gsl_##i * gsl_##j##x; \
    const double dl6y = gsl_##i##y * gsl_##j + gsl_##i * gsl_##j##y; \
    const double dl6z = gsl_##i##z * gsl_##j + gsl_##i * gsl_##j##z; \
                                                                     \
    const double tmp  = w6*dl6  + w4*dl4  + w2*dl2 + w0*dl0;         \
    const double tmpx = w6*dl6x + w4*dl4x + w2*dl2x;                 \
    const double tmpy = w6*dl6y + w4*dl4y + w2*dl2y;                 \
    const double tmpz = w6*dl6z + w4*dl4z + w2*dl2z;                 \
    const double tmpsum = tmpx*x + tmpy*y + tmpz*z;                  \
    M [j][i] = M[i][j] =  fl  *tmp;                                  \
    Mx[j][i] = Mx[i][j] = fl_x*tmp + fl_ri*(tmpx - x*tmpsum);        \
    My[j][i] = My[i][j] = fl_y*tmp + fl_ri*(tmpy - y*tmpsum);        \
    Mz[j][i] = Mz[i][j] = fl_z*tmp + fl_ri*(tmpz - z*tmpsum);        \
  } while (0)
/* End of bond matrix macros */


/*
  Construction of bond matrix, and its derivatives
  with respect to the coordinates
*/
void PairMGPT::hamltn_5_raw(const double xin,const double yin,const double zin,
                             double M [8][8],double Mx[8][8],
                             double My[8][8],double Mz[8][8],
                             double *fl_deriv_sum_p) {

  const double r = sqrt(xin*xin + yin*yin + zin*zin),ri = 1.0/r;
  const double x = xin*ri,y = yin*ri,z = zin*ri;

  // d-d
  // call delndd(x,y,z)
  const double x2 = x*x,y2 = y*y,z2 = z*z;
  const double xy = x*y,xz = x*z,yz = y*z;

  const double sr3 = sqrt(3.0),sr3i = 1.0/sr3;
  const double frac_1_3 = 1.0/3.0,frac_2_3 = 2.0/3.0,frac_4_3 = 4.0/3.0;

  const double ddl_1 = ddl[1],ddl_2 = ddl[2],ddl_3 = ddl[3];
  const double w4 = ddl_1 - frac_4_3*ddl_2 + frac_1_3*ddl_3;
  const double w2 = ddl_2 - ddl_3;
  const double w0 = ddl_2;


  //del4
  double gsl_1 ,gsl_2 ,gsl_3 ,gsl_4 ,gsl_5;
  double gsl_1x,gsl_2x,gsl_3x,gsl_4x,gsl_5x;
  double gsl_1y,gsl_2y,gsl_3y,gsl_4y,gsl_5y;
  double gsl_1z,gsl_2z,gsl_3z,gsl_4z,gsl_5z;

  double dl2,dl2x,dl2y,dl2z;

  double fl,fl_x,fl_y,fl_z,fl_ri;
  double fl_rp,fl_p1,fl_r0,fl_al;

  gsl_1  = 0.5*(3.0*z2 - 1.0);
  gsl_1x = 0.0;
  gsl_1y = 0.0;
  gsl_1z = 3.0*z;

  gsl_2  = sr3*xz;
  gsl_2x = sr3*z;
  gsl_2y = 0.0;
  gsl_2z = sr3*x;

  gsl_3  = sr3*yz;
  gsl_3x = 0.0;
  gsl_3y = sr3*z;
  gsl_3z = sr3*y;

  gsl_4  = sr3*(x2 - y2)*0.5;
  gsl_4x = sr3*x;
  gsl_4y = -sr3*y;
  gsl_4z = 0.0;

  gsl_5  = sr3*xy;
  gsl_5x = sr3*y;
  gsl_5y = sr3*x;
  gsl_5z = 0.0;

  // Compute bond length potential
  fl_deriv_new(r,ri,x,y,z,fl,fl_x,fl_y,fl_z , fl_rp,fl_p1,fl_r0,fl_al);
  fl_ri = fl*ri;

  *fl_deriv_sum_p =
    fl_rp*splinepot.drp  + fl_p1*splinepot.dp1 +
    fl_r0*splinepot.dr00 + fl_al*splinepot.dal;

  // del2

  //del2.m[1][1] = z2 - 2.0/3.0;
  dl2 = z2 - frac_2_3;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 2*z;
  MAKE_ELEMENT_5(1,1);

  //del2.m[1][2] = xz/sr3;
  dl2  = xz*sr3i;
  dl2x = z*sr3i;
  dl2y = 0.0;
  dl2z = x*sr3i;
  MAKE_ELEMENT_5(1,2);

  //del2.m[1][3] = yz/sr3;
  dl2  = yz*sr3i;
  dl2x = 0.0;
  dl2y = z*sr3i;
  dl2z = y*sr3i;
  MAKE_ELEMENT_5(1,3);

  //del2.m[1][4] = -(x2 - y2)*sr3i;
  dl2  = -(x2 - y2)*sr3i;
  dl2x = -2.0*sr3i*x;
  dl2y =  2.0*sr3i*y;
  dl2z = 0.0;
  MAKE_ELEMENT_5(1,4);

  //del2.m[1][5] = -2.0*xy*sr3i;
  dl2 = -2.0*xy*sr3i;
  dl2x = -2.0*y*sr3i;
  dl2y = -2.0*x*sr3i;
  dl2z =  0.0;
  MAKE_ELEMENT_5(1,5);

  //del2.m[2][2] = -y2;
  dl2  = -y2;
  dl2x =  0.0;
  dl2y = -2.0*y;
  dl2z =  0.0;
  MAKE_ELEMENT_5(2,2);

  //del2.m[2][3] = xy;
  dl2  = xy;
  dl2x = y;
  dl2y = x;
  dl2z = 0.0;
  MAKE_ELEMENT_5(2,3);

  //del2.m[2][4] = xz;
  dl2  = xz;
  dl2x = z;
  dl2y = 0.0;
  dl2z = x;
  MAKE_ELEMENT_5(2,4);

  //del2.m[2][5] = yz;
  dl2  = yz;
  dl2x = 0.0;
  dl2y = z;
  dl2z = y;
  MAKE_ELEMENT_5(2,5);

  //del2.m[3][3] = -x2;
  dl2  = -x2;
  dl2x = -2.0*x;
  dl2y =  0.0;
  dl2z =  0.0;
  MAKE_ELEMENT_5(3,3);

  //del2.m[3][4] = -yz;
  dl2  = -yz;
  dl2x =  0.0;
  dl2y = -z;
  dl2z = -y;
  MAKE_ELEMENT_5(3,4);

  //del2.m[3][5] = xz;
  dl2  = xz;
  dl2x = z;
  dl2y = 0.0;
  dl2z = x;
  MAKE_ELEMENT_5(3,5);

  //del2.m[4][4] = -z2;
  dl2  = -z2;
  dl2x =  0.0;
  dl2y =  0.0;
  dl2z = -2.0*z;
  MAKE_ELEMENT_5(4,4);

  //del2.m[4][5] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;
  MAKE_ELEMENT_5(4,5);

  //del2.m[5][5] = -z2;
  dl2  = -z2;
  dl2x =  0.0;
  dl2y =  0.0;
  dl2z = -2.0*z;
  MAKE_ELEMENT_5(5,5);
}


void PairMGPT::hamltn_7_raw(const double xin,const double yin,const double zin,
                  double M [8][8],double Mx[8][8],
                  double My[8][8],double Mz[8][8],
                             double *fl_deriv_sum_p) {

  const double r = sqrt(xin*xin + yin*yin + zin*zin),ri = 1.0/r;
  const double x = xin*ri,y = yin*ri,z = zin*ri;

  // d-d
  // call delndd(x,y,z)
  const double x2 = x*x,y2 = y*y,z2 = z*z;
  const double xy = x*y,xz = x*z,yz = y*z;
  const double x4 = x2*x2,y4 = y2*y2;

  //const double sr3 = sqrt(3.0);//,sr3i = 1.0/sr3;
  //const double frac_1_3 = 1.0/3.0,frac_2_3 = 2.0/3.0,frac_4_3 = 4.0/3.0;

  const double sr01   = sqrt(0.1);
  const double sr015  = sqrt(0.15);
  const double sr024  = sqrt(0.24);
  const double sr0375 = sqrt(0.375);
  const double sr06   = sqrt(0.6);
  const double sr0625 = sqrt(0.625);
  const double sr09   = sqrt(0.9);
  const double sr15   = sqrt(1.5);
  const double sr24   = sqrt(2.4);
  const double sr36   = sqrt(3.6);
  const double sr375  = sqrt(3.75);
  const double sr96   = sqrt(9.6);
  const double sr150  = sqrt(15.0);



  const double ddl_1 = ddl[1],ddl_2 = ddl[2],ddl_3 = ddl[3],ddl_4 = ddl[4];
  const double w6 = ddl_1 - 1.5*ddl_2 + 0.6*ddl_3 - 0.1*ddl_4;
  const double w4 = 0.625*ddl_2 - ddl_3 + 0.375*ddl_4;
  const double w2 = 0.625*(ddl_2 - ddl_4);
  const double w0 = 0.625*ddl_2 + 0.375*ddl_4;


  //del6
  double gsl_1 ,gsl_2 ,gsl_3 ,gsl_4 ,gsl_5 ,gsl_6, gsl_7;
  double gsl_1x,gsl_2x,gsl_3x,gsl_4x,gsl_5x,gsl_6x,gsl_7x;
  double gsl_1y,gsl_2y,gsl_3y,gsl_4y,gsl_5y,gsl_6y,gsl_7y;
  double gsl_1z,gsl_2z,gsl_3z,gsl_4z,gsl_5z,gsl_6z,gsl_7z;

  double dl2,dl2x,dl2y,dl2z;
  double dl4,dl4x,dl4y,dl4z;
  double t1;

  double fl,fl_x,fl_y,fl_z,fl_ri;
  double fl_rp,fl_p1,fl_r0,fl_al;

  //gslf[1] = 0.5*(5.0*n2 - 3.0)*n;
  gsl_1  = 0.5*(5.0*z2 - 3.0)*z;
  gsl_1x = 0.0;
  gsl_1y = 0.0;
  gsl_1z = 7.5*z2 - 1.5;

  //gslf[2] = sr0375*(5.0*n2 - 1.0)*l;
  gsl_2  = sr0375*(5.0*z2 - 1.0)*x;
  gsl_2x = sr0375*(5.0*z2 - 1.0);
  gsl_2y = 0.0;
  gsl_2z = sr0375*10.0*xz;

  //gslf[3] = sr0375*(5.0*n2 - 1.0)*m;
  gsl_3  = sr0375*(5.0*z2 - 1.0)*y;
  gsl_3x = 0.0;
  gsl_3y = sr0375*(5.0*z2 - 1.0);
  gsl_3z = sr0375*10.0*yz;

  //gslf[4] = sr375*(l2 - m2)*n;
  gsl_4  =  sr375*(x2 - y2)*z;
  gsl_4x =  2.0*sr375*xz;
  gsl_4y = -2.0*sr375*yz;
  gsl_4z =  sr375*(x2 - y2);

  //gslf[5] = sr150*lm*n;
  gsl_5  = sr150*xy*z;
  gsl_5x = sr150*yz;
  gsl_5y = sr150*xz;
  gsl_5z = sr150*xy;

  //gslf[6] = sr0625*(l2 - 3.0*m2)*l;
  gsl_6  =  sr0625*(x2 - 3.0*y2)*x;
  gsl_6x =  3.0*sr0625*(x2 - y2);
  gsl_6y = -6.0*sr0625*xy;
  gsl_6z =  0.0;

  //gslf[7] = sr0625*(3.0*l2 - m2)*m;
  gsl_7  = sr0625*(3.0*x2 - y2)*y;
  gsl_7x = 6.0*sr0625*xy;
  gsl_7y = 3.0*sr0625*(x2 - y2);
  gsl_7z = 0.0;


  // Compute bond length potential
  fl_deriv_new(r,ri,x,y,z,fl,fl_x,fl_y,fl_z , fl_rp,fl_p1,fl_r0,fl_al);
  fl_ri = fl*ri;

  *fl_deriv_sum_p =
    fl_rp*splinepot.drp  + fl_p1*splinepot.dp1 +
    fl_r0*splinepot.dr00 + fl_al*splinepot.dal;

  // del2f

  //del2f.m[1][1] = 0.4*(3.0*n2 - 1.0);
  dl2  = 0.4*(3.0*z2 - 1.0);
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 2.4*z;

  //del4f.m[1][1] = 0.60*(5.0*n2 - 4.0)*n2;
  dl4  = 0.60*(5.0*z2 - 4.0)*z2;
  dl4x = 0.0;
  dl4y = 0.0;
  dl4z = 0.60*(20.0*z2 - 8.0)*z;

  MAKE_ELEMENT_7(1,1);


  //del2f.m[1][2] = sr024*ln;
  dl2  = sr024*xz;
  dl2x = sr024*z;
  dl2y = 0.0;
  dl2z = sr024*x;

  //del4f.m[1][2] = sr024*(5.0*n2 - 2.0)*ln;
  dl4  = sr024*(5.0*z2 - 2.0)*xz;
  dl4x = sr024*(5.0*z2 - 2.0)*z;
  dl4y = 0.0;
  dl4z = sr024*(15.0*z2 - 2.0)*x;

  MAKE_ELEMENT_7(1,2);


  //del2f.m[1][3] = sr024*mn;
  dl2  = sr024*yz;
  dl2x = 0.0;
  dl2y = sr024*z;
  dl2z = sr024*y;

  //del4f.m[1][3] = sr024*(5.0*n2 - 2.0)*mn;
  dl4  = sr024*(5.0*z2 - 2.0)*yz;
  dl4x = 0.0;
  dl4y = sr024*(5.0*z2 - 2.0)*z;
  dl4z = sr024*(15.0*z2 - 2.0)*y;

  MAKE_ELEMENT_7(1,3);


  //del2f.m[1][4] = -sr06*(l2 - m2);
  dl2  = -sr06*(x2 - y2);
  dl2x = -2.0*sr06*x;
  dl2y =  2.0*sr06*y;
  dl2z = 0.0;

  //del4f.m[1][4] = -sr06*(l2 - m2)*n2;
  dl4  = -sr06*(x2 - y2)*z2;
  dl4x = -2.0*sr06*x*z2;
  dl4y =  2.0*sr06*y*z2;
  dl4z = -2.0*sr06*(x2 - y2)*z;

  MAKE_ELEMENT_7(1,4);


  //del2f.m[1][5] = -sr24*lm;
  dl2  = -sr24*xy;
  dl2x = -sr24*y;
  dl2y = -sr24*x;
  dl2z =  0.0;

  //del4f.m[1][5] = -sr24*lm*n2;
  dl4  = -sr24*xy*z2;
  dl4x = -sr24*y*z2;
  dl4y = -sr24*x*z2;
  dl4z = -2.0*sr24*xy*z;

  MAKE_ELEMENT_7(1,5);


  //del2f.m[1][6] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[1][6] = sr36*(3.0*m2 - l2)*ln;
  dl4  = sr36*(3.0*y2 - x2)*xz;
  dl4x = 3.0*sr36*(y2 - x2)*z;
  dl4y = 6.0*sr36*y*xz;
  dl4z = sr36*(3.0*y2 - x2)*x;

  MAKE_ELEMENT_7(1,6);


  //del2f.m[1][7] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[1][7] = -sr36*(3.0*l2 - m2)*mn;
  dl4  = -sr36*(3.0*x2 - y2)*yz;
  dl4x = -6.0*sr36*x*yz;
  dl4y = -3.0*sr36*(x2 - y2)*z;
  dl4z = -sr36*(3.0*x2 - y2)*y;

  MAKE_ELEMENT_7(1,7);


  //del2f.m[2][2] = 0.3*(1.0 - 4.0*m2 + n2);
  dl2  =  0.3 - 1.2*y2 + 0.3*z2;
  dl2x =  0.0;
  dl2y = -2.4*y;
  dl2z =  0.6*z;

  //del4f.m[2][2] = -0.4*l*l - 2.5*(m2 - 0.6*l2)*n2;
  dl4  = -0.4*x2 - 2.5*(y2 - 0.6*x2)*z2;
  dl4x = -0.8*x + 3.0*x*z2;
  dl4y = -5.0*y*z2;
  dl4z = -5.0*(y2 - 0.6*x2)*z;

  MAKE_ELEMENT_7(2,2);


  //del2f.m[2][3] = 1.2*lm;
  dl2  = 1.2*xy;
  dl2x = 1.2*y;
  dl2y = 1.2*x;
  dl2z = 0.0;

  //del4f.m[2][3] = 0.4*(10.0*n2 - 1.0)*lm;
  dl4  = (4.0*z2 - 0.4)*xy;
  dl4x = (4.0*z2 - 0.4)*y;
  dl4y = (4.0*z2 - 0.4)*x;
  dl4z = 8.0*z*xy;

  MAKE_ELEMENT_7(2,3);


  //del2f.m[2][4] = sr09*ln;
  dl2  = sr09*xz;
  dl2x = sr09*z;
  dl2y = 0.0;
  dl2z = sr09*x;

  //del4f.m[2][4] = sr01*(6.0*n2 - 8.0*m2 - 1.0)*ln;
  dl4  =  sr01*(6.0*z2 - 8.0*y2 - 1.0)*xz;
  dl4x =  sr01*(6.0*z2 - 8.0*y2 - 1.0)*z;
  dl4y = -16.0*sr01*y*xz;
  dl4z =  sr01*(18.0*z2 - 8.0*y2 - 1.0)*x;

  MAKE_ELEMENT_7(2,4);


  //del2f.m[2][5] = sr09*mn;
  dl2  = sr09*yz;
  dl2x = 0.0;
  dl2y = sr09*z;
  dl2z = sr09*y;

  //del4f.m[2][5] = sr01*(2.0*n2 - 8.0*m2 + 3.0)*mn;
  dl4  = sr01*(2.0*z2 - 8.0*y2 + 3.0)*yz;
  dl4x = 0.0;
  dl4y = sr01*(2.0*z2 - 24.0*y2 + 3.0)*z;
  dl4z = sr01*(6.0*z2 - 8.0*y2 + 3.0)*y;

  MAKE_ELEMENT_7(2,5);


  //del2f.m[2][6] = -sr015*(l2 - m2);
  dl2  = -sr015*(x2 - y2);
  dl2x = -2.0*sr015*x;
  dl2y =  2.0*sr015*y;
  dl2z =  0.0;

  //del4f.m[2][6] = sr375*(l2 - m2 - 1.4*l4 + 1.2*l2*m2 + m4);
  dl4  = sr375*(x2 - y2 - 1.4*x4 + 1.2*x2*y2 + y4);
  dl4x = sr375*(2.0 - 5.6*x2 + 2.4*y2)*x;
  dl4y = sr375*(-2.0 + 2.4*x2 + 4.0*y2)*y;
  dl4z = 0.0;

  MAKE_ELEMENT_7(2,6);


  //del2f.m[2][7] = -sr06*lm;
  dl2  = -sr06*xy;
  dl2x = -sr06*y;
  dl2y = -sr06*x;
  dl2z =  0.0;

  //del4f.m[2][7] = sr96*(n2 - l2 + 0.25)*lm;
  dl4  = sr96*(z2 - x2 + 0.25)*xy;
  dl4x = sr96*(z2 - 3.0*x2 + 0.25)*y;
  dl4y = sr96*(z2 - x2 + 0.25)*x;
  dl4z = 2.0*sr96*z*xy;

  MAKE_ELEMENT_7(2,7);


  //del2f.m[3][3] = 0.30*(1.0 - 4.0*l2 + n2);
  dl2  =  0.3 - 1.2*x2 + 0.3*z2;
  dl2x = -2.4*x;
  dl2y =  0.0;
  dl2z =  0.6*z;

  //del4f.m[3][3] = -0.4*m2 - 2.5*(l2 - 0.6*m2)*n2;
  dl4  = -0.4*y2 - 2.5*(x2 - 0.6*y2)*z2;
  dl4x = -5.0*x*z2;
  dl4y = y*(3.0*z2 - 0.8);
  dl4z = -5.0*(x2 - 0.6*y2)*z;

  MAKE_ELEMENT_7(3,3);


  //del2f.m[3][4] = -sr09*mn;
  dl2  = -sr09*yz;
  dl2x =  0.0;
  dl2y = -sr09*z;
  dl2z = -sr09*y;

  //del4f.m[3][4] = -sr01*(6.0*n2 - 8.0*l2 - 1.0)*mn;
  dl4  = -sr01*(6.0*z2 - 8.0*x2 - 1.0)*yz;
  dl4x =  16.0*sr01*x*yz;
  dl4y = -sr01*(6.0*z2 - 8.0*x2 - 1.0)*z;
  dl4z = -sr01*(18.0*z2 - 8.0*x2 - 1.0)*y;

  MAKE_ELEMENT_7(3,4);


  //del2f.m[3][5] = sr09*ln;
  dl2  = sr09*xz;
  dl2x = sr09*z;
  dl2y = 0.0;
  dl2z = sr09*x;

  //del4f.m[3][5] = sr01*(2.0*n2 - 8.0*l2 + 3.0)*ln;
  dl4  = sr01*(2.0*z2 - 8.0*x2 + 3.0)*xz;
  dl4x = sr01*(2.0*z2 - 24.0*x2 + 3.0)*z;
  dl4y = 0.0;
  dl4z = sr01*(6.0*z2 - 8.0*x2 + 3.0)*x;

  MAKE_ELEMENT_7(3,5);


  //del2f.m[3][6] = sr06*lm;
  dl2  = sr06*xy;
  dl2x = sr06*y;
  dl2y = sr06*x;
  dl2z = 0.0;

  //del4f.m[3][6] = sr96*(m2 - n2 - 0.25)*lm;
  dl4  =  sr96*(y2 - z2 - 0.25)*xy;
  dl4x =  sr96*(y2 - z2 - 0.25)*y;
  dl4y =  sr96*(3.0*y2 - z2 - 0.25)*x;
  dl4z = -2.0*sr96*z*xy;

  MAKE_ELEMENT_7(3,6);


  //del2f.m[3][7] = -sr015*(l2 - m2);
  dl2  = -sr015*(x2 - y2);
  dl2x = -2.0*sr015*x;
  dl2y =  2.0*sr015*y;
  dl2z =  0.0;

  //del4f.m[3][7] = sr375*(l2 - m2 + 1.4*m4 - 1.2*l2*m2 - l4);
  dl4  = sr375*(x2 - y2 + 1.4*y4 - 1.2*x2*y2 - x4);
  dl4x = sr375*(2.0 - 2.4*y2 - 4.0*x2)*x;
  dl4y = sr375*(-2.0 + 5.6*y2 - 2.4*x2)*y;
  dl4z = 0.0;

  MAKE_ELEMENT_7(3,7);


  //del2f.m[4][4] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[4][4] = (2.0 - 3.0*n2)*n2 - 4.0*l2*m2;
  dl4  = (2.0 - 3.0*z2)*z2 - 4.0*x2*y2;
  dl4x = -8.0*x*y2;
  dl4y = -8.0*x2*y;
  dl4z =  (4.0 - 12.0*z2)*z;

  MAKE_ELEMENT_7(4,4);

  //del2f.m[4][5] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[4][5] = 2.0*(l2 - m2)*lm;
  dl4  = 2.0*(x2 - y2)*xy;
  dl4x = 2.0*(3.0*x2 - y2)*y;
  dl4y = 2.0*(x2 - 3.0*y2)*x;
  dl4z = 0.0;

  MAKE_ELEMENT_7(4,5);


  //del2f.m[4][6] = sr15*ln;
  dl2  = sr15*xz;
  dl2x = sr15*z;
  dl2y = 0.0;
  dl2z = sr15*x;

  //del4f.m[4][6] = -sr15*(2.0*n2 - 1.0)*ln;
  dl4  = -sr15*(2.0*z2 - 1.0)*xz;
  dl4x = -sr15*(2.0*z2 - 1.0)*z;
  dl4y =  0.0;
  dl4z = -sr15*(6.0*z2 - 1.0)*x;

  MAKE_ELEMENT_7(4,6);


  //del2f.m[4][7] = sr15*mn;
  dl2  = sr15*yz;
  dl2x = 0.0;
  dl2y = sr15*z;
  dl2z = sr15*y;

  //del4f.m[4][7] = -sr15*(2.0*n2 - 1.0)*mn;
  dl4  = -sr15*(2.0*z2 - 1.0)*yz;
  dl4x =  0.0;
  dl4y = -sr15*(2.0*z2 - 1.0)*z;
  dl4z = -sr15*(6.0*z2 - 1.0)*y;

  MAKE_ELEMENT_7(4,7);


  //del2f.m[5][5] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[5][5] = -pow((2.0*n2 - 1.0),2) + 4.0*l2*m2;
  t1 = 2.0*z2 - 1.0;
  dl4  = -t1*t1 + 4.0*x2*y2;
  dl4x =  8.0*x*y2;
  dl4y =  8.0*x2*y;
  dl4z = -8.0*t1*z;

  MAKE_ELEMENT_7(5,5);


  //del2f.m[5][6] = -sr15*mn;
  dl2  = -sr15*yz;
  dl2x =  0.0;
  dl2y = -sr15*z;
  dl2z = -sr15*y;

  //del4f.m[5][6] = sr15*(2.0*n2 - 1.0)*mn;
  dl4  = sr15*(2.0*z2 - 1.0)*yz;
  dl4x = 0.0;
  dl4y = sr15*(2.0*z2 - 1.0)*z;
  dl4z = sr15*(6.0*z2 - 1.0)*y;

  MAKE_ELEMENT_7(5,6);

  //del2f.m[5][7] = sr15*ln;
  dl2  = sr15*xz;
  dl2x = sr15*z;
  dl2y = 0.0;
  dl2z = sr15*x;

  //del4f.m[5][7] = -sr15*(2.0*n2 - 1.0)*ln;
  dl4  = -sr15*(2.0*z2 - 1.0)*xz;
  dl4x = -sr15*(2.0*z2 - 1.0)*z;
  dl4y =  0.0;
  dl4z = -sr15*(6.0*z2 - 1.0)*x;

  MAKE_ELEMENT_7(5,7);


  //del2f.m[6][6] = -(3.0*n2 - 1.0)/2.0;
  dl2  =  0.5 - 1.5*z2;
  dl2x =  0.0;
  dl2y =  0.0;
  dl2z = -3.0*z;

  //del4f.m[6][6] = 1.5*(n2 - 1.0)*n2;
  dl4  = (1.5*z2 - 1.5)*z2;
  dl4x = 0.0;
  dl4y = 0.0;
  dl4z = (6.0*z2 - 3.0)*z;

  MAKE_ELEMENT_7(6,6);


  //del2f.m[6][7] = 0.0;
  dl2  = 0.0;
  dl2x = 0.0;
  dl2y = 0.0;
  dl2z = 0.0;

  //del4f.m[6][7] = 0.0;
  dl4  = 0.0;
  dl4x = 0.0;
  dl4y = 0.0;
  dl4z = 0.0;

  MAKE_ELEMENT_7(6,7);


  //del2f.m[7][7] = -(3.0*n2 - 1.0)/2.0;
  dl2  =  0.5 - 1.5*z2;
  dl2x =  0.0;
  dl2y =  0.0;
  dl2z = -3.0*z;

  //del4f.m[7][7] = 1.5*(n2 - 1.0)*n2;
  dl4  = (1.5*z2 - 1.5)*z2;
  dl4x = 0.0;
  dl4y = 0.0;
  dl4z = (6.0*z2 - 3.0)*z;

  MAKE_ELEMENT_7(7,7);
}

/************************************************************************/

/* ----------------------------------------------------------------------
 * Fast Model Generalized Pseudopotential Theory (MGPT) interatomic
 * potential routine.
 *
 * Copyright (2015) Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Tomas Oppelstrup (oppelstrup2@llnl.gov) and John Moriarty
 * (moriarty2@llnl.gov)
 * LLNL-CODE-674031  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (as published by the
 * Free Software Foundation) version 2, dated June 1991.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
 * GNU General Public License for more details.
 *
 * LLNL Preamble Notice
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was performed under the auspices
 * of the DOE by Lawrence Livermore National Laboratory under Contract No.
 * DE-AC52-07NA27344.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process,
 * or services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
------------------------------------------------------------------------- */
