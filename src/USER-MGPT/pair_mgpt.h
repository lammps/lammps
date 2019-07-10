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
   Contributing authors: Tomas Oppelstrup, LLNL (oppelstrup2@llnl.gov)
   and John Moriarty, LLNL (moriarty2@llnl.gov)

   Fast MGPT algorithm developed by Tomas Oppelstrup (2015) based on the
   matrix MGPT v4.4 FORTRAN routine of John Moriarty (2006) as converted
   to C++ for LAMMPS application by Jamie Marian and Alexander Stukowski
   (2011).  See LLNL copyright notice at bottom of this file.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mgpt,PairMGPT)

#else

#ifndef LMP_PAIR_MGPT_H
#define LMP_PAIR_MGPT_H

#include <new>
#include <cmath>
#include <cstdlib>
#include <cassert>

#include "pair.h"
#include "domain.h"

#include "mgpt_readpot.h"
#include "mgpt_linalg.h"


namespace LAMMPS_NS {


class PairMGPT : public Pair {
  mgpt_linalg linalg;
public:

  class Doublet {
  public:
    int i,j;
  public:
    Doublet(const Doublet &t) : i(t.i),j(t.j) {}
    Doublet(int ii,int jj) : i(ii < jj ? ii:jj),j(ii < jj ? jj : ii) {}

    Doublet operator=(const Doublet &t) {
      i = t.i;
      j = t.j;
      return *this;
    }
    int operator==(const Doublet &b) const {
      return (i == b.i) && (j == b.j);
    }
    int hash() const { return i*333331 + j*331; }
  };

  template<typename T,typename K> class Hash {

    class Link {
    public:
      T data;
      Link *next;
      K key;
      int hits;
      Link(const K &k,Link *n) : next(n),key(k),hits(1) {}

      static void *operator new(std::size_t sz) {
        const size_t align = 32;
        size_t x = (size_t) (void *) ::operator new(sz+align);
        size_t y = (x + align) - ((x+align)&(align-1));
        assert(sizeof(void *) <= align);
        assert((x & (sizeof(void *)-1)) == 0);
        ((void **) y)[-1] = (void *) x;
        return (void *) y;
      }
      static void operator delete(void *ptr) {
        ::operator delete(((void **) ptr)[-1]);
      }
    };

    int isprime(int x) {
      if(x%2 == 0)
        return 0;
      else {
        int k = 3;
        while(k*k <= x) {
          if(x%k == 0) return 0;
          k = k+2;
        }
        return 1;
    }
  }

    int size,used;
    Link **table;

    int maxlength,nstep,nsearch;
  public:

    class Iterator {
      Hash &H;
      int idx;
      Link *p;
    public:
      Iterator(Hash &HH) : H(HH),idx(-1),p(0) { next(); }
      Iterator(Hash &HH,int iidx,Link *pp) : H(HH),idx(iidx),p(pp) {}
      void next() {
        if(idx >= H.Size()) return;
        if(p != 0) p = p->next;
        if(p == 0) {
          do {
            idx = idx+1;
            if(idx >= H.Size()) return;
            p = H.table[idx];
          } while(p == 0);
        }
      }
      K *key() { return &p->key; }
      T *data() { return &p->data; }
      Link *link() { return p; }

      int operator==(const Iterator &a) {
        return idx==a.idx && p==a.p;
      }
      int operator!=(const Iterator &a) {
        return !(*this == a);
      }
    };

    Hash(int sz) {
      while(!isprime(sz)) sz = sz + 1;
      size = sz;
      used = 0;


      table = new Link *[size];
      for(int i = 0; i<size; i++)
        table[i] = 0;

      /* Counters for statistics */
      maxlength = 0;
      nstep = 0;
      nsearch = 0;
    }

    ~Hash() {
      for(int i = 0; i<size; i++) {
        Link *p = table[i];
        while(p != 0) {
          Link *q = p->next;
        delete p;
        p = q;
        }
      }
      delete[] table;
    }

    Iterator begin() { return Iterator(*this); }
    Iterator end() { return Iterator(*this,size,0); }

    int Size() { return size; }
    int Used() { return used; }
    int NSearch() { return nsearch; }
    int MaxLength() { return maxlength; }
    int NStep() { return nstep; }
    T * Insert(const K &key) {
      int idx = key.hash() % size;
      if(idx < 0) idx = idx + size;
      if(idx >= size || idx < 0) {
        printf("(1) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
        exit(1);
      }

      used = used + 1;
      if(1) {
        table[idx] = new Link(key,table[idx]);
        return &table[idx]->data;
      } else { /* This is for threading... and incomplete */
        typedef Link *LinkPtr;
        LinkPtr ptr = table[idx],last = 0,dataptr = new Link(key,0);

        while(ptr != 0) {
          last = ptr;
          ptr = ptr->next;
        }
        *((volatile LinkPtr *) &(last->next)) = dataptr;
        return &(dataptr->data);
      }
    }
    void Remove(const K &key) {
      int idx = key.hash() % size;
      Link *p,*last = 0;
      int count = 1;
      if(idx < 0) idx = idx + size;
      if(idx >= size || idx < 0) {
        printf("(2) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
        exit(1);
      }

      p = table[idx];
      while(p != 0 && !(p->key == key)) {
        last = p;
        p = p->next;
        count = count + 1;
      }

      if(p != 0) {
        used = used - 1;
        if(last == 0)
          table[idx] = p->next;
        else
          last->next = p->next;
        delete p;
      }

      if(count > maxlength)
        maxlength = count;
      nsearch = nsearch + 1;
      nstep = nstep + count;
    }
    T * Lookup(const K &key) {
      int idx = key.hash() % size;
      Link *p;
      int count = 1;
      if(idx < 0) idx = idx + size;
      if(idx >= size || idx < 0) {
        printf("(3) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
        exit(1);
      }


      p = table[idx];
      while(p != 0 && !(p->key == key)) {
        p = p->next;
        count = count + 1;
      }

      if(count > maxlength)
        maxlength = count;
      nsearch = nsearch + 1;
      nstep = nstep + count;

      if(p != 0) p->hits++;

      return (p == 0) ? 0 : &p->data;
    }
  };

 public:
  PairMGPT(class LAMMPS *);
  ~PairMGPT();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);

 private:

  void read_files(const char* parminFile, const char* potinFile, double vol);
  void allocate();

  struct Matrix {
    static int sz;

    double m[8][8];

    int align_check() {
      return ((((unsigned long long int) m) & 31) > 0);
    }
    void zero() {
      for(int i = 0; i<8; i++)
        for(int j = 0; j<8; j++)
          m[i][j] = 0.0;
    }

    void operator=(const Matrix &A) {
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          m[i][j] = A.m[i][j];
    }
    void operator=(double x) {
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          m[i][j] = x;
    }
    Matrix operator+(const Matrix &B) const {
      Matrix s;
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          s.m[i][j] = m[i][j] + B.m[i][j];
      return s;
    }
    Matrix operator-(const Matrix &B) const {
      Matrix s;
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          s.m[i][j] = m[i][j] - B.m[i][j];
      return s;
    }
    Matrix operator-() const {
      Matrix s;
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          s.m[i][j] =  -m[i][j];
      return s;
    }
    Matrix operator*(double x) const {
      Matrix P;
      for(int i = 1; i<=sz; i++)
        for(int j = 0; j<=sz; j++)
          P.m[i][j] = m[i][j] * x;
      return P;
    }
    Matrix operator/(double x) const {
      return (*this) * (1.0/x);
    }
    Matrix transpose() const {
      Matrix T;
      for(int i = 1; i<=sz; i++)
        for(int j = 1; j<=sz; j++)
          T.m[j][i] = m[i][j];
      return T;
    }
  };
  Matrix transpose(const Matrix &A) { return A.transpose(); }

  /* Preprocessor stuff to set alignment requirements on bonda_data
     and triplet_data structures. Without alignmnt, optimized algebra
     routines can fail.
  */
#ifdef __GNUC__
  #define PREALIGN struct
  #define POSTALIGN __attribute__((__aligned__(32)))
#elif defined(__bgq__)
  #define PREALIGN __align(32) struct
  #define POSTALIGN
#elif defined(__INTEL_COMPILER)
  /* This will probably not be used, since the Intel compiler also defines __GNUC__ */
  #define PREALIGN struct __declspec(align(32))
  #define POSTALIGN
#else
  /* Try GNU syntax anyway... */
  #define PREALIGN struct
  #define POSTALIGN __attribute__((__aligned__(32)))
  /*
  #define PREALIGN struct
  #define POSTALIGN
  */
#endif

  PREALIGN /*struct*/ bond_data {
    Matrix H,Hx,Hy,Hz;
    double fl_deriv_sum;
    double pad[3];
  } POSTALIGN;
  PREALIGN /*struct*/ triplet_data {
    Matrix  H1H2;
    Matrix  H1xH2 , H1yH2 , H1zH2;
    Matrix  H1H2x , H1H2y , H1H2z;

    int align_check() {
      return
        (H1H2.align_check()  << 0) |
        (H1xH2.align_check() << 1) |
        (H1yH2.align_check() << 2) |
        (H1zH2.align_check() << 3) |
        (H1H2x.align_check() << 4) |
        (H1H2y.align_check() << 5) |
        (H1H2z.align_check() << 6) ;
    }

    void zero() {
      H1H2.zero();
      H1xH2.zero(); H1yH2.zero(); H1zH2.zero();
      H1H2x.zero(); H1H2y.zero(); H1H2z.zero();
    }
  } POSTALIGN;

#undef PREALIGN
#undef POSTALIGN

  void make_bond(const double xx[][3],int i,int j,bond_data *bptr);
  void make_triplet(bond_data *ij_bond,bond_data *ik_bond,triplet_data *triptr);
  triplet_data *get_triplet(const double xx[][3],int i,int j,int k,
                            Hash<bond_data,Doublet> *bhash,triplet_data *twork,
                            double *dvir_ij_p,double *dvir_ik_p);

  int c1_outside(const double a[3],
                 int triclinic,const double alpha[3]) {
    const double stol = 1e-5;

    if(triclinic) {
      for(int p = 0; p<3; p++) {
        double cog = a[p];
        if(cog < domain->sublo_lamda[p]-0.5*rmax*alpha[p]-stol) return 1;
        if(cog > domain->subhi_lamda[p]+0.5*rmax*alpha[p]+stol) return 1;
      }

    } else {
      double rout = 0.0;


      for(int p = 0; p<3; p++) {
        double cog = a[p];
        if(cog < domain->sublo[p]-0.5*rmax-stol) return 1;
        if(cog > domain->subhi[p]+0.5*rmax+stol) return 1;

        if(cog < domain->sublo[p]-stol) {
          double t = cog - (domain->sublo[p]-stol);
          rout = rout + t*t;
        } else if(cog > domain->subhi[p]+stol) {
          double t = cog - (domain->subhi[p]+stol);
          rout = rout + t*t;
        }

      }

      if(rout > 0.25*rmax*rmax)
        return 1;
    }

    return 0;
  }
  int c2_outside(const double a[3],const double b[3],
                 int triclinic,const double alpha[3]) {
    const double stol = 1e-5;

    if(triclinic) {
      for(int p = 0; p<3; p++) {
        double cog = 0.5*(a[p] + b[p]);
        if(cog < domain->sublo_lamda[p]-0.5*rcrit*alpha[p]-stol) return 1;
        if(cog > domain->subhi_lamda[p]+0.5*rcrit*alpha[p]+stol) return 1;
      }
    } else {
      double rout = 0.0;

      for(int p = 0; p<3; p++) {
        double cog = 0.5*(a[p] + b[p]);
        if(cog < domain->sublo[p]-0.5*rcrit-stol) return 1;
        if(cog > domain->subhi[p]+0.5*rcrit+stol) return 1;

        if(cog < domain->sublo[p]-stol) {
          double t = cog - (domain->sublo[p]-stol);
          rout = rout + t*t;
        } else if(cog > domain->subhi[p]+stol) {
          double t = cog - (domain->subhi[p]+stol);
          rout = rout + t*t;
        }

      }

      if(rout > 0.25*rcrit*rcrit)
        return 1;
    }

    return 0;
  }
  double get_weight(const int triclinic,
                    const double a[3] = 0,const double b[3] = 0,
                    const double c[3] = 0,const double d[3] = 0) {
    const double
      *s0 = triclinic ? domain->sublo_lamda : domain->sublo,
      *s1 = triclinic ? domain->subhi_lamda : domain->subhi;
    double weight = 1.0;
    const double stol = 1e-5;

    for(int p = 0; p<3; p++) {
      double cog = 0.0,q,w,n = 0.0;
      if(a != 0) { cog = cog + a[p]; n = n + 1; }
      if(b != 0) { cog = cog + b[p]; n = n + 1; }
      if(c != 0) { cog = cog + c[p]; n = n + 1; }
      if(d != 0) { cog = cog + d[p]; n = n + 1; }
      cog = cog * (1.0/n);

      if(cog < 0.5*(s0[p]+s1[p])) q = cog - s0[p];
      else q = s1[p] - cog;

      w = q*(0.5/stol) + 0.5;
      if(w > 1.0) w = 1.0;
      if(w < 0.0) w = 0.0;
      weight = weight * w;
    }
    return weight;
  }

  void force_debug_3t(double xx[][3],
                      int i0,int j0,int k0,
                      int i ,int j ,int k ,
                      double dfix,double dfiy,double dfiz,
                      double dfjx,double dfjy,double dfjz,
                      double dfkx,double dfky,double dfkz);

  void force_debug_3v(double xx[][3],
                      int i0,int j0,int k0,
                      int i ,int j ,int k ,
                      double dfix,double dfiy,double dfiz,
                      double dfjx,double dfjy,double dfjz,
                      double dfkx,double dfky,double dfkz);

  void force_debug_4(double xx[][3],
                     int i0,int j0,int k0,int m0,
                     int i ,int j ,int k ,int m ,
                     double dfix,double dfiy,double dfiz,
                     double dfjx,double dfjy,double dfjz,
                     double dfkx,double dfky,double dfkz,
                     double dfmx,double dfmy,double dfmz);

  double numderiv3t(double xx[][3],int i,int j,int k,int p);
  double numderiv3v(double xx[][3],int i,int j,int k,int p,int ipert);
  double numderiv4(double xx[][3],int i,int j,int k,int m,int p);
  void compute_x(const int *nnei,const int * const *nlist,
                 double *e_s,double *e_p,double *e_t,double *e_q,
                 int evflag,int newton_pair);

  /* Reimplementation of bond matrix computation */
  void fl_deriv_new(double r,double ri,double xhat,double yhat,double zhat,
                    double &fl_0,double &fl_x,double &fl_y,double &fl_z,
                    double &fl_rp,double &fl_p1,double &fl_r0,double &fl_al);
  void hamltn_5_raw(const double xin,const double yin,const double zin,
                    double M [8][8],double Mx[8][8],
                    double My[8][8],double Mz[8][8],
                    double *fl_deriv_sum_p);
  void hamltn_7_raw(const double xin,const double yin,const double zin,
                    double M [8][8],double Mx[8][8],
                    double My[8][8],double Mz[8][8],
                    double *fl_deriv_sum_p);
  /* * */
  // Old matrix routines, only used in force debug routines.

  /// This function calculates the matrix product of ha and hb.
  inline Matrix prodmat(const Matrix& ha, const Matrix& hb) const {
    Matrix h;
    for(int l = 1; l <= lmax; l++) {
      for(int n = 1; n <= lmax; n++) {
        h.m[l][n] = 0.0;
        for(int m = 1; m <= lmax; m++)
          h.m[l][n] += ha.m[l][m] * hb.m[m][n];
      }
    }
    return h;
  }

  /// This function calculates the trace of the matrix product of ha and hb.
  inline double trace(const Matrix& ha, const Matrix& hb) const {
    double zquan = 0.0;
    for(int n = 1; n <= lmax; n++) {
      double cquan = 0.0;
      for(int m = 1; m <= lmax; m++)
        cquan += ha.m[n][m] * hb.m[m][n];
      zquan += cquan;
    }
    return zquan;
  }
  /* * */

  inline void transprod(const Matrix& a,const Matrix& b,Matrix &c) const
    {
      int i,j,k;

      if(lmax == 5) {
        const int n = 5;
        for(i = 1; i<=n; i++)
          for(j = 1; j<=n; j++) {
            double s = 0.0;
            for(k = 1; k<=n; k++)
              s = s + a.m[i][k]*b.m[j][k];
            c.m[i][j] = s;
              }
      } else if(lmax == 7) {
        const int n = 7;
        for(i = 1; i<=n; i++)
          for(j = 1; j<=n; j++) {
            double s = 0.0;
            for(k = 1; k<=n; k++)
              s = s + a.m[i][k]*b.m[j][k];
            c.m[i][j] = s;
          }
      } else {
        const int n = lmax;
        for(i = 1; i<=n; i++)
          for(j = 1; j<=n; j++) {
            double s = 0.0;
            for(k = 1; k<=n; k++)
              s = s + a.m[i][k]*b.m[j][k];
            c.m[i][j] = s;
          }
      }
    }

  inline double transtrace_s(const float (*A)[8],const float (*B)[8]) const {
    const int n = lmax;
    double s = 0.0;
    int i,j;

    for(i = 0; i<n; i++)
      for(j = 1; j<=n; j++)
        s = s + A[i][j]*B[i][j];
    return s;
  }
  inline double transtrace(const Matrix& a,const Matrix& b) const
    {
      int i,k;
      double s = 0.0;

      if(linalg.single)
        return transtrace_s((const float (*)[8]) &a.m[1][0],(const float (*)[8]) &b.m[1][0]);

      //printf("Calling transtrace... That is shit\n");

      if(lmax == 5) {
        const int n = 5;
        for(i = 1; i<=n; i++)
          for(k = 1; k<=n; k++)
            s = s + a.m[i][k]*b.m[i][k];
      } else if(lmax == 7) {
        const int n = 7;
        for(i = 1; i<=n; i++)
          for(k = 1; k<=n; k++)
            s = s + a.m[i][k]*b.m[i][k];
      } else {
        const int n = lmax;
        for(i = 1; i<=n; i++)
          for(k = 1; k<=n; k++)
            s = s + a.m[i][k]*b.m[i][k];
      }
      return s;
    }

  double cutoff;
  //double vpair[601];
  //double ktan[601];
  //double dvdvol[601];
  //double delr;
  //double rp,p1,al,r0,vc,vd,ve,pc,pd,pe;
  double rmax;
  double rcrit;
  //double evol0, pvol0, pot_input_vol;
  //double epsr, epsa;

  int thrion;
  int fourion;

  //int mode;
  double ddl[5];
  int lmax, lang;
  Matrix del0;
  double anorm3, anorm4;

  // Flag indicating whether volumetric pressure should be used.
  // Volumetric pressure means that terms emanating from the
  // derivative of the energy with respect to the potential atomic
  // volume parameter is included.
  int volpres_flag,nbody_flag;

  potdata splinepot;
};

}

#endif
#endif

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
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2, dated June 1991.
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
