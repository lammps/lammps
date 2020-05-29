/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(smtbq,PairSMTBQ)

#else

#ifndef LMP_PAIR_SMTBQ_H
#define LMP_PAIR_SMTBQ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSMTBQ : public Pair {
 public:
  PairSMTBQ(class LAMMPS *);
  virtual ~PairSMTBQ();
  virtual void compute(int, int);

  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();


protected:
  struct Param {
    double sto,n0,ne,chi,dj;
    double R,dzeta;   //Rayon slater
    double cutsq;
    double qform, masse;                   // Charge formelle
    char *nom;
  };

  struct Intparam {
    double a,p,ksi,q,r0,dc1,dc2;
    double abuck,rhobuck,neig_cut,aOO,bOO,r1OO,r2OO;
    char *typepot,*mode;
    int intsm;
  };

  double rmin,dr,ds;            // table parameter
  int kmax;
  bigint Qstep;                 // Frequency of charge resolution
  double precision;             // accuracy of convergence
  int loopmax;                  // max of iteration

  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  char *QEqMode;                // name of QEqMode
  char *Bavard;                 // Verbose parameter
  char *writepot;               // write or not the electronegativity component
  char *writeenerg;         // write or not the energy component
  char *QInitMode;              // mode of initialization of charges
  double zlim1QEq;              // z limit for QEq equilibration
  double zlim2QEq;              // z limit for QEq equilibration
  double QOxInit;               // Initial charge for oxygen atoms (if the charge is not specified)
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  int maxintparam;          // max # of interaction sets
  int maxintsm;             // max # of interaction SM
  double r1Coord,r2Coord;
  Param *params;                // parameter set for an I atom
  Intparam *intparams;          // parameter set for an I interaction

  int nmax,*nQEqall,*nQEqaall,*nQEqcall;
  double *qf,*q1,*q2,Nevery,Neverypot;

// Coulombian interaction
   double *esm, **fafb, **dfafb, *fafbOxOxSurf, *dfafbOxOxSurf, *fafbTiOxSurf,*dfafbTiOxSurf;
   double *potqn, *dpotqn, Vself, *Zsm,*coord, *fafbOxOxBB, *dfafbOxOxBB,*fafbTiOxBB, *dfafbTiOxBB ;
   int **intype, **coultype;
   int *NCo;
   double coordOxBulk,coordOxSurf,ROxSurf,coordOxBB,ROxBB;

// Covalent interaction
   double *ecov, *potmad, *potself, *potcov; //, *chimet;
   double **tabsmb,**dtabsmb, **tabsmr, **dtabsmr, *sbcov, *sbmet;
   double ncov;

// Neighbor Table
   int nteam,cluster,*hybrid;
   int *nvsm, **vsm, *flag_QEq;

// Parallelisation
   int me, nproc;
   double *tab_comm;

// GAMMAS function
   double *fct;

// HERE its routines
// =====================================
   void allocate();
   virtual void read_file(char *);

   void tabsm();
   void tabqeq();

   void potqeq(int, int, double, double, double,
               double &, int, double &);
   void pot_ES (int, int, double, double &);
   void pot_ES2 (int, int, double, double &);

   double self(Param *, double);
   double qfo_self(Param *, double);

   virtual void repulsive(Intparam *, double, int, int, double &, int, double &);
   virtual void rep_OO (Intparam *, double, double &, int, double &);
   virtual void Attr_OO (Intparam *, double, double &, int, double &);
   virtual void attractive(Intparam *, double, int, int, double, int, double );
   void f_att(Intparam *, int, int, double, double &) ;
   void pot_COV(Param *, int, double &);
   double potmet(Intparam *, double, int, double, int, double);

   double fcoupure(double, double, double);
   double fcoupured(double, double, double);

   double fcoup2(double, double, double);
   double Intfcoup2(double, double, double);
   double Primfcoup2(double, double, double);

   void groupBulkFromSlab_QEq();
   void groupQEqAllParallel_QEq();
   void groupQEqAll_QEq();
   void groupSurface_QEq();
   void QForce_charge(int);
   void Charge();
   void Init_charge(int*, int*, int*);
   void CheckEnergyVSForce();

// ===========================================
// Communication pack
   int pack_forward_comm (int, int*, double*, int, int*);
   void unpack_forward_comm (int, int, double*);
   int pack_reverse_comm (int, int, double*);
   void unpack_reverse_comm (int, int*, double*);
   void forward (double*); void reverse (double*);
   void forward_int (int*); void reverse_int (int*);

   int Tokenize( char* , char***  );

   inline double vec3_dot(const double x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
   }

   template <class T> const T& min ( const T& a, const T& b ) {
   return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
   }

// =============================================
// Gammas function
   void gammas(double &, double &, double &, double &, double &, double &,
                  double &, double &, double &, double &, double &, double &,
                  double &, double &, double &);

   void css(double &, double, double, double, double, double,
                    double &, double &, double &, double &,double &,
                    double &, double &, double &, double &);

   double coeffs(int, int, int);

   double binm(int, int);
   void caintgs (double, int, double *);
   void cbintgs( double, int, double *);

// =====================================
// short range neighbor list

  void add_pages(int howmany = 1);
//  void Short_neigh();
  int maxpage, pgsize, oneatom, **pages;
  double cutmin;
};

}

#endif
#endif
