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

#ifdef PAIR_CLASS

PairStyle(comb,PairComb)

#else

#ifndef LMP_PAIR_COMB_H
#define LMP_PAIR_COMB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairComb : public Pair {
 public:
  PairComb(class LAMMPS *);
  virtual ~PairComb();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

  double yasu_char(double *, int &);

 private:
  struct Param {
    double lam11,lam12,lam21,lam22;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga1,biga2,bigb1,bigb2;
    double bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    double plp1,plp3,plp6,a123,aconf;
    double rlm1,rlm2;
    double romiga,romigb,romigc,romigd,addrep;
    double QU1,QL1,DU1,DL1,Qo1,dQ1,aB1,bB1,nD1,bD1;
    double QU2,QL2,DU2,DL2,Qo2,dQ2,aB2,bB2,nD2,bD2;
    double chi,dj,dk,dl,dm,esm1,esm2,cmn1,cmn2,cml1,cml2;
    double coulcut, lcut, lcutsq, hfocor;
    int ielement,jelement,kelement;
    int powermint;
  };
  
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  double precision;
  Param *params;                // parameter set for an I-J-K interaction

  int nmax;
  double *qf;

  double *esm, **fafb, **dfafb, **ddfafb, **phin, **dphin, **erpaw;
  double *charge;
  int **intype, *typeno;
  int *NCo, cor_flag, cuo_flag, cuo_flag1, cuo_flag2;
  double **bbij;

  void allocate();
  virtual void read_file(char *);
  void setup();
  virtual void repulsive(Param *, double, double &, int, 
			 double &, double, double);
  double zeta(Param *, double, double, double *, double *);
  void force_zeta(Param *, int, int, int, double, double, double, double, 
		  double &, double &, double &);
  void attractive(Param *, double, double, double, double *, double *,
		  double *, double *, double *);
  double elp(Param *, double, double, double *, double *);
  void flp(Param *, double, double, double *, double *, double *, 
	   double *, double *);
  double comb_fc(double, Param *);
  double comb_fc_d(double, Param *);
  double comb_fc2(double);
  double comb_fc2_d(double);
  double comb_fc3(double);
  double comb_fc3_d(double);
  virtual double comb_fa(double, Param *, double,double);
  virtual double comb_fa_d(double, Param *, double,double);
  double comb_bij(double, Param *);
  double comb_bij_d(double, Param *);
  double comb_gijk(double, Param *);
  double comb_gijk_d(double, Param *);
  void comb_zetaterm_d(double, double *, double, double *, double,
			       double *, double *, double *, Param *);
  void costheta_d(double *, double, double *, double,
		  double *, double *, double *);
  double self(Param *, double, double);
  void sm_table();
  void potal_calc(double &, double &, double &);
  void tri_point(double, int &, int &, int &, double &, double &, 
		 double &, int &);
  void direct(int,int,int,int,double,double,double,double,double,double,
	double,double,double,double &,double &);
  void field(Param *,double,double,double,double &,double &);
  double qfo_self(Param *, double, double);
  void qfo_short(Param *, int, int, double, double, double,
		  double &, double &);
  void qfo_direct (int, int, int, int, double, double, double, double, 
	double, double &);
  void qfo_field(Param *, double,double ,double ,double &, double &);
  void qsolve(double *);
  void Over_cor(Param *, double, int, double &, double &);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_comm(int , int *, double *, int, int *);
  void unpack_comm(int , int , double *);

  // Short range neighbor list
  void add_pages(int );
  void Short_neigh();
  int npage, maxpage, pgsize, oneatom, **pages;
  int *sht_num, **sht_first;	// short-range neighbor list
  double cutmin;

  // vector functions, inline for efficiency

  inline double vec3_dot(double *x, double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }
  inline void vec3_add(double *x, double *y, double *z) {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }
  inline void vec3_scale(double k, double *x, double *y) {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }
  inline void vec3_scaleadd(double k, double *x, double *y, double *z) {
    z[0] = k*x[0]+y[0];  z[1] = k*x[1]+y[1];  z[2] = k*x[2]+y[2];
  }
};

}

#endif
#endif
