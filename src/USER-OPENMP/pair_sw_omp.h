/* -*- c++ -*- -----------------------------------------------------------
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

PairStyle(sw/omp,PairSWOMP)

#else

#ifndef LMP_PAIR_SW_OMP_H
#define LMP_PAIR_SW_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairSWOMP : public PairOMP {
 public:
  PairSWOMP(class LAMMPS *);
  ~PairSWOMP();
  virtual void compute(int, int);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);

  virtual double init_one(int, int);
  virtual void init_style();

 protected:
  struct Param {
    double epsilon,sigma;
    double littlea,lambda,gamma,costheta;
    double biga,bigb;
    double powerp,powerq;
    double tol;
    double cut,cutsq;
    double sigma_gamma,lambda_epsilon,lambda_epsilon2;
    double c1,c2,c3,c4,c5,c6;
    int ielement,jelement,kelement;
  };
  
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

 protected:
  void allocate();
  void read_file(char *);
  void setup();
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
/* ---------------------------------------------------------------------- */

  template <int EFLAG>
  void twobody(const Param *param, const double rsq, double &fforce, double &eng)
  {
    double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

    r = sqrt(rsq);
    rinvsq = 1.0/rsq;
    rp = pow(r,-param->powerp);
    rq = pow(r,-param->powerq);
    rainv = 1.0 / (r - param->cut);
    rainvsq = rainv*rainv*r; 
    expsrainv = exp(param->sigma * rainv);
    fforce = (param->c1*rp - param->c2*rq +
	      (param->c3*rp -param->c4*rq) * rainvsq) * expsrainv * rinvsq;
    if (EFLAG) eng = (param->c5*rp - param->c6*rq) * expsrainv;
};

/* ---------------------------------------------------------------------- */
  template <int EFLAG>
  void threebody(const Param *paramij, const Param *paramik, const Param *paramijk, 
		 const double rsq1, const double rsq2,
		 double *delr1, double *delr2, 
		 double *fj, double *fk, double &eng)
    {
      double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
      double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
      double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
      double facang,facang12,csfacang,csfac1,csfac2;

      r1 = sqrt(rsq1);
      rinvsq1 = 1.0/rsq1;
      rainv1 = 1.0/(r1 - paramij->cut);
      gsrainv1 = paramij->sigma_gamma * rainv1;
      gsrainvsq1 = gsrainv1*rainv1/r1; 
      expgsrainv1 = exp(gsrainv1);

      r2 = sqrt(rsq2);
      rinvsq2 = 1.0/rsq2;
      rainv2 = 1.0/(r2 - paramik->cut);
      gsrainv2 = paramik->sigma_gamma * rainv2;
      gsrainvsq2 = gsrainv2*rainv2/r2; 
      expgsrainv2 = exp(gsrainv2);

      rinv12 = 1.0/(r1*r2);
      cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
      delcs = cs - paramijk->costheta;
      delcssq = delcs*delcs;

      facexp = expgsrainv1*expgsrainv2;

      // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
      //          facexp*delcssq;

      facrad = paramijk->lambda_epsilon * facexp*delcssq;
      frad1 = facrad*gsrainvsq1;
      frad2 = facrad*gsrainvsq2;
      facang = paramijk->lambda_epsilon2 * facexp*delcs;
      facang12 = rinv12*facang;
      csfacang = cs*facang;
      csfac1 = rinvsq1*csfacang;
  
      fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
      fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
      fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;
  
      csfac2 = rinvsq2*csfacang;

      fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
      fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
      fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;

      if (EFLAG) eng = facrad;
    };
};

}

#endif
#endif
