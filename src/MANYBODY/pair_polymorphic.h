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

PairStyle(polymorphic,PairPolymorphic)

#else

#ifndef LMP_PAIR_POLYMORPHIC_H
#define LMP_PAIR_POLYMORPHIC_H

#include "pair.h"

namespace LAMMPS_NS {

//===========================================
class C1function {
 public:
  C1function(){};
  virtual ~C1function() {};
  virtual double value(double x)=0;
  virtual double derivative(double x)=0;
};
//===========================================
class C1constant : public C1function {
 public:
  C1constant(double C): C_(C) {};
  C1constant(): C_(0) {};
  virtual double value(double x) { return C_; }
  virtual double derivative(double x) { return 0.; }
 protected:
  double C_;
};
//===========================================
class C1exponential : public C1function {
 public:
  C1exponential(double C, double lambda): C_(C),lambda_(lambda) {};
  C1exponential(): C_(0),lambda_(0) {};
  virtual double value(double x) { return C_*exp(lambda_*x); }
  virtual double derivative(double x) { return lambda_*C_*exp(lambda_*x); }
 protected:
  double C_,lambda_;
};
//===========================================
class C1sine : public C1function {
 public:
  C1sine(double C, double lambda): C_(C),lambda_(lambda) {};
  C1sine(): C_(0),lambda_(0) {};
  virtual double value(double x) { return C_*sin(lambda_*x); }
  virtual double derivative(double x) { return lambda_*C_*cos(lambda_*x); }
 protected:
  double C_,lambda_;
};
//===========================================
class C1cosine : public C1function {
 public:
  C1cosine(double C, double lambda): C_(C),lambda_(lambda) {};
  C1cosine(): C_(0),lambda_(0) {};
  virtual double value(double x) { return C_*cos(lambda_*x); }
  virtual double derivative(double x) { return -lambda_*C_*sin(lambda_*x); }
 protected:
  double C_,lambda_;
};
//===========================================
class C1tabularFunction : public C1function {
 public:
  C1tabularFunction():size(0),xmin(0),xmax(0) {
    resize(0);
  }
  C1tabularFunction(int n):size(0),xmin(0),xmax(0) {
    resize(n);
  }
  C1tabularFunction(int n, double x1, double x2):size(0),xmin(x1),xmax(x2) {
    resize(n);
  }
  virtual ~C1tabularFunction() {
    resize(0);
  }
  virtual double value(double x) {
    double y,dy;
    this->value(x,y,1,dy,0);
    return y;
  }
  virtual double derivative(double x) {
    double y,dy;
    this->value(x,y,0,dy,1);
    return dy;
  }
  void value(double x, double &y, int ny, double &y1, int ny1)
  {
    double ps = (x - xmin) * rdx;
    int ks = ps;
    if (ks > size-2) ks = size-2;
    ps = ps - ks;
    if (ps > 1.0) ps = 1.0;
    if (ny) y = ((ys3[ks]*ps + ys2[ks])*ps + ys1[ks])*ps + ys[ks];
    if (ny1) y1 = (ys6[ks]*ps + ys5[ks])*ps + ys4[ks];
  }
  void set_xrange(double x1, double x2) {
    xmin = x1;
    xmax = x2;
  }
  void set_values(int n, double x1, double x2, double * values)
  {
    resize(n);
    xmin = x1;
    xmax = x2;
    memcpy(ys,values,n*sizeof(double));
    initialize();
  }

  void print_value()
  {
     printf("%d %f %f %f \n",size,xmin,xmax,rdx);
     printf(" \n");
     for (int i = 0; i < size; i++) {
        printf("%f %f \n",xs[i],ys[i]);
     }
  }

  double get_xmin() {
    return xmin;
  }
  double get_xmax() {
    return xmax;
  }
  double get_rdx() {
    return rdx;
  }

 protected:
  void resize(int n) {
    if (size == 0) {
      xs = NULL;
      ys = NULL;
      ys1 = NULL;
      ys2 = NULL;
      ys3 = NULL;
      ys4 = NULL;
      ys5 = NULL;
      ys6 = NULL;
    }
    if (n == 0) {
      if (xs) delete [] xs;
      if (ys) delete [] ys;
      if (ys1) delete [] ys1;
      if (ys2) delete [] ys2;
      if (ys3) delete [] ys3;
      if (ys4) delete [] ys4;
      if (ys5) delete [] ys5;
      if (ys6) delete [] ys6;
      size = 0;
    }
    else if (n != size) {
      size = n;
      if (xs) delete [] xs;
      xs = new double[n];
      if (ys) delete [] ys;
      ys = new double[n];
      if (ys1) delete [] ys1;
      ys1 = new double[n];
      if (ys2) delete [] ys2;
      ys2 = new double[n];
      if (ys3) delete [] ys3;
      ys3 = new double[n];
      if (ys4) delete [] ys4;
      ys4 = new double[n];
      if (ys5) delete [] ys5;
      ys5 = new double[n];
      if (ys6) delete [] ys6;
      ys6 = new double[n];
    }
  }
  void initialize() {
    int n = size;
    rdx = (xmax-xmin)/(n-1.0);
    for (int i = 0; i < n; i++) {
      xs[i] = xmin+i*rdx;
    }
    rdx = 1.0 / rdx;
    ys1[0] = ys[1] - ys[0];
    ys1[1] = 0.5 * (ys[2] - ys[0]);
    ys1[n-2] = 0.5 * (ys[n-1] - ys[n-3]);
    ys1[n-1] = ys[n-1] - ys[n-2];
    for (int i = 2; i < n-2; i++) {
      ys1[i]=((ys[i-2]-ys[i+2])+ 8.0*(ys[i+1]-ys[i-1]))/12.0;
    }
    for (int i = 0; i < n-1; i++) {
      ys2[i]=3.0*(ys[i+1]-ys[i])-2.0*ys1[i]-ys1[i+1];
      ys3[i]=ys1[i]+ys1[i+1]-2.0*(ys[i+1]-ys[i]);
    }
    ys2[n-1]=0.0;
    ys3[n-1]=0.0;
    for (int i = 0; i < n; i++) {
      ys4[i]=ys1[i]*rdx;
      ys5[i]=2.0*ys2[i]*rdx;
      ys6[i]=3.0*ys3[i]*rdx;
    }
  }
  int size;
  double xmin,xmax,rdx;
  double * ys, * ys1, * ys2, * ys3, * ys4, * ys5, * ys6;
  double * xs;
};

//===========================================
class PairPolymorphic : public Pair {
 public:
  PairPolymorphic(class LAMMPS *);
  virtual ~PairPolymorphic();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  struct PairParameters {
    double cut;
    double cutsq;
    bool xi; // "indicator"
    class C1function * U;
    class C1function * V;
    class C1function * W;
    class C1function * P;
    class C1function * F;
    PairParameters() {
      cut = 0.0;
      cutsq = 0.0;
      xi =  true;
      U = NULL;
      V = NULL;
      W = NULL;
      P = NULL;
      F = NULL;
    };
  };
  struct TripletParameters {
    class C1function * G;
    TripletParameters() {
      G = NULL;
    };
  };

  bool eta; // global indicator
  int nx,nr,ng; // table sizes
  double maxX;

  // parameter sets
  PairParameters    * pairParameters;    // for I-J interaction
  TripletParameters * tripletParameters; // for I-J-K interaction

  char **elements;              // names of unique elements
  int **elem2param;             // map: element pairs to parameters
  int ***elem3param;            // map: element triplets to parameters
  int *type_map;                // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int npair,ntriple;

  void allocate();
  void read_line(FILE *, char *);
  void read_array(FILE *, int, double *);
  virtual void read_file(char *, char**);
  class C1function * create_function(char *, FILE *);

  void setup();

  void attractive(PairParameters *, TripletParameters *, double, double,
                  double, double *, double *, double *, double *, double *);

  void ters_zetaterm_d(double, double *, double, double *, double, double *,
                       double *, double *, PairParameters *, TripletParameters *);
  void costheta_d(double *, double, double *, double,
		  double *, double *, double *);

  // inlined functions for efficiency

  inline double vec3_dot(const double x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  inline void vec3_add(const double x[3], const double y[3],
		       double * const z) const {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }

  inline void vec3_scale(const double k, const double x[3],
			 double y[3]) const {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }

  inline void vec3_scaleadd(const double k, const double x[3],
			    const double y[3], double * const z) const {
    z[0] = k*x[0]+y[0];
    z[1] = k*x[1]+y[1];
    z[2] = k*x[2]+y[2];
  }
};

}

#endif
#endif
