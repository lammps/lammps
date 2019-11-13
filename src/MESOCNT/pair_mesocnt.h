#ifdef PAIR_CLASS

PairStyle(mesocnt, PairMesoCNT)

#else 

#ifndef LMP_PAIR_MESOCNT_H
#define LMP_PAIR_MESOCNT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMesoCNT : public Pair {
 public:
  PairMesoCNT(class LAMMPS *);
  ~PairMesoCNT();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  int uinf_points,gamma_points,phi_points,usemi_points;
  double ang,angrec,e,erec,funit;
  double r,r2,d,rc,rc2,rc0;
  double d_ang,rc_ang,rc2_ang;
  double sig,eps,n_sig,comega,ctheta;
  double hstart_gamma,hstart_uinf,
	 hstart_phi,psistart_phi,hstart_usemi,xistart_usemi;
  double delh_gamma,delh_uinf,delh_phi,delpsi_phi,delh_usemi,delxi_usemi;
  
  double **uinf_coeff,**gamma_coeff,****phi_coeff,****usemi_coeff;
  double *p1,*p2,*param;
  double **flocal,**fglobal,**basis;

  void allocate();
  void sort(int *, int);
  void read_file(char *, double *, double &, double &, int);
  void read_file(char *, double **, double &, double &, 
		  double &, double &, int);

  void spline_coeff(double *, double **, double, int);
  void spline_coeff(double **, double ****, double, double, int);

  double spline(double, double, double, double **, int);
  double dspline(double, double, double, double **, int);
  double spline(double, double, double, double, double, double, 
		  double ****, int);
  double dxspline(double, double, double, double, double, double, 
		  double ****, int);
  double dyspline(double, double, double, double, double, double, 
		  double ****, int);

  void geominf(const double *, const double *, const double *, 
		  const double *, double *, double **);
  void geomsemi(const double *, const double *, const double *,
		  const double *, const double *, double *, double **);
  double weight(const double *, const double *, const double *,
		  const double *);

  double uinf(const double *);
  double usemi(const double *);
  void finf(const double *, double **);
  void fsemi(const double *, double **);

  // inlined functions for efficiency

  inline int heaviside(double x) {
    if (x > 0) return 1;
    else return 0;
  }

  inline double s(double x) {
    return heaviside(-x) + heaviside(x)*heaviside(1-x)*(1 - x*x*(3 - 2*x));
  }

  inline double ds(double x) {
    return 6 * heaviside(x) * heaviside(1-x) * x * (x-1);
  }

  inline double s5(double x) {
    double x2 = x * x;
    return heaviside(-x) + heaviside(1-x)*(1 - x2*x*(6*x2 - 15*x + 10));
  }

  inline double ds5(double x) {
    double x2 = x * x;
    return -30 * heaviside(x) * heaviside(1-x) * x2 * (x2 - 2*x + 1);
  }
};

}

#endif
#endif
