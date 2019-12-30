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
  int nlocal_size,reduced_neigh_size;
  int n;
  int *reduced_nlist,*numchainlist;
  int **reduced_neighlist,**nchainlist,**endlist;
  int ***chainlist;

  double ang,ang_inv,e,e_inv,funit;
  double r,rsq,d,rc,rcsq,rc0,cutoff,cutoffsq;
  double r_ang,rsq_ang,d_ang,rc_ang,rcsq_ang,cutoff_ang,cutoffsq_ang;
  double sig,comega,ctheta;
  double hstart_uinf,hstart_gamma,
	 hstart_phi,psistart_phi,hstart_usemi,xistart_usemi;
  double delh_uinf,delh_gamma,delh_phi,delpsi_phi,delh_usemi,delxi_usemi;
  
  double p1[3],p2[3],p[3],m[3];
  double *param,*w,*wnode;
  double **dq_w;
  double ***q1_dq_w,***q2_dq_w;
  double **uinf_coeff,**gamma_coeff,****phi_coeff,****usemi_coeff;
  double **flocal,**fglobal,**basis;

  char *uinf_file,*gamma_file,*phi_file,*usemi_file;

  void allocate();
  void bond_neigh();
  void neigh_common(int, int, int &, int *);
  void chain_split(int *, int, int &, int **, int *, int *);
  void sort(int *, int);
  void read_file(const char *, double *, double &, double &, int);
  void read_file(const char *, double **, double &, double &, 
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

  void geometry(const double *, const double *, const double *,
		            const double *, const double *, 
                double *, double *, double *, double **);
  void weight(const double *, const double *, const double *,
		          const double *, double &, double *, double *, 
              double *, double *);

  void finf(const double *, double &, double **);
  void fsemi(const double *, double &, double &, double **);

  // inlined functions for efficiency

  inline double heaviside(double x) {
    if (x > 0) return 1.0;
    else return 0.0;
  }

  inline double s(double x) {
    return heaviside(-x) + heaviside(x)*heaviside(1-x)*(1 - x*x*(3 - 2*x));
  }

  inline double ds(double x) {
    return 6 * heaviside(x) * heaviside(1-x) * x * (x-1);
  }

  inline double s5(double x) {
    double x2 = x * x;
    return heaviside(-x) 
	    + heaviside(x)*heaviside(1-x)*(1 - x2*x*(6*x2 - 15*x + 10));
  }

  inline double ds5(double x) {
    double x2 = x * x;
    return -30 * heaviside(x) * heaviside(1-x) * x2 * (x2 - 2*x + 1);
  }
};

}

#endif
#endif
