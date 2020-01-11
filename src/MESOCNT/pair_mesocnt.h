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

  double ang,ang_inv,eunit,funit;
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

// inlined matrix functions

namespace MathExtra {

  // set matrix to zero

  inline void zeromat3(double m[3][3])
  {
    m[0][0] = m[0][1] = m[0][2] = 0.0;
    m[1][0] = m[1][1] = m[1][2] = 0.0;
    m[2][0] = m[2][1] = m[2][2] = 0.0;
  }

  inline void zeromat3(double **m)
  {
    m[0][0] = m[0][1] = m[0][2] = 0.0;
    m[1][0] = m[1][1] = m[1][2] = 0.0;
    m[2][0] = m[2][1] = m[2][2] = 0.0;
  }

  // add two matrices

  inline void plus3(const double m[3][3], double **m2,
                               double **ans)
  {
    ans[0][0] = m[0][0]+m2[0][0];
    ans[0][1] = m[0][1]+m2[0][1];
    ans[0][2] = m[0][2]+m2[0][2];
    ans[1][0] = m[1][0]+m2[1][0];
    ans[1][1] = m[1][1]+m2[1][1];
    ans[1][2] = m[1][2]+m2[1][2];
    ans[2][0] = m[2][0]+m2[2][0];
    ans[2][1] = m[2][1]+m2[2][1];
    ans[2][2] = m[2][2]+m2[2][2];
  }

  // subtract two matrices

  inline void minus3(const double m[3][3], const double m2[3][3],
                               double ans[3][3])
  {
    ans[0][0] = m[0][0]-m2[0][0];
    ans[0][1] = m[0][1]-m2[0][1];
    ans[0][2] = m[0][2]-m2[0][2];
    ans[1][0] = m[1][0]-m2[1][0];
    ans[1][1] = m[1][1]-m2[1][1];
    ans[1][2] = m[1][2]-m2[1][2];
    ans[2][0] = m[2][0]-m2[2][0];
    ans[2][1] = m[2][1]-m2[2][1];
    ans[2][2] = m[2][2]-m2[2][2];
  }

  inline void minus3(double **m, const double m2[3][3],
                               double ans[3][3])
  {
    ans[0][0] = m[0][0]-m2[0][0];
    ans[0][1] = m[0][1]-m2[0][1];
    ans[0][2] = m[0][2]-m2[0][2];
    ans[1][0] = m[1][0]-m2[1][0];
    ans[1][1] = m[1][1]-m2[1][1];
    ans[1][2] = m[1][2]-m2[1][2];
    ans[2][0] = m[2][0]-m2[2][0];
    ans[2][1] = m[2][1]-m2[2][1];
    ans[2][2] = m[2][2]-m2[2][2];
  }

  // compute outer product of two vectors

  inline void outer3(const double *v1, const double *v2,
                                double ans[3][3])
  {
    ans[0][0] = v1[0]*v2[0]; ans[0][1] = v1[0]*v2[1]; ans[0][2] = v1[0]*v2[2];
    ans[1][0] = v1[1]*v2[0]; ans[1][1] = v1[1]*v2[1]; ans[1][2] = v1[1]*v2[2];
    ans[2][0] = v1[2]*v2[0]; ans[2][1] = v1[2]*v2[1]; ans[2][2] = v1[2]*v2[2];
  }
}

#endif
#endif
