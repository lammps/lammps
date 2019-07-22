#ifdef PAIR_CLASS

PairStyle(mesocnt, PairMesoCNT)

#else

#ifndef LMP_PAIR_MESOCNT_H
#define LMP_PAIR_MESOCNT_H

#include "pair.h"
#include <vector>

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
  int me, n, gamma_points, pot_points;
  int uinf_points, usemi_points, phi_points;
  int redlist_size,chain_size,end_size;
  int *redlist,*nchain,*end;
  int **chain;
  double cutoff, cutoffsq;
  double angstrom, angstromrec, qelectron, qelectronrec, forceunit;
  double diameter_angstrom, cutoff_angstrom, cutoffsq_angstrom;
  double sigma, epsilon, n_sigma, radius, radiussq, diameter, rc, rc0, comega, ctheta;
  double start_gamma, start_uinf, startxi_usemi, starth_phi;
  double del_gamma, del_uinf, delxi_usemi, delh_phi;
  double *starth_usemi, *startzeta_phi;
  double *delh_usemi, *delzeta_phi;
  double *gamma_data, *uinf_data;
  double *m, *p1, *p2, *param, *wvector;
  double *dw1, *dw2, *sumdw1, *sumdw2;
  double *fchain11, *fchain12, *fchain21, *fchain22, *fchain1, *fchain2;
  double **dw1q1, **dw1q2, **dw2q1, **dw2q2;
  double **sumdw1q1, **sumdw1q2, **sumdw2q1, **sumdw2q2;
  double **sumdw1p1, **sumdw1p2, **sumdw2p1, **sumdw2p2;
  double **dp1_dr1, **dp2_dr1, **dp1_dr2, **dp2_dr2;
  double **usemi_data, **phi_data;
  double **gamma_coeff, **uinf_coeff;
  double **flocal,**fglobal,**basis;
  double ***usemi_coeff, ****phi_coeff;
  char *gamma_file, *uinf_file, *usemi_file, *phi_file;

  void allocate();
  
  double spline(double, double, double, double **, int);
  double spline(double, double, double *, double, double *, 
		  double, double ***, int);
  double spline(double, double, double, double, double, 
		  double, double ****, int);
  double dspline(double, double, double, double **, int);
  double dxspline(double, double, double *, double, double *, 
		  double, double ***, int);
  double dyspline(double, double, double *, double, double *, 
		  double, double ***, int);
  double dxspline(double, double, double, double, double, 
		  double, double ****, int);
  double dyspline(double, double, double, double, double, 
		  double, double ****, int);

  void spline_coeff(double *, double **, int);
  void spline_coeff(double **, double ***, int);
  void spline_coeff(double **, double ****, double, double, int);
  
  void read_file(char *, double *, double *, double *, int);
  void read_file(char *, double **, double *, double *, 
		  double *, double *, int);

  double uinf(double *);
  double usemi(double *);
  void finf(double *, double *, double *, double *, double *, double *, double **);
  void fsemi(double *, double *, double *, double *, double *, double *, double **);

  void geominf(const double *, const double *, const double *, 
		  const double *, double *, double *, double **);
  void geomsemi(const double *, const double *, const double *, const double *, 
		  const double *, double *, double *, double **);
  void weight(const double *, const double *,
		  const double *, const double *, double *);
  int heaviside(double);
  double s(double);
  double ds(double);
  double s5(double);
  double ds5(double);
  void sort(int *, int);

  void zero(double **);
  void plus(const double * const *, const double * const *, double **);
  void scalar(double, double**);
  void outer(const double *, const double *, double **);
  void matrix_vector(const double * const *, const double *, double *);
  void trans_matrix_vector(const double * const *, const double *, double *);
};

}
#endif
#endif
