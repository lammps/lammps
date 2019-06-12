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
  int n, gamma_points, pot_points;
  int redlist_size,chain_size,end_size;
  int *redlist,*nchain,*end;
  int **chain;
  double cutoff, cutoffsq;
  double angstrom, angstromrec, qelectron, qelectronrec, forceunit;
  double sigma, epsilon, n_sigma, radius, radiussq, diameter, rc, rc0, comega, ctheta;
  double start_gamma, start_uinf, startxi_usemi, starth_phi;
  double del_gamma, del_uinf, delxi_usemi, delh_phi;
  double *starth_usemi, *startzeta_phi;
  double *delh_usemi, *delzeta_phi;
  double *gamma_data, *uinf_data;
  double *p1, *p2, *param;
  double **usemi_data, **phi_data;
  double **gamma_coeff, **uinf_coeff;
  double **flocal,**basis;
  double ***usemi_coeff, ***phi_coeff;
  char *gamma_file, *uinf_file, *usemi_file, *phi_file;
  
  void allocate();
  
  double spline(double, double, double, double **, int);
  double spline(double, double, double *, double, double *, 
		  double, double ***, int);
  double dspline(double, double, double, double **, int);
  double dxspline(double, double, double *, double, double *, 
		  double, double ***, int);
  double dyspline(double, double, double *, double, double *, 
		  double, double ***, int);

  void spline_coeff(double *, double **, int);
  void spline_coeff(double **, double ***, int);
  
  void read_file(char *, double *, double *, double *, int);
  void read_file(char *, double **, double *, double *, 
		  double *, double *, int);

  double uinf(double *);
  double usemi(double *);
  void finf(double *, double **);
  void fsemi(double *, double **);

  void geom(const double *, const double *, const double *, 
		  const double *, double *, double **);
  double weight(const double *, const double *,
		  const double *, const double *);
  int heaviside(double);
  double s(double x);
};

}
#endif
#endif
