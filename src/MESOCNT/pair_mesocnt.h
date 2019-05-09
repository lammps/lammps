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
  void write_restart(FILE *);
  void read_restart(FILE *);

  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  int n, gamma_points, pot_points;
  double sigma, epsilon, n_sigma, radius;
  double start_gamma, start_uinf, startxi_usemi, starth_phi;
  double del_gamma, del_uinf, delxi_usemi, delh_phi;
  double *delh_usemi, *delzeta_phi;
  double *gamma_data, *uinf_data;
  double **usemi_data, **phi_data;
  double **gamma_coeff, **uinf_coeff;
  double ***usemi_coeff, ***phi_coeff;
  char *gamma_file, *uinf_file, *usemi_file, *phi_file;
  
  void allocate();
  
  double spline(double, double, double, double **, int);
  double spline(double, double, double *, double, double *, double ***, int);
  double dspline(double, double, double, double **, int);
  double dxspline(double, double, double *, double, double *, double ***, int);
  double dyspline(double, double, double *, double, double *, double ***, int);

  void spline_coeff(double *, double **, int);
  void spline_coeff(double **, double ***, int);
  
  void read_file(char *, double *, double *, int);
  void read_file(char *, double **, double *, double *, int);
};

}
#endif
#endif
