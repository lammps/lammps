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
  double del_gamma, del_u_inf, delxi_u_semi, delh_phi;
  double *delh_u_semi, *delzeta_phi;
  double *gamma_data, *u_inf_data;
  double **u_semi_data, **phi_data;
  double **gamma_coeff, **u_inf_coeff;
  double ***u_semi_coeff, ***phi_coeff;
  char *gamma_file, *u_inf_file, *u_semi_file, *phi_file;
  
  void allocate();
  
  double eval_gamma(double);
  double eval_u_inf(double);
  double eval_u_semi(double, double);
  double eval_u_semi(double, double);
  double d_gamma(double);
  double d_u_inf(double);
  double dh_u_semi(double, double);
  double dxi_u_semi(double, double);
  double dh_u_semi(double, double);
  double dzeta_u_semi(double, double);

  int heaviside(double);
  int sgn(double);

  double distance(double, double, double, double, double, double, double);
  double cutoff(double, double, double);

  void read_file(char *, double *);
  void read_file(char *, double **);
}

}
#endif
#endif
