/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(pinn, PairPINN)

#else

#ifndef LMP_PAIR_PINN_H
#define LMP_PAIR_PINN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPINN : public Pair {
  public:
  PairPINN(class LAMMPS *lmp);
  ~PairPINN();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  struct Layer {
    int nnodes;
    double *Weights; // a matrix to store weights.
    double *Biases;  // a vector to store biases.
    double *fdot;    // a vector to store derivatives of logistic function w.r.t. their arguments.
  };

  protected:
  struct Param {
    int nlayers;                  // # of layers of neural network
    Layer *layers;                // layers of neural network
    double cut;                   // cutoff (rc)
    double cutsq;                 // square of cutoff
    double cut_range;             // cutoff range (hc)
    double giref;                 // shift to Gis
    int gimethod;                 // method to compute Gis
    int actfunc;                  // activation function to use
    double gwidth;                // Gaussian width used in Gis
    int nLPOrders;                // # of orders of Legendre Polynomials
    int *LPOrders;                // Legendre Polynomial orders
    int nsigmas;                  // # of positions of gaussian functions
    double *sigmas;               // positions of gaussians
    double *mass;
  };

  Param *params;
  int *map;                      // mapping from atom types to elements
  int **map_gi;                  // mapping from element pairs, j-k, to gis
  int nelements;                 // # of unique elements
  char **elements;               // names of unique elements
  double cutmax;
  double cutmaxsq;
  double *baseline_bo_params;    // baseline BO parameters
  int nbaseline_bo_params;       // # of baseline BO parameters
  int nBase;                     // flag: 0 - base BOP parameters are not used
                                 //       1 - base BOP parameters are unsed
  int pinn_kind_flag1;
  int pinn_kind_flag2;

  // bond-order parameters
  double **big_a;       // A_ij
  double **alpha;       // alpha_ij
  double **big_b;       // B_ij
  double **beta;        // beta_ij
  double ***small_h;    // h_ijk
  double *sigma;        // sigma_i
  double ***small_a;    // a_ijk
  double ***lambda;     // lamda_ijk

  // derivatives wrt bond-order parameters
  double **deriv_wrt_big_a;
  double **deriv_wrt_alpha;
  double **deriv_wrt_big_b;
  double **deriv_wrt_beta;
  double ***deriv_wrt_small_h;
  double *deriv_wrt_sigma;
  double ***deriv_wrt_small_a;
  double ***deriv_wrt_lambda;

  // specifics of BO
  double *Sij;
  double *Zij;
  double *bij;
  double **Sijk;
  double **CSijk;

  void allocate();
  void read_file(char *);
  void setup_params();

  void eval_nnet(double *, int, double *, int);
  void eval_nnet_d(Layer *, double *, int, double *, int);
  double ann_fc(double, Param *);
  double ann_fc_d(double, Param *);
  double ann_fc(double, double, double);
  double ann_fc_d(double, double, double);
  double ann_fs(double, int, Param *);
  double pinn_fs(const double, const int);
  double ann_fs_d(double, int, Param *);
  double pinn_fs_d(const double, const int);
  void compute_atomic_LSP(const int, double *);

  double bop_fr(int, int, double);
  double bop_fr_d(int, int, Param *, double);
  double bop_fa(int, int, double);
  double bop_fa_d(int, int, Param *, double);
  double screen_ijk(int, int, int, Param *, double);
  double bop_zeta(int, int, int, Param *, double, double, double);
  double dot_product(double *, double *, int);
  void force_Sijk_ri(double, double *, double *, double, double, double,
                     double, double *, double *, double *);
  void dcostheta_ijk_dri(double cstheta, double *drj,
                         double *drk, double *fj, double *fk);

  void unpack_vec_to_bop_sets(const double *, const int);

  void force_costheta_ijk_ri(double, double, double *, double *,
                             double *, double *, double *);
  void iderivs_wrt_bops(const int, const double);
  void pack_derivs_wrt_bop_sets_to_vec(double *, const int);

  void vec_mult_mat(double *, double *, double *, const int, const int);

};

}

#endif
#endif

