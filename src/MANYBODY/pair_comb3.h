/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(comb3,PairComb3);
// clang-format on
#else

#ifndef LMP_PAIR_COMB3_H
#define LMP_PAIR_COMB3_H

#include "pair.h"

namespace LAMMPS_NS {

class PairComb3 : public Pair {
 public:
  PairComb3(class LAMMPS *);
  ~PairComb3() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  virtual double combqeq(double *, int &);
  double enegtot;

  static constexpr int NPARAMS_PER_LINE = 74;

 protected:
  // general potential parameters
  struct Param {
    int ielement, jelement, kelement, powermint;
    int ielementgp, jelementgp, kelementgp;        //element group
    int ang_flag, pcn_flag, rad_flag, tor_flag;    //angle, coordination,radical, torsion flag
    double lami, lambda, alfi, alpha1, alpha2, alpha3, beta;
    double pcos6, pcos5, pcos4, pcos3, pcos2, pcos1, pcos0;
    double gamma, powerm, powern, bigA, bigB1, bigB2, bigB3;
    double bigd, bigr, cut, cutsq, c1, c2, c3, c4;
    double p6p0, p6p1, p6p2, p6p3, p6p4, p6p5, p6p6;
    double ptork1, ptork2;
    double addrepr, addrep, vdwflag;
    double QU, QL, DU, DL, Qo, dQ, aB, bB, nD, bD, qmin, qmax;
    double chi, dj, dk, dl, dm, esm, cmn1, cmn2, pcmn1, pcmn2;
    double coulcut, lcut, lcutsq;
    double veps, vsig, pcna, pcnb, pcnc, pcnd, polz, curl, pcross;
    double paaa, pbbb;
    double curlcut1, curlcut2, curl0;
  };

  // general setups
  double PI, PI2, PI4, PIsq;               // PIs
  double cutmin;                           // min cutoff for all elements
  double cutmax;                           // max cutoff for all elements
  double precision;                        // tolerance for QEq convergence
  Param *params;                           // parameter set for an I-J-K interaction
  int debug_eng1, debug_eng2, debug_fq;    // logic controlling debugging outputs
  int pack_flag;

  // Short range neighbor list
  void Short_neigh();
  int pgsize, oneatom;
  int *sht_num, **sht_first;
  MyPage<int> *ipage;

  // loop up tables and flags
  int nmax, **intype;
  int pol_flag, polar;
  double *qf, **bbij, *charge, *NCo;
  double *esm, **fafb, **dfafb, **ddfafb, **phin, **dphin, **erpaw;
  double **vvdw, **vdvdw;
  double **afb, **dafb;
  double **dpl, bytes;
  double *xcctmp, *xchtmp, *xcotmp;

  // additional carbon parameters
  int cflag;
  int nsplpcn, nsplrad, nspltor;
  int maxx, maxy, maxz, maxxc, maxyc, maxconj;
  int maxxcn[4];
  double vmaxxcn[4], dvmaxxcn[4];
  int ntab;
  double iin2[16][2], iin3[64][3];
  double brad[4], btor[4], bbtor, ptorr;
  double fi_tor[3], fj_tor[3], fk_tor[3], fl_tor[3];
  double radtmp, fi_rad[3], fj_rad[3], fk_rad[3];

  double ccutoff[6], ch_a[7];

  //COMB3-v18 arrays for CHO
  // We wanna dynamic arrays
  // C angle arrays, size = ntab+1
  double pang[20001];
  double dpang[20001];
  double ddpang[20001];

  //coordination spline arrays
  double pcn_grid[4][5][5][5];
  double pcn_gridx[4][5][5][5];
  double pcn_gridy[4][5][5][5];
  double pcn_gridz[4][5][5][5];
  double pcn_cubs[4][4][4][4][64];

  //coordination spline arrays
  double rad_grid[3][5][5][11];
  double rad_gridx[3][5][5][11];
  double rad_gridy[3][5][5][11];
  double rad_gridz[3][5][5][11];
  double rad_spl[3][4][4][10][64];

  //torsion spline arrays
  double tor_grid[1][5][5][11];
  double tor_gridx[1][5][5][11];
  double tor_gridy[1][5][5][11];
  double tor_gridz[1][5][5][11];
  double tor_spl[1][4][4][10][64];

  // initialization functions
  void allocate();
  void read_lib();
  void setup_params();
  virtual void read_file(char *);

  // cutoff functions
  double comb_fc(double, Param *);
  double comb_fc_d(double, Param *);
  double comb_fc_curl(double, Param *);
  double comb_fc_curl_d(double, Param *);
  double comb_fccc(double);
  double comb_fccc_d(double);
  double comb_fcch(double);
  double comb_fcch_d(double);
  double comb_fccch(double);
  double comb_fccch_d(double);
  double comb_fcsw(double);

  // short range terms
  void attractive(Param *, Param *, Param *, double, double, double, double, double, double, double,
                  double *, double *, double *, double *, double *, int, double);
  virtual void comb_fa(double, Param *, Param *, double, double, double &, double &);
  virtual void repulsive(Param *, Param *, double, double &, int, double &, double, double);

  // bond order terms
  double comb_bij(double, Param *, double, int, double);
  double comb_gijk(double, Param *, double);
  void comb_gijk_d(double, Param *, double, double &, double &);
  double zeta(Param *, Param *, double, double, double *, double *, int, double);
  void comb_bij_d(double, Param *, double, int, double &, double &, double &, double &, double &,
                  double &, double);
  void coord(Param *, double, int, double &, double &, double &, double &, double &, double);
  void comb_zetaterm_d(double, double, double, double, double, double *, double, double *, double,
                       double *, double *, double *, Param *, Param *, Param *, double);
  void costheta_d(double *, double, double *, double, double *, double *, double *);
  void force_zeta(Param *, Param *, double, double, double, double &, double &, double &, double &,
                  double &, double &, double &, double &, double &, double &, double &, double &,
                  double &, int, double &, double, double, int, int, int, double, double, double);
  void cntri_int(int, double, double, double, int, int, int, double &, double &, double &, double &,
                 Param *);

  // Legendre polynomials
  void selfp6p(Param *, Param *, double, double &, double &);
  double ep6p(Param *, Param *, double, double, double *, double *, double &);
  void fp6p(Param *, Param *, double, double, double *, double *, double *, double *, double *);

  // long range q-dependent terms
  double self(Param *, double);
  void tables();
  void potal_calc(double &, double &, double &);
  void tri_point(double, int &, int &, int &, double &, double &, double &);
  void vdwaals(int, int, int, int, double, double, double, double, double &, double &);
  void direct(Param *, Param *, int, int, int, double, double, double, double, double, double,
              double, double, double &, double &, int, int);
  void field(Param *, Param *, double, double, double, double &, double &);
  int heaviside(double);
  double switching(double);
  double switching_d(double);
  double chicut1, chicut2;

  // radical terms
  double rad_init(double, Param *, int, double &, double);
  void rad_calc(double, Param *, Param *, double, double, int, int, double, double);
  void rad_int(int, double, double, double, int, int, int, double &, double &, double &, double &);
  void rad_forceik(Param *, double, double *, double, double);
  void rad_force(Param *, double, double *, double);

  // torsion terms
  double bbtor1(int, Param *, Param *, double, double, double, double *, double *, double *,
                double);    //modified by TAO
  void tor_calc(double, Param *, Param *, double, double, int, int, double, double);
  void tor_int(int, double, double, double, int, int, int, double &, double &, double &, double &);
  void tor_force(int, Param *, Param *, double, double, double, double, double *, double *,
                 double *);    //modified by TAO

  // charge force terms
  double qfo_self(Param *, double);
  void qfo_short(Param *, Param *, double, double, double, double &, double &, int, int, int);
  void qfo_direct(Param *, Param *, int, int, int, double, double, double, double, double, double &,
                  double &, double, double, int, int);
  void qfo_field(Param *, Param *, double, double, double, double &, double &);
  void qfo_dipole(double, int, int, int, int, double, double *, double, double, double, double &,
                  double &, int, int);
  void qsolve(double *);

  // dipole - polarization terms
  double dipole_self(Param *, int);
  void dipole_init(Param *, Param *, double, double *, double, int, int, int, double, double,
                   double, double, double, int, int);
  void dipole_calc(Param *, Param *, double, double, double, double, double, int, int, int, double,
                   double, double, double, double, int, int, double &, double &, double *);

  // communication functions
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
};
}    // namespace LAMMPS_NS

#endif
#endif
