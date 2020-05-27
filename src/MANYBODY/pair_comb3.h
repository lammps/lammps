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

PairStyle(comb3,PairComb3)

#else

#ifndef LMP_PAIR_COMB3_H
#define LMP_PAIR_COMB3_H

#include "pair.h"

namespace LAMMPS_NS {

class PairComb3 : public Pair {
 public:
  PairComb3(class LAMMPS *);
  virtual ~PairComb3();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();
  virtual double combqeq(double *, int &);
  double enegtot;

 // general potential parameters
 protected:
  struct Param {
    int ielement,jelement,kelement,powermint;
    int ielementgp,jelementgp,kelementgp;       //element group
    int ang_flag,pcn_flag,rad_flag,tor_flag;    //angle, coordination,radical, torsion flag
    double lami,lambda,alfi,alpha1,alpha2,alpha3,beta;
    double pcos6,pcos5,pcos4,pcos3,pcos2,pcos1,pcos0;
    double gamma,powerm,powern,bigA,bigB1,bigB2,bigB3;
    double bigd,bigr,cut,cutsq,c1,c2,c3,c4;
    double p6p0,p6p1,p6p2,p6p3,p6p4,p6p5,p6p6;
    double ptork1,ptork2;
    double addrepr,addrep, vdwflag;
    double QU,QL,DU,DL,Qo,dQ,aB,bB,nD,bD,qmin,qmax;
    double chi,dj,dk,dl,dm,esm,cmn1,cmn2,pcmn1,pcmn2;
    double coulcut, lcut, lcutsq;
    double veps, vsig, pcna, pcnb, pcnc, pcnd, polz, curl, pcross;
    double paaa, pbbb;
    double curlcut1, curlcut2, curl0;
  };

  // general setups
  int nelements;                        // # of unique elements
  int ***elem2param;                    // mapping from element triplets to parameters
  int *map;                             // mapping from atom types to elements
  int nparams;                          // # of stored parameter sets
  int maxparam;                         // max # of parameter sets
  double PI,PI2,PI4,PIsq;               // PIs
  double cutmin;                        // min cutoff for all elements
  double cutmax;                        // max cutoff for all elements
  double precision;                     // tolerance for QEq convergence
  char **elements;                      // names of unique elements
  Param *params;                        // parameter set for an I-J-K interaction
  int debug_eng1, debug_eng2, debug_fq; // logic controlling debugging outputs
  int pack_flag;

  // Short range neighbor list
  void Short_neigh();
  int pgsize, oneatom;
  int *sht_num, **sht_first;
  MyPage<int> *ipage;

  // loop up tables and flags
  int nmax, **intype;
  int  pol_flag, polar;
  double *qf, **bbij, *charge, *NCo;
  double *esm, **fafb, **dfafb, **ddfafb, **phin, **dphin, **erpaw;
  double **vvdw, **vdvdw;
  double **afb, **dafb;
  double **dpl, bytes;
  double *xcctmp, *xchtmp, *xcotmp;

  // additional carbon parameters
  int cflag;
  int nsplpcn,nsplrad,nspltor;
  int maxx,maxy,maxz,maxxc,maxyc,maxconj;
  int maxxcn[4];
  double vmaxxcn[4],dvmaxxcn[4];
  int ntab;
  double iin2[16][2],iin3[64][3];
  double brad[4], btor[4], bbtor, ptorr;
  double fi_tor[3], fj_tor[3], fk_tor[3], fl_tor[3];
  double radtmp, fi_rad[3], fj_rad[3], fk_rad[3];

  double ccutoff[6],ch_a[7];

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
  void attractive(Param *, Param *, Param *, double, double, double, double,
        double, double, double, double *, double *, double *,
        double *, double *, int, double);
  virtual void comb_fa(double, Param *, Param *, double, double,
          double &, double &);
  virtual void repulsive(Param *, Param *,double, double &, int,
         double &, double, double);

  // bond order terms
  double comb_bij(double, Param *, double, int, double);
  double comb_gijk(double, Param *, double);
  void comb_gijk_d(double, Param *, double, double &, double &);
  double zeta(Param *, Param *, double, double, double *, double *, int, double);
  void comb_bij_d(double, Param *, double, int, double &,
          double &, double &, double &, double &, double &, double);
  void coord(Param *, double, int, double &, double &,
          double &, double &, double &, double);
  void comb_zetaterm_d(double, double, double, double, double,
        double *, double, double *, double, double *, double *,
        double *, Param *, Param *, Param *, double);
  void costheta_d(double *, double, double *, double,
          double *, double *, double *);
  void force_zeta(Param *, Param *, double, double, double, double &,
        double &, double &, double &, double &, double &, double &,
        double &, double &, double &, double &, double &, double &,
        int, double &, double,double, int, int, int,
        double , double , double);
  void cntri_int(int, double, double, double, int, int, int,
        double &, double &, double &, double &, Param *);

  // Legendre polynomials
  void selfp6p(Param *, Param *, double, double &, double &);
  double ep6p(Param *, Param *, double, double, double *, double * ,double &);
  void fp6p(Param *, Param *, double, double, double *, double *, double *,
          double *, double *);

  // long range q-dependent terms
  double self(Param *, double);
  void tables();
  void potal_calc(double &, double &, double &);
  void tri_point(double, int &, int &, int &, double &, double &,
         double &);
  void vdwaals(int,int,int,int,double,double,double,double,
          double &, double &);
  void direct(Param *, Param *, int,int,int,double,double,
        double,double,double,double, double,double,double &,double &,
         int, int);
  void field(Param *, Param *,double,double,double,double &,double &);
  int heaviside(double);
  double switching(double);
  double switching_d(double);
  double chicut1, chicut2;

  // radical terms
  double rad_init(double, Param *, int, double &, double);
  void rad_calc(double, Param *, Param *, double, double, int,
                  int, double, double);
  void rad_int(int , double, double, double, int, int, int,
        double &, double &, double &, double &);
  void rad_forceik(Param *,  double, double *, double, double);
  void rad_force(Param *,  double, double *,  double);

  // torsion terms
  double bbtor1(int, Param *, Param *, double, double, double,
        double *, double *, double *, double);                     //modified by TAO
  void tor_calc(double, Param *, Param *, double, double, int,
                  int, double, double);
  void tor_int(int , double, double, double, int, int, int,
        double &, double &, double &, double &);
  void tor_force(int, Param *, Param *, double, double, double,
        double, double *, double *, double *);                //modified by TAO

  // charge force terms
  double qfo_self(Param *, double);
  void qfo_short(Param *, Param *, double, double, double,
        double &, double &, int, int, int);
  void qfo_direct(Param *, Param *, int, int, int, double,
        double, double, double, double, double &, double &,
        double, double, int, int);
  void qfo_field(Param *, Param *,double,double ,double ,double &, double &);
  void qfo_dipole(double, int, int, int, int, double, double *, double,
        double, double, double &, double &, int, int);
  void qsolve(double *);

  // dipole - polarization terms
  double dipole_self(Param *, int);
  void dipole_init(Param *, Param *, double, double *, double,
        int, int, int, double, double, double, double, double, int , int);
  void dipole_calc(Param *, Param *, double, double, double, double, double,
        int, int, int, double, double, double, double, double, int , int,
        double &, double &, double *);

  // communication functions
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  // vector functions, inline for efficiency
  inline double vec3_dot(double *x, double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }
  inline void vec3_add(double *x, double *y, double *z) {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }
  inline void vec3_scale(double k, double *x, double *y) {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }
  inline void vec3_scaleadd(double k, double *x, double *y, double *z) {
    z[0] = k*x[0]+y[0];  z[1] = k*x[1]+y[1];  z[2] = k*x[2]+y[2];
  }
 };
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style COMB3 requires atom IDs

This is a requirement to use the COMB3 potential.

E: Pair style COMB3 requires newton pair on

See the newton command.  This is a restriction to use the COMB3
potential.

E: Pair style COMB3 requires atom attribute q

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open COMB3 lib.comb3 file

The COMB3 library file cannot be opened.  Check that the path and name
are correct.

E: Cannot open COMB3 potential file %s

The specified COMB3 potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in COMB3 potential file

Incorrect number of words per line in the potential file.

E: Illegal COMB3 parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

E: Neighbor list overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

E: Error in vdw spline: inner radius > outer radius

A pre-tabulated spline is invalid.  Likely a problem with the
potential parameters.

*/
