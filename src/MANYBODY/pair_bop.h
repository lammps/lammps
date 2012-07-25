/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   The this work follows the formulation from (a) D.G. Pettifor, et al., Mat.
   Sci. and Eng. A365, 2-13, (2004) and (b) D.A. Murdick, et al., Phys.
   Rev. B 73, 045206 (2006). (c) D.K. Ward, et al., Phys. Rev. B 85, 115206
   (2012)

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bop,PairBOP)

#else

#ifndef LMP_PAIR_BOP_H
#define LMP_PAIR_BOP_H

#include "pair.h"
#include "update.h"

namespace LAMMPS_NS {

class PairBOP : public Pair {
 public:
  PairBOP(class LAMMPS *);
  virtual ~PairBOP();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 private:
  int me;
  int maxneigh;                    // maximum size of neighbor list on this processor
  int update_list;                 // check for changing maximum size of neighbor list
  int maxbopn;                     // maximum size of bop neighbor list for allocation
  int maxnall;                     // maximum size of bop neighbor list for allocation
  int *map;                        // mapping from atom types to elements
  int nelements;                   // # of unique elments
  int nr;                     // increments for the BOP potential
  int nBOt;                   // second BO increments
  int bop_types;                   // number of elments in potential
  int npairs;                      // number of element pairs
  char **elements;                   // names of unique elements
  int ***elem2param;
  int nparams;
  int bop_step;
  int allocate_pi;
  int allocate_sigma;
  int allocate_neigh;
  int nb_pi,nb_sg;

  int *BOP_index;                 // index for neighbor list position
  int neigh_total;                // total number of neighbors stored
  int *cos_index;                 // index for neighbor cosine if not using on the fly
  int *neigh_flag;                 // index for neighbor cosine if not using on the fly
  int cos_total;                  // number of cosines stored if not using on the fly
  int neigh_ct;                  //  limit for large arrays

/*Parameters variables*/

  int ncutoff,nfunc;
  int a_flag;
  double *pi_a,*pro_delta,*pi_delta;
  double *pi_p,*pi_c,*sigma_r0,*pi_r0,*phi_r0;
  double *sigma_rc,*pi_rc,*phi_rc,*r1,*sigma_beta0;
  double *pi_beta0,*phi0,*sigma_n,*pi_n,*phi_m;
  double *sigma_nc,*pi_nc,*phi_nc;
  double *pro,*sigma_delta,*sigma_c,*sigma_a;
  double ***sigma_g0,***sigma_g1,***sigma_g2,***sigma_g3;
  double ***sigma_g4,*sigma_f,*sigma_k,*small3;
  double small1,small2,small3g,small4,small5,small6,small7;
  double which,alpha,alpha1,beta1,gamma1,alpha2,beta2,alpha3;
  double beta3,rsmall,rbig,rcore;
  char **words;

  double cutmax;                //max cutoff for all elements
  int otfly;                        //Defines whether to do on the fly
                                //calculations of angles and distances
                                //on the fly will slow down calculations
                                //but requires less memory on = 1, off=0

  int table;                        //determines the method for reading in
                                //potential parameters a preset table
                                //or generate the tables using a spline

/* Neigh variables */

  double *rcut,*dr,*rdr;
  double **disij,*rij;

/*Triple variables */

  double *cosAng,***dcosAng,***dcAng;

/*Double variables */

  double *betaS,*dBetaS,*betaP;
  double *dBetaP,*repul,*dRepul;

/*Sigma variables */

  int **itypeSigBk,*nSigBk;
  double *sigB;
  double *sigB1;


/*Pi variables */

  int **itypePiBk,*nPiBk;
  double *piB;

/*Grids1 variables */

  double **pBetaS,**pBetaS1,**pBetaS2,**pBetaS3;
  double **pBetaS4,**pBetaS5,**pBetaS6;

/*Grids2 variables */

  double **pBetaP,**pBetaP1,**pBetaP2,**pBetaP3;
  double **pBetaP4,**pBetaP5,**pBetaP6;

/*Grids3 variables */

  double **pRepul,**pRepul1,**pRepul2,**pRepul3;
  double **pRepul4,**pRepul5,**pRepul6;

/*Grids4 variables */

  double **FsigBO,**FsigBO1,**FsigBO2,**FsigBO3;
  double **FsigBO4,**FsigBO5,**FsigBO6;
  double dBO,rdBO;

/* End of BOP variables */

  double **rcmin,**rcmax,**rcmaxp;
  struct B_PI{
    double dAA[3];
    double dBB[3];
    double dPiB[3];
    int temp;
    int i;
    int j;
  };
  B_PI *bt_pi;

  struct B_SG{
    double dAA[3];
    double dBB[3];
    double dCC[3];
    double dDD[3];
    double dEE[3];
    double dEE1[3];
    double dFF[3];
    double dAAC[3];
    double dBBC[3];
    double dCCC[3];
    double dDDC[3];
    double dEEC[3];
    double dFFC[3];
    double dGGC[3];
    double dUT[3];
    double dSigB1[3];
    double dSigB[3];
    int temp;
    int i;
    int j;
  };
  B_SG *bt_sg;

  void setPbetaS();
  void setPbetaP();
  void setPrepul();
  void setSign();
  void gneigh();
  void theta();
  void theta_mod();
  void sigmaBo();
  void PiBo();
  void sigmaBo_otf();
  void PiBo_otf();
  void sigmaBo_noa();
  void sigmaBo_noa_otf();
  void memory_theta_create();
  void memory_theta_destroy();
  void memory_theta_grow();
  double cutoff(double, double, int, double);
  double betaSfunc(int, double);
  double dBetaSfunc(int, double, double, double);
  double betaPfunc(int, double);
  double dBetaPfunc(int, double, double, double);
  double repulfunc(int, double);
  double dRepulfunc(int, double, double, double);

  void read_file(char *);
  void read_table(char *);
  void allocate();
  void create_pi(int);
  void create_sigma(int);
  void destroy_pi();
  void destroy_sigma();
  void grow_pi(int,int);
  void grow_sigma(int,int);
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

E: Pair style BOP requires atom IDs

This is a requirement to use the BOP potential.

E: Pair style BOP requires newton pair on

See the newton command.  This is a restriction to use the BOP
potential.

E: Pair style bop requires comm ghost cutoff at least 3x larger than %g

Use the communicate ghost command to set this.  See the pair bop
doc page for more details.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Too many atom pairs for pair bop

The number of atomic pairs exceeds the expected number.  Check your
atomic structure to ensure that it is realistic.

E: Too many atom triplets for pair bop

The number of three atom groups for angle determinations exceeds the
expected number.  Check your atomic structrure to ensure that it is
realistic.

E: Cannot open BOP potential file %s

The specified BOP potential file cannot be opened.  Check that the
path and name are correct.

*/
