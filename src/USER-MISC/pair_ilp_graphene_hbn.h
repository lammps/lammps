/* ----------------------------------------------------------------------
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

PairStyle(ilp/graphene/hbn,PairILPGrapheneHBN)

#else

#ifndef LMP_PAIR_ILP_GRAPHENE_HBN_H
#define LMP_PAIR_ILP_GRAPHENE_HBN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBN : public Pair {
 public:
  PairILPGrapheneHBN(class LAMMPS *);
  virtual ~PairILPGrapheneHBN();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void ILP_neigh();
  void calc_normal();
  void calc_FRep(int, int);
  void calc_FvdW(int, int);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  int me;
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  MyPage<int> *ipage;              // neighbor list pages
  int *ILP_numneigh;                // # of pair neighbors for each atom
  int **ILP_firstneigh;             // ptr to 1st neighbor of each atom
  int tap_flag;                    // flag to turn on/off taper function

  struct Param {
    double z0,alpha,epsilon,C,delta,d,sR,reff,C6,S;
    double delta2inv,seff,lambda,rcut;
    int ielement,jelement;
  };
  Param *params;       // parameter set for I-J interactions
  char **elements;     // names of unique elements
  int **elem2param;    // mapping from element pairs to parameters
  int *map;            // mapping from atom types to elements
  int nelements;       // # of unique elements
  int nparams;         // # of stored parameter sets
  int maxparam;        // max # of parameter sets
  int nmax;            // max # of atoms

  double cut_global;
  double cut_normal;
  double **cut;
  double **cutILPsq;    // mapping the cutoff from element pairs to parameters
  double **offset;
  double **normal;
  double ***dnormdri;
  double ****dnormal;

  void read_file( char * );
  void allocate();

  /* ----Calculate the long-range cutoff term */
  inline double calc_Tap(double r_ij, double Rcut) {
    double Tap,r;
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    if(r >= 1.0) {Tap = 0.0;}
    else {
      Tap = Tap_coeff[7] * r + Tap_coeff[6];
      Tap = Tap * r  + Tap_coeff[5];
      Tap = Tap * r  + Tap_coeff[4];
      Tap = Tap * r  + Tap_coeff[3];
      Tap = Tap * r  + Tap_coeff[2];
      Tap = Tap * r  + Tap_coeff[1];
      Tap = Tap * r  + Tap_coeff[0];
    }
    return(Tap);
}

  /* ----Calculate the derivatives of long-range cutoff term */
  inline double calc_dTap(double r_ij, double Rcut) {
    double dTap,r;
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    if(r >= 1.0) {dTap = 0.0;}
    else {
      dTap = 7.0*Tap_coeff[7] * r + 6.0*Tap_coeff[6];
      dTap = dTap * r  + 5.0*Tap_coeff[5];
      dTap = dTap * r  + 4.0*Tap_coeff[4];
      dTap = dTap * r  + 3.0*Tap_coeff[3];
      dTap = dTap * r  + 2.0*Tap_coeff[2];
      dTap = dTap * r  + Tap_coeff[1];
      dTap = dTap/Rcut;
    }
    return(dTap);
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

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/

