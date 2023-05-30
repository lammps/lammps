/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
// clang-format off
PairStyle(bop,PairBOP);
// clang-format on
#else

#ifndef LMP_PAIR_BOP_H
#define LMP_PAIR_BOP_H

#include "pair.h"

namespace LAMMPS_NS {
class TabularFunction;

class PairBOP : public Pair {

 public:
  PairBOP(class LAMMPS *);
  ~PairBOP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

 private:
  struct PairParameters {
    double cutB, cutBsq, cutL, cutLsq;
    TabularFunction *betaS;
    TabularFunction *betaP;
    TabularFunction *rep;
    TabularFunction *cphi;
    TabularFunction *bo;
    PairParameters();
    ~PairParameters();
  };

  struct PairList1 {
    double r, dis[3];
    double betaS, dBetaS, betaP, dBetaP, rep, dRep;
    PairList1(){};
  };

  struct PairList2 {
    double r, dis[3];
    double rep, dRep;
    PairList2(){};
  };

  struct TripleList {
    double G, dG, cosAng, dCosAngi[3], dCosAngj[3], dCosAngk[3];
    TripleList(){};
  };

  struct B_SG {
    double dAA[3];
    double dBB[3];
    double dCC[3];
    double dDD[3];
    double dEE1[3];
    double dFF[3];
    double dAAC[3];
    double dSigB1[3];
    double dSigB[3];
    int temp;
    int i;
    int j;
  };

  struct B_PI {
    double dAA[3];
    double dBB[3];
    double dPiB[3];
    int temp;
    int i;
    int j;
  };

  PairParameters *pairParameters;
  TabularFunction *tripletParameters;

  // Parameters variables

  double small1, small2, small3g, small4, small5, small6, small7, *pi_p;
  double *sigma_c, *sigma_a, *pi_c, *pi_a, *sigma_delta, *pi_delta;
  double *sigma_f, *sigma_k, *small3;
  double *pro_delta, *pro;

  int bop_types;          // number of elments in potential file
  int npairs;             // number of element pairs
  int ntriples;           // number of all triples
  char **bop_elements;    // names of elements in potential file
  double *bop_masses;     // masses of elements in potential file
  double bytes;

  int otfly;    // = 1 faster, more memory, = 0 slower, less memory

  PairList1 *pairlist1;
  PairList2 *pairlist2;
  TripleList *triplelist;

  B_SG *bt_sg;
  B_PI *bt_pi;

  int *BOP_index;       // index for neighbor list position
  int *BOP_total;       // index for neighbor list position
  int *BOP_index2;      // index for neighbor list position
  int *BOP_total2;      // index for neighbor list position
  int *neigh_index;     // index for neighbor list position
  int *neigh_index2;    // index for neighbor list position
  int atomlimit;        // current size of atom based list
  int neighlimit;       // current size of neighbor based list
  int neighlimit2;      // current size of neighbor based list
  int neineilimit;      // current size of triple based list
  int sglimit;          // current size of bt_sg
  int pilimit;          // current size of bt_pi
  int *cos_index;       // index for neighbor cosine if not using on the fly
  double cutmax;

#if defined(LMP_BOP_WRITE_TABLES)
  void write_tables(int);
#endif

  void gneigh();
  void angle(double, double *, double, double *, double &, double *, double *);
  double SigmaBo(int, int);
  double PiBo(int, int);
  void read_table(char *);
  void allocate();
  void memory_sg(int);
  void memory_pi(int);
  void initial_sg(int);
  void initial_pi(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
