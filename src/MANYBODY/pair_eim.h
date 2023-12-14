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
PairStyle(eim,PairEIM);
// clang-format on
#else

#ifndef LMP_PAIR_EIM_H
#define LMP_PAIR_EIM_H

#include "pair.h"

#include <map>
#include <utility>

namespace LAMMPS_NS {

class PairEIM : public Pair {
 public:
  PairEIM(class LAMMPS *);
  ~PairEIM() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

  struct Setfl {
    double division, rbig, rsmall;
    int nr;
    int *ielement, *tp;
    double *mass, *negativity, *ra, *ri, *Ec, *q0;
    double *rcutphiA, *rcutphiR, *Eb, *r0, *alpha, *beta, *rcutq, *Asigma, *rq, *rcutsigma, *Ac,
        *zeta, *rs;
    double dr, cut;
    double ***Fij, ***Gij, ***phiij;
    double **cuts;
  };

 protected:
  double **cutforcesq, cutmax;
  int nmax;
  double *rho, *fp;
  int rhofp;

  Setfl *setfl;

  // potentials as array data

  int nr;
  int nFij, nGij, nphiij;
  double **Fij, **Gij, **phiij;
  int **type2Fij, **type2Gij, **type2phiij;

  // potentials in spline form used for force computation

  double dr, rdr;
  double *negativity, *q0;
  double ***Fij_spline, ***Gij_spline, ***phiij_spline;

  void allocate();
  void array2spline();
  void interpolate(int, double, double *, double **, double);

  double funccutoff(double, double, double);
  double funcphi(int, int, double);
  double funcsigma(int, int, double);
  double funccoul(int, int, double);

  void read_file(char *);
  void deallocate_setfl();
  void file2array();
};

class EIMPotentialFileReader : protected Pointers {
  std::string filename;
  static const int MAXLINE = 1024;
  char line[MAXLINE];
  double conversion_factor;

  void parse(FILE *fp);
  char *next_line(FILE *fp);
  std::pair<std::string, std::string> get_pair(const std::string &a, const std::string &b);

 public:
  EIMPotentialFileReader(class LAMMPS *lmp, const std::string &filename,
                         const int auto_convert = 0);

  void get_global(PairEIM::Setfl *setfl);
  void get_element(PairEIM::Setfl *setfl, int i, const std::string &name);
  void get_pair(PairEIM::Setfl *setfl, int ij, const std::string &elemA, const std::string &elemB);

 private:
  // potential parameters
  double division;
  double rbig;
  double rsmall;

  struct ElementData {
    int ielement;
    double mass;
    double negativity;
    double ra;
    double ri;
    double Ec;
    double q0;
  };

  struct PairData {
    double rcutphiA;
    double rcutphiR;
    double Eb;
    double r0;
    double alpha;
    double beta;
    double rcutq;
    double Asigma;
    double rq;
    double rcutsigma;
    double Ac;
    double zeta;
    double rs;
    int tp;
  };

  std::map<std::string, ElementData> elements;
  std::map<std::pair<std::string, std::string>, PairData> pairs;
};

}    // namespace LAMMPS_NS

#endif
#endif
