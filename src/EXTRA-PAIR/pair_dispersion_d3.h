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
PairStyle(dispersion/d3,PairDispersionD3);
// clang-format on
#else

#ifndef LMP_PAIR_DISPERSION_D3_H
#define LMP_PAIR_DISPERSION_D3_H

#include "pair.h"

namespace LAMMPS_NS {

class PairDispersionD3 : public Pair {

 public:
  PairDispersionD3(class LAMMPS *);
  ~PairDispersionD3() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  int pack_reverse_comm(int, int, double *) override;

  void unpack_forward_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  int nmax;

  double rthr;      // R^2 distance to cutoff for D3_calculation
  double cn_thr;    // R^2 distance to cutoff for CN_calculation

  std::string damping_type;              // damping function type
  double s6, s8, s18, rs6, rs8, rs18;    // XC parameters
  double a1, a2, alpha, alpha6, alpha8;

  double *r2r4;        // scale r4/r2 values of the atoms by sqrt(Z)
  double *rcov;        // covalent radii
  int *mxci;           // How large the grid for c6 interpolation
  double **r0ab;       // cut-off radii for all element pairs
  double *****c6ab;    // C6 for all element pairs
  double *cn;          // Coordination numbers
  double *dc6;         // dC6i(iat) saves dE_dsp/dCN(iat)

  int communicationStage;    // communication stage

  void allocate();
  virtual void set_funcpar(std::string &);

  void calc_coordination_number();

  int find_atomic_number(std::string &);
  std::vector<int> is_int_in_array(int *, int, int);

  void read_r0ab(int *, int);
  void set_limit_in_pars_array(int &, int &, int &, int &);
  void read_c6ab(int *, int);

  double *get_dC6(int, int, double, double);
};
}    // namespace LAMMPS_NS
#endif
#endif
