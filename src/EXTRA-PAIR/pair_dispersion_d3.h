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

// global ad hoc parameters of the D3 model, shared by the plain pair style and
// its accelerated variants.  They live in a named namespace, not at global
// scope, because this is a style header and gets included all over the place.

namespace LAMMPS_NS::DispersionD3 {

  static constexpr double K1 = 16.0;

  /*  reasonable choices for k3 are between 3 and 5 :
      this gives smooth curves with maxima around the integer values
      k3=3 give for CN=0 a slightly smaller value than computed
      for the free atom. This also yields to larger CN for atoms
      in larger molecules but with the same chemical environment
      which is physically not right.
      values >5 might lead to bumps in the potential.
  */

  static constexpr double K3 = -4.0;

  static constexpr double AUTOANG = 0.52917725;     // atomic units (Bohr) to Angstrom
  static constexpr double AUTOEV = 27.21140795;     // atomic units (Hartree) to eV

  // conversion factor for the tabulated C6 reference values
  static constexpr double AUTOANG3 = AUTOANG * AUTOANG * AUTOANG;
  static constexpr double AUTOANG6 = AUTOANG3 * AUTOANG3;
}    // namespace LAMMPS_NS::DispersionD3

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

  int dampingCode;                         // Which damping function to use
  double s6, s8, s18, rs6, rs8, rs18;      // XC parameters
  double a1, a2, alpha, alpha6, alpha8;    // XC parameters

  double *r2r4;        // scale r4/r2 values of the atoms by sqrt(Z)
  double *rcov;        // covalent radii
  int *mxci;           // How large the grid for c6 interpolation
  double **r0ab;       // cut-off radii for all element pairs
  double *****c6ab;    // C6 for all element pairs
  int max_mxci;        // Maximum grid size of the C_i coefficient
  double *cn;          // Coordination numbers
  double *dc6;         // dC6i(iat) saves dE_dsp/dCN(iat)

  int communicationStage;    // communication stage

  double memory_usage() override;
  virtual void allocate();
  virtual void set_funcpar(std::string &);

  virtual void calc_coordination_number();

  int find_atomic_number(std::string &);
  std::vector<int> is_int_in_array(int *, int, int);

  void read_r0ab(int *, int);
  void set_limit_in_pars_array(int &, int &, int &, int &);
  void read_c6ab(int *, int);

  // writes {C6, dC6/dCN_i, dC6/dCN_j} to c6_res; must not use static storage,
  // it is called concurrently from the threaded and device variants
  void get_dC6(int, int, double, double, double *);
};
}    // namespace LAMMPS_NS
#endif
#endif
