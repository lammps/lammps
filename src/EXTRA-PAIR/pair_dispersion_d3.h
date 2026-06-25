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
#include <cmath>
#include <type_traits>

namespace LAMMPS_NS {
template <typename T> inline T pow2(T x)
{ return x * x; }

template <typename T> inline T pow4(T x)
{
  const T x2 = x * x;
  return x2 * x2;
}

template <typename T> inline T pow6(T x)
{ return pow2(x) * pow4(x); }

template <typename T> inline T pow8(T x)
{
  const T x4 = pow4(x);
  return x4 * x4;
}

template <typename T> inline T pow_int(T base, int exponent)
{
  if (exponent == 0) return static_cast<T>(1);

  if (exponent < 0) return static_cast<T>(1) / pow_int(base, -exponent);

  T result = static_cast<T>(1);
  while (exponent > 0) {
    if (exponent & 1) result *= base;
    base *= base;
    exponent >>= 1;
  }
  return result;
}

template <typename T> inline bool is_integer_value(T x)
{
  const T nearest_int = std::round(x);
  return x == nearest_int;
}

template <typename T> inline T pow_general(T base, T exponent)
{
  if (is_integer_value(exponent)) return pow_int(base, static_cast<int>(std::round(exponent)));

  using std::pow;
  return pow(base, exponent);
}

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
