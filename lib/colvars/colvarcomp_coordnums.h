// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.
//

#ifndef COLVARCOMP_COORDNUM_H
#define COLVARCOMP_COORDNUM_H

#include <memory>

#include "colvar.h"
#include "colvarcomp.h"
#include "colvarmodule.h"


/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::coordnum : public colvar::cvc {
public:

  coordnum();
  virtual ~coordnum();
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

  enum {
    ef_null = 0,
    ef_gradients = 1,
    ef_use_pairlist = (1 << 9),
    ef_rebuild_pairlist = (1 << 10)
  };

  /// Compute the switching function (1-l2**(en/2))/(1-l2**(ed/2)) for a given squared scaled distance
  /// @param[in] l2 Square norm of (Dx/r0x, Dy/r0y, Dz/r0z)
  /// @param[out] dFdl2 Derivative of the result with respect to l2
  /// @param en Numerator exponent
  /// @param ed Denominator exponent
  /// @param pairlist_tol Pairlist tolerance
  template <int flags>
  static cvm::real switching_function(cvm::real const &l2, cvm::real &dFdl2, int en, int ed,
                                      cvm::real pairlist_tol);

  /// Main kernel for the coordination number
  template <int flags>
  static cvm::real compute_pair_coordnum(cvm::rvector const &inv_r0_vec,
                                         cvm::rvector const &inv_r0sq_vec, int en, int ed,
                                         const cvm::real a1x, const cvm::real a1y, const cvm::real a1z,
                                         const cvm::real a2x, const cvm::real a2y, const cvm::real a2z,
                                         cvm::real &g1x, cvm::real &g1y, cvm::real &g1z,
                                         cvm::real &g2x, cvm::real &g2y, cvm::real &g2z,
                                         cvm::real pairlist_tol, cvm::real pairlist_tol_l2_max,
                                         cvm::system_boundary_conditions const &bc);

  /// Workhorse function
  template <bool use_group1_com, bool use_group2_com, int flags> int compute_coordnum();

  /// Workhorse function
  template <bool use_group1_com, bool use_group2_com, int flags> void main_loop();

protected:
  /// First atom group
  cvm::atom_group *group1 = nullptr;
  /// Second atom group
  cvm::atom_group *group2 = nullptr;

  /// Cutoff distances along each dimension
  cvm::rvector r0_vec;

  /// Inverse of r0_vec
  cvm::rvector inv_r0_vec;

  /// Square of inv_r0_vec
  cvm::rvector inv_r0sq_vec;

  /// Set r0_vec and related fields
  void update_cutoffs(cvm::rvector const &r0_vec_i);

  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 12;

  /// The number of pairwise distances being calculated
  size_t num_pairs = 0;

  /// If true, group1 will be treated as a single atom
  bool b_group1_center_only = false;

  /// If true, group2 will be treated as a single atom
  bool b_group2_center_only = false;

  /// Tolerance for the pair list
  cvm::real tolerance = 0.0;

  /// Value of the squared scaled distance (l^2) that matches the given tolerance
  cvm::real tolerance_l2_max = 1.0e20;

  /// Recompute the value of tolerance_l2_max
  void compute_tolerance_l2_max();

  /// Frequency of update of the pair list
  int pairlist_freq = 100;

  /// Pair list
  std::unique_ptr<bool []> pairlist;

};


/// \brief Colvar component: self-coordination number within a group
/// (colvarvalue::type_scalar type, range [0:N*(N-1)/2])
class colvar::selfcoordnum : public colvar::coordnum {
public:

  selfcoordnum();
  virtual void calc_value();
  virtual void calc_gradients();

  /// Workhorse function
  template <int flags> void selfcoordnum_sequential_loop();

  /// Main workhorse function
  template <int flags> int compute_selfcoordnum();
};


/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::groupcoordnum : public colvar::coordnum {
public:
  groupcoordnum();
  virtual ~groupcoordnum() {}
  virtual void calc_value();
  virtual void calc_gradients();
};


/// \brief Colvar component: hydrogen bond, defined as the product of
/// a colvar::coordnum and 1/2*(1-cos((180-ang)/ang_tol))
/// (colvarvalue::type_scalar type, range [0:1])
class colvar::h_bond : public colvar::cvc {
public:
  /// Constructor for atoms already allocated
  h_bond(cvm::atom_group::simple_atom const &acceptor, cvm::atom_group::simple_atom const &donor,
         cvm::real r0, int en, int ed);
  h_bond();
  virtual ~h_bond() {}
  virtual int init(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();

protected:
  /// Cutoff distances along each dimension
  cvm::rvector r0_vec;
  /// Integer exponent of the function numerator
  int en = 6;
  /// Integer exponent of the function denominator
  int ed = 8;
};


template <int flags>
inline cvm::real colvar::coordnum::switching_function(cvm::real const &l2, cvm::real &dFdl2,
                                               int en, int ed,
                                               cvm::real pairlist_tol)
{
  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = cvm::integer_power(l2, en2);
  cvm::real const xd = cvm::integer_power(l2, ed2);
  cvm::real const eps_l2 = 1.0e-7;
  cvm::real const h = l2 - 1.0;
  cvm::real const en2_r = (cvm::real) en2;
  cvm::real const ed2_r = (cvm::real) ed2;
  cvm::real func_no_pairlist;

  if (std::abs(h) < eps_l2) {
    // Order-2 Taylor expansion: c0 + c1*h + c2*h^2
    cvm::real const c0 = en2_r / ed2_r;
    cvm::real const c1 = (en2_r * (en2_r - ed2_r)) / (2.0 * ed2_r);
    cvm::real const c2 = (en2_r * (en2_r - ed2_r) * (2.0 * en2_r - ed2_r - 3.0)) / (12.0 * ed2_r);
    func_no_pairlist = c0 + h * (c1 + h * c2);
  } else {
    func_no_pairlist = (1.0 - xn) / (1.0 - xd);
  }

  cvm::real func, inv_one_pairlist_tol;
  if (flags & ef_use_pairlist) {
    inv_one_pairlist_tol = 1 / (1.0-pairlist_tol);
    func = (func_no_pairlist - pairlist_tol) * inv_one_pairlist_tol;
  } else {
    func = func_no_pairlist;
  }

  // If the value is too small and we are correcting for the tolerance, the result is negative
  // and we need to exclude it rather than let it contribute to the sum or the gradients.
  if (func < 0)
    return 0;

  if (flags & ef_gradients) {
    // Logarithmic derivative: 1st-order Taylor expansion around l2 = 1
    cvm::real log_deriv;
    if (std::abs(h) < eps_l2) {
      cvm::real const g0 = 0.5 * (en2_r - ed2_r);
      cvm::real const g1 = ((en2_r - ed2_r) * (en2_r + ed2_r - 6.0)) / 12.0;
      log_deriv = g0 + h * g1;
    } else {
      log_deriv = (ed2_r * xd / ((1.0 - xd) * l2)) - (en2_r * xn / ((1.0 - xn) * l2));
    }
    dFdl2 = (flags & ef_use_pairlist) ?
      func_no_pairlist * inv_one_pairlist_tol * log_deriv :
      func * log_deriv;
  }

  return func;
}


template<int flags>
inline cvm::real colvar::coordnum::compute_pair_coordnum(cvm::rvector const &inv_r0_vec,
                                                         cvm::rvector const &inv_r0sq_vec,
                                                         int en,
                                                         int ed,
                                                         const cvm::real a1x,
                                                         const cvm::real a1y,
                                                         const cvm::real a1z,
                                                         const cvm::real a2x,
                                                         const cvm::real a2y,
                                                         const cvm::real a2z,
                                                         cvm::real& g1x,
                                                         cvm::real& g1y,
                                                         cvm::real& g1z,
                                                         cvm::real& g2x,
                                                         cvm::real& g2y,
                                                         cvm::real& g2z,
                                                         cvm::real pairlist_tol,
                                                         cvm::real pairlist_tol_l2_max,
                                                         cvm::system_boundary_conditions const &bc)
{
  const cvm::atom_pos pos1{a1x, a1y, a1z};
  const cvm::atom_pos pos2{a2x, a2y, a2z};

  cvm::rvector const diff = bc.position_distance(pos1, pos2);
  cvm::rvector const scal_diff(diff.x * inv_r0_vec.x,
                               diff.y * inv_r0_vec.y,
                               diff.z * inv_r0_vec.z);
  cvm::real const l2 = scal_diff.norm2();
  if (flags & ef_use_pairlist) {
    if (l2 > pairlist_tol_l2_max) {
      // Exit if the distance is such that F(l2) < pairlist_tol
      return 0.0;
    }
  }

  cvm::real dFdl2 = 0.0;
  cvm::real F = switching_function<flags>(l2, dFdl2, en, ed, pairlist_tol);

  if ((flags & ef_gradients) && (F > 0.0)) {
    cvm::rvector const dl2dx((2.0 * inv_r0sq_vec.x) * diff.x,
                             (2.0 * inv_r0sq_vec.y) * diff.y,
                             (2.0 * inv_r0sq_vec.z) * diff.z);

    const cvm::rvector G = dFdl2*dl2dx;
    g1x += -1.0*G.x;
    g1y += -1.0*G.y;
    g1z += -1.0*G.z;
    g2x +=      G.x;
    g2y +=      G.y;
    g2z +=      G.z;
  }

  return F;
}

#endif // COLVARCOMP_COORDNUM_H
