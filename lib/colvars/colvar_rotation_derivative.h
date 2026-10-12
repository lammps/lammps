#ifndef COLVAR_ROTATION_DERIVATIVE
#define COLVAR_ROTATION_DERIVATIVE

#include "colvar_gpu_support.h"
#include "colvartypes.h"
#include <type_traits>
#include <cstring>
#include <array>

#ifndef _noalias
#if defined(__INTEL_COMPILER) || (defined(__PGI) && !defined(__NVCOMPILER))
#define _noalias restrict
#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER) || defined(__NVCOMPILER)
#define _noalias __restrict
#else
#define _noalias
#endif
#endif

/// \brief Helper enum class for specifying options in rotation_derivative::prepare_derivative
enum class rotation_derivative_dldq {
  /// Require the derivative of the leading eigenvalue with respect to the atom coordinates
  use_dl = 1 << 0,
  /// Require the derivative of the leading eigenvector with respect to the atom coordinates
  use_dq = 1 << 1
};

inline COLVARS_HOST_DEVICE constexpr rotation_derivative_dldq operator|(rotation_derivative_dldq Lhs, rotation_derivative_dldq Rhs) {
  return static_cast<rotation_derivative_dldq>(
    static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Lhs) |
    static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Rhs));
}

inline COLVARS_HOST_DEVICE constexpr bool operator&(rotation_derivative_dldq Lhs, rotation_derivative_dldq Rhs)
{
  return (static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Lhs) &
          static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Rhs));
}

/// \brief Helper class for calculating the derivative of rotation
// template <typename T1, typename T2, bool soa = false>
struct rotation_derivative {
  /// \brief Reference to the rotation
  const cvm::rotation &m_rot;
  /// \brief Reference to the atom positions of group 1
  // const std::vector<cvm::real> &m_pos1;
  cvm::ag_vector_real_t::const_iterator pos1x;
  cvm::ag_vector_real_t::const_iterator pos1y;
  cvm::ag_vector_real_t::const_iterator pos1z;
  /// \brief Reference to the atom positions of group 2
  // const std::vector<cvm::real> &m_pos2;
  cvm::ag_vector_real_t::const_iterator pos2x;
  cvm::ag_vector_real_t::const_iterator pos2y;
  cvm::ag_vector_real_t::const_iterator pos2z;
  /// \brief Number of atoms in group1 (used in SOA)
  size_t m_num_atoms_pos1;
  /// \brief Number of atoms in group1 (used in SOA)
  size_t m_num_atoms_pos2;
  /// \brief Temporary variable that will be updated if prepare_derivative called
  cvm::real tmp_Q0Q0[4][4];
  cvm::real tmp_Q0Q0_L[4][4][4];
  /*! @brief Constructor of the cvm::rotation::derivative class for SOA
    *  @param[in]  rot   The cvm::rotation object (must have called
    *                    `calc_optimal_rotation` before calling
    *                    `calc_derivative_wrt_group1` and
    *                    `calc_derivative_wrt_group2`)
    *  @param[in]  pos1  The atom positions of group 1
    *  @param[in]  pos2  The atom positions of group 2
    *  @param[in]  num_atoms_pos1 The number of atoms in group1
    *  @param[in]  num_atoms_pos2 The number of atoms in group2
    */
  rotation_derivative(
    const cvm::rotation &rot,
    const cvm::ag_vector_real_t &pos1,
    const cvm::ag_vector_real_t &pos2,
    const size_t num_atoms_pos1,
    const size_t num_atoms_pos2):
      m_rot(rot),
      pos1x(pos1.cbegin()),
      pos1y(pos1x + num_atoms_pos1),
      pos1z(pos1y + num_atoms_pos1),
      pos2x(pos2.cbegin()),
      pos2y(pos2x + num_atoms_pos2),
      pos2z(pos2y + num_atoms_pos2),
      m_num_atoms_pos1(num_atoms_pos1),
      m_num_atoms_pos2(num_atoms_pos2) {}
  /*! @brief This function must be called before `calc_derivative_wrt_group1`
    *         and `calc_derivative_wrt_group2` in order to prepare the tmp_Q0Q0
    *        and tmp_Q0Q0_L.
    *  @param[in] require_dl_dq Require the calculation of the derivatives of L or/and Q
    *                           with respect to atoms.
    */
  void prepare_derivative(rotation_derivative_dldq require_dl_dq) {
    if (require_dl_dq & rotation_derivative_dldq::use_dl) {
      const auto &Q0 = m_rot.S_eigvec[0];
      tmp_Q0Q0[0][0] = Q0[0] * Q0[0];
      tmp_Q0Q0[0][1] = Q0[0] * Q0[1];
      tmp_Q0Q0[0][2] = Q0[0] * Q0[2];
      tmp_Q0Q0[0][3] = Q0[0] * Q0[3];
      tmp_Q0Q0[1][0] = Q0[1] * Q0[0];
      tmp_Q0Q0[1][1] = Q0[1] * Q0[1];
      tmp_Q0Q0[1][2] = Q0[1] * Q0[2];
      tmp_Q0Q0[1][3] = Q0[1] * Q0[3];
      tmp_Q0Q0[2][0] = Q0[2] * Q0[0];
      tmp_Q0Q0[2][1] = Q0[2] * Q0[1];
      tmp_Q0Q0[2][2] = Q0[2] * Q0[2];
      tmp_Q0Q0[2][3] = Q0[2] * Q0[3];
      tmp_Q0Q0[3][0] = Q0[3] * Q0[0];
      tmp_Q0Q0[3][1] = Q0[3] * Q0[1];
      tmp_Q0Q0[3][2] = Q0[3] * Q0[2];
      tmp_Q0Q0[3][3] = Q0[3] * Q0[3];
    }
    if (require_dl_dq & rotation_derivative_dldq::use_dq) {
      const auto &Q0 = m_rot.S_eigvec[0];
      const auto &Q1 = m_rot.S_eigvec[1];
      const auto &Q2 = m_rot.S_eigvec[2];
      const auto &Q3 = m_rot.S_eigvec[3];
      cvm::real const L0 = m_rot.S_eigval[0];
      cvm::real const L1 = m_rot.S_eigval[1];
      cvm::real const L2 = m_rot.S_eigval[2];
      cvm::real const L3 = m_rot.S_eigval[3];

      tmp_Q0Q0_L[0][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[0];

      tmp_Q0Q0_L[0][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[0];

      tmp_Q0Q0_L[0][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[0];

      tmp_Q0Q0_L[0][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[0][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[0];

      tmp_Q0Q0_L[1][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[1];

      tmp_Q0Q0_L[1][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[1];

      tmp_Q0Q0_L[1][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[1];

      tmp_Q0Q0_L[1][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[1][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[1];

      tmp_Q0Q0_L[2][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[2];

      tmp_Q0Q0_L[2][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[2];

      tmp_Q0Q0_L[2][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[2];

      tmp_Q0Q0_L[2][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[2][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[2];

      tmp_Q0Q0_L[3][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[3][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[3][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[3][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[3];
      tmp_Q0Q0_L[3][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[3];
    }
  }
  /*! @brief Actual implementation of the derivative calculation
    *  @param[in]  ds  The derivative of matrix S with respect to an atom of
    *                  either group 1 or group 2
    *  @param[out] dl0_out The output of derivative of L
    *  @param[out] dq0_out The output of derivative of Q
    *  @param[out] ds_out  The output of derivative of overlap matrix S
    */
  template <bool use_dl, bool use_dq, bool use_ds>
  void calc_derivative_impl(
    const cvm::rvector (&ds)[4][4],
    cvm::rvector* _noalias const dl0_out,
    std::array<cvm::rvector, 4>* _noalias const dq0_out,
    std::array<std::array<cvm::rvector, 4>, 4>* _noalias const ds_out) const {
    if (use_ds) {
      // this code path is for debug_gradients, so not necessary to unroll the loop
      *ds_out = std::array<std::array<cvm::rvector, 4>, 4>();
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          (*ds_out)[i][j] = ds[i][j];
        }
      }
    }
    if (use_dl) {
      /* manually loop unrolling of the following loop:
        dl0_1.reset();
        for (size_t i = 0; i < 4; i++) {
          for (size_t j = 0; j < 4; j++) {
            dl0_1 += Q0[i] * ds_1[i][j] * Q0[j];
          }
        }
      */
      *dl0_out = tmp_Q0Q0[0][0] * ds[0][0] +
                 tmp_Q0Q0[0][1] * ds[0][1] +
                 tmp_Q0Q0[0][2] * ds[0][2] +
                 tmp_Q0Q0[0][3] * ds[0][3] +
                 tmp_Q0Q0[1][0] * ds[1][0] +
                 tmp_Q0Q0[1][1] * ds[1][1] +
                 tmp_Q0Q0[1][2] * ds[1][2] +
                 tmp_Q0Q0[1][3] * ds[1][3] +
                 tmp_Q0Q0[2][0] * ds[2][0] +
                 tmp_Q0Q0[2][1] * ds[2][1] +
                 tmp_Q0Q0[2][2] * ds[2][2] +
                 tmp_Q0Q0[2][3] * ds[2][3] +
                 tmp_Q0Q0[3][0] * ds[3][0] +
                 tmp_Q0Q0[3][1] * ds[3][1] +
                 tmp_Q0Q0[3][2] * ds[3][2] +
                 tmp_Q0Q0[3][3] * ds[3][3];
    }
    if (use_dq) {
      // we can skip this check if a fixed-size array is used
      // if (dq0_out->size() != 4) dq0_out->resize(4);
      /* manually loop unrolling of the following loop:
        dq0_1.reset();
        for (size_t p = 0; p < 4; p++) {
          for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
              dq0_1[p] +=
                (Q1[i] * ds_1[i][j] * Q0[j]) / (L0-L1) * Q1[p] +
                (Q2[i] * ds_1[i][j] * Q0[j]) / (L0-L2) * Q2[p] +
                (Q3[i] * ds_1[i][j] * Q0[j]) / (L0-L3) * Q3[p];
            }
          }
        }
      */
      (*dq0_out)[0] = tmp_Q0Q0_L[0][0][0] * ds[0][0] +
                      tmp_Q0Q0_L[0][0][1] * ds[0][1] +
                      tmp_Q0Q0_L[0][0][2] * ds[0][2] +
                      tmp_Q0Q0_L[0][0][3] * ds[0][3] +
                      tmp_Q0Q0_L[0][1][0] * ds[1][0] +
                      tmp_Q0Q0_L[0][1][1] * ds[1][1] +
                      tmp_Q0Q0_L[0][1][2] * ds[1][2] +
                      tmp_Q0Q0_L[0][1][3] * ds[1][3] +
                      tmp_Q0Q0_L[0][2][0] * ds[2][0] +
                      tmp_Q0Q0_L[0][2][1] * ds[2][1] +
                      tmp_Q0Q0_L[0][2][2] * ds[2][2] +
                      tmp_Q0Q0_L[0][2][3] * ds[2][3] +
                      tmp_Q0Q0_L[0][3][0] * ds[3][0] +
                      tmp_Q0Q0_L[0][3][1] * ds[3][1] +
                      tmp_Q0Q0_L[0][3][2] * ds[3][2] +
                      tmp_Q0Q0_L[0][3][3] * ds[3][3];

      (*dq0_out)[1] = tmp_Q0Q0_L[1][0][0] * ds[0][0] +
                      tmp_Q0Q0_L[1][0][1] * ds[0][1] +
                      tmp_Q0Q0_L[1][0][2] * ds[0][2] +
                      tmp_Q0Q0_L[1][0][3] * ds[0][3] +
                      tmp_Q0Q0_L[1][1][0] * ds[1][0] +
                      tmp_Q0Q0_L[1][1][1] * ds[1][1] +
                      tmp_Q0Q0_L[1][1][2] * ds[1][2] +
                      tmp_Q0Q0_L[1][1][3] * ds[1][3] +
                      tmp_Q0Q0_L[1][2][0] * ds[2][0] +
                      tmp_Q0Q0_L[1][2][1] * ds[2][1] +
                      tmp_Q0Q0_L[1][2][2] * ds[2][2] +
                      tmp_Q0Q0_L[1][2][3] * ds[2][3] +
                      tmp_Q0Q0_L[1][3][0] * ds[3][0] +
                      tmp_Q0Q0_L[1][3][1] * ds[3][1] +
                      tmp_Q0Q0_L[1][3][2] * ds[3][2] +
                      tmp_Q0Q0_L[1][3][3] * ds[3][3];

      (*dq0_out)[2] = tmp_Q0Q0_L[2][0][0] * ds[0][0] +
                      tmp_Q0Q0_L[2][0][1] * ds[0][1] +
                      tmp_Q0Q0_L[2][0][2] * ds[0][2] +
                      tmp_Q0Q0_L[2][0][3] * ds[0][3] +
                      tmp_Q0Q0_L[2][1][0] * ds[1][0] +
                      tmp_Q0Q0_L[2][1][1] * ds[1][1] +
                      tmp_Q0Q0_L[2][1][2] * ds[1][2] +
                      tmp_Q0Q0_L[2][1][3] * ds[1][3] +
                      tmp_Q0Q0_L[2][2][0] * ds[2][0] +
                      tmp_Q0Q0_L[2][2][1] * ds[2][1] +
                      tmp_Q0Q0_L[2][2][2] * ds[2][2] +
                      tmp_Q0Q0_L[2][2][3] * ds[2][3] +
                      tmp_Q0Q0_L[2][3][0] * ds[3][0] +
                      tmp_Q0Q0_L[2][3][1] * ds[3][1] +
                      tmp_Q0Q0_L[2][3][2] * ds[3][2] +
                      tmp_Q0Q0_L[2][3][3] * ds[3][3];

      (*dq0_out)[3] = tmp_Q0Q0_L[3][0][0] * ds[0][0] +
                      tmp_Q0Q0_L[3][0][1] * ds[0][1] +
                      tmp_Q0Q0_L[3][0][2] * ds[0][2] +
                      tmp_Q0Q0_L[3][0][3] * ds[0][3] +
                      tmp_Q0Q0_L[3][1][0] * ds[1][0] +
                      tmp_Q0Q0_L[3][1][1] * ds[1][1] +
                      tmp_Q0Q0_L[3][1][2] * ds[1][2] +
                      tmp_Q0Q0_L[3][1][3] * ds[1][3] +
                      tmp_Q0Q0_L[3][2][0] * ds[2][0] +
                      tmp_Q0Q0_L[3][2][1] * ds[2][1] +
                      tmp_Q0Q0_L[3][2][2] * ds[2][2] +
                      tmp_Q0Q0_L[3][2][3] * ds[2][3] +
                      tmp_Q0Q0_L[3][3][0] * ds[3][0] +
                      tmp_Q0Q0_L[3][3][1] * ds[3][1] +
                      tmp_Q0Q0_L[3][3][2] * ds[3][2] +
                      tmp_Q0Q0_L[3][3][3] * ds[3][3];
    }
  }
  /*! @brief Calculate the derivatives of S, the leading eigenvalue L and
    *         the leading eigenvector Q with respect to `m_pos1`
    *  @param[in]  ia        The index the of atom
    *  @param[out] dl0_1_out The output of derivative of L with respect to
    *                        ia-th atom of group 1
    *  @param[out] dq0_1_out The output of derivative of Q with respect to
    *                        ia-th atom of group 1
    *  @param[out] ds_1_out  The output of derivative of overlap matrix S with
    *                        respect to ia-th atom of group 1
    */
  template <bool use_dl, bool use_dq, bool use_ds>
  void calc_derivative_wrt_group1(
    size_t ia, cvm::rvector* _noalias const dl0_1_out = nullptr,
    std::array<cvm::rvector, 4>* _noalias const dq0_1_out = nullptr,
    std::array<std::array<cvm::rvector, 4>, 4>* _noalias const ds_1_out = nullptr) const {
      // if (dl0_1_out == nullptr && dq0_1_out == nullptr) return;
      const cvm::real a2x = *(pos2x + ia);
      const cvm::real a2y = *(pos2y + ia);
      const cvm::real a2z = *(pos2z + ia);
      const cvm::rvector ds_1[4][4] = {
        {{ a2x,  a2y,  a2z}, { 0.0, a2z,  -a2y}, {-a2z,  0.0,  a2x}, { a2y, -a2x,  0.0}},
        {{ 0.0,  a2z, -a2y}, { a2x, -a2y, -a2z}, { a2y,  a2x,  0.0}, { a2z,  0.0,  a2x}},
        {{-a2z,  0.0,  a2x}, { a2y,  a2x,  0.0}, {-a2x,  a2y, -a2z}, { 0.0,  a2z,  a2y}},
        {{ a2y, -a2x,  0.0}, { a2z,  0.0,  a2x}, { 0.0,  a2z,  a2y}, {-a2x, -a2y,  a2z}}};
      calc_derivative_impl<use_dl, use_dq, use_ds>(ds_1, dl0_1_out, dq0_1_out, ds_1_out);
    }
  /*! @brief Calculate the derivatives of S, the leading eigenvalue L and
    *         the leading eigenvector Q with respect to `m_pos2`
    *  @param[in]  ia        The index the of atom
    *  @param[out] dl0_2_out The output of derivative of L with respect to
    *                        ia-th atom of group 2
    *  @param[out] dq0_2_out The output of derivative of Q with respect to
    *                        ia-th atom of group 2
    *  @param[out] ds_2_out  The output of derivative of overlap matrix S with
    *                        respect to ia-th atom of group 2
    */
  template <bool use_dl, bool use_dq, bool use_ds>
  void calc_derivative_wrt_group2(
    size_t ia, cvm::rvector* _noalias const dl0_2_out = nullptr,
    std::array<cvm::rvector, 4>* _noalias const dq0_2_out = nullptr,
    std::array<std::array<cvm::rvector, 4>, 4>* _noalias const ds_2_out = nullptr) const {
    // if (dl0_2_out == nullptr && dq0_2_out == nullptr) return;
    const cvm::real a1x = *(pos1x + ia);
    const cvm::real a1y = *(pos1y + ia);
    const cvm::real a1z = *(pos1z + ia);
    const cvm::rvector ds_2[4][4] = {
      {{ a1x,  a1y,  a1z}, { 0.0, -a1z,  a1y}, { a1z,  0.0, -a1x}, {-a1y,  a1x,  0.0}},
      {{ 0.0, -a1z,  a1y}, { a1x, -a1y, -a1z}, { a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}},
      {{ a1z,  0.0, -a1x}, { a1y,  a1x,  0.0}, {-a1x,  a1y, -a1z}, { 0.0,  a1z,  a1y}},
      {{-a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}, { 0.0,  a1z,  a1y}, {-a1x, -a1y,  a1z}}};
    calc_derivative_impl<use_dl, use_dq, use_ds>(ds_2, dl0_2_out, dq0_2_out, ds_2_out);
  }

  /*! @brief Project the force on \f$q_i\f$ (or the gradient of
   *  \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$) to the force on the correlation
   *  matrix \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, where \f$\mathbf{r}_1\f$ and
   *  \f$\mathbf{r}_2\f$ are \f$N\times 3\f$ matrices containing the atom
   *  positions of atom group 1 and atom group 2 after centering on origin.
   *
   *  This function can be only called after
   *  prepare_derivative(rotation_derivative_dldq::use_dq). It mulitplies the input with
   *  \f$\frac{\mathrm{d}q_i}{\mathrm{d}C}\f$.
   *
   *  @tparam i The i-th component of the quaternion, must be in the range [0, 3]
   *  @param[in] f_on_q The force on \f$q_i\f$ or the gradient of \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$
   *
   *  @return A 3x3 matrix. The matrix element at row \f$j\f$ and column \f$k\f$
   *  is \f$f_{q_i}\frac{\mathrm{d}q_i}{\mathrm{d}C_{jk}}\f$.
   */
  template <int i>
  inline cvm::rmatrix project_force_to_C_from_dxdqi(cvm::real f_on_q) const {
    static_assert((i < 4) && (i >= 0), "i must be in [0, 3] in project_force_to_C_from_dxdqi.");
    cvm::rmatrix result;
    result.xx = f_on_q * ( tmp_Q0Q0_L[i][0][0] + tmp_Q0Q0_L[i][1][1] - tmp_Q0Q0_L[i][2][2] - tmp_Q0Q0_L[i][3][3] );
    result.xy = f_on_q * ( tmp_Q0Q0_L[i][0][3] + tmp_Q0Q0_L[i][1][2] + tmp_Q0Q0_L[i][2][1] + tmp_Q0Q0_L[i][3][0] );
    result.xz = f_on_q * (-tmp_Q0Q0_L[i][0][2] + tmp_Q0Q0_L[i][1][3] - tmp_Q0Q0_L[i][2][0] + tmp_Q0Q0_L[i][3][1] );
    result.yx = f_on_q * (-tmp_Q0Q0_L[i][0][3] + tmp_Q0Q0_L[i][1][2] + tmp_Q0Q0_L[i][2][1] - tmp_Q0Q0_L[i][3][0] );
    result.yy = f_on_q * ( tmp_Q0Q0_L[i][0][0] - tmp_Q0Q0_L[i][1][1] + tmp_Q0Q0_L[i][2][2] - tmp_Q0Q0_L[i][3][3] );
    result.yz = f_on_q * ( tmp_Q0Q0_L[i][0][1] + tmp_Q0Q0_L[i][1][0] + tmp_Q0Q0_L[i][2][3] + tmp_Q0Q0_L[i][3][2] );
    result.zx = f_on_q * ( tmp_Q0Q0_L[i][0][2] + tmp_Q0Q0_L[i][1][3] + tmp_Q0Q0_L[i][2][0] + tmp_Q0Q0_L[i][3][1] );
    result.zy = f_on_q * (-tmp_Q0Q0_L[i][0][1] - tmp_Q0Q0_L[i][1][0] + tmp_Q0Q0_L[i][2][3] + tmp_Q0Q0_L[i][3][2] );
    result.zz = f_on_q * ( tmp_Q0Q0_L[i][0][0] - tmp_Q0Q0_L[i][1][1] - tmp_Q0Q0_L[i][2][2] + tmp_Q0Q0_L[i][3][3] );
    return result;
  }

  /*! @brief Project the force on \f$\mathbf{q}\f$ (or the gradient of
   *  \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$) to the force on the correlation
   *  matrix \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, where \f$\mathbf{r}_1\f$ and
   *  \f$\mathbf{r}_2\f$ are \f$N\times 3\f$ matrices containing the atom
   *  positions of atom group 1 and atom group 2 after centering on origin.
   *
   *  See also project_force_to_C_from_dxdqi().
   *
   *  @tparam dim4_array_t The type of force acting on \f$\mathbf{q}\f$.
   *
   *  @param[in] sum_dxdq The force on \f$\mathbf{q}\f$ or the gradient vector
   *  \f$(\frac{\mathrm{d}x}{\mathrm{d}q_0}, \frac{\mathrm{d}x}{\mathrm{d}q_1}, \frac{\mathrm{d}x}{\mathrm{d}q_2}, \frac{\mathrm{d}x}{\mathrm{d}q_3})\f$.
   *
   *  @return A 3x3 matrix. The matrix element at row \f$j\f$ and column \f$k\f$
   *  is \f$\sum_{i=0}^3 f_{q_i}\frac{\mathrm{d}q_i}{\mathrm{d}C_{jk}}\f$.
   *
   */
  template <typename dim4_array_t>
  inline cvm::rmatrix project_force_to_C_from_dxdq(const dim4_array_t& sum_dxdq) const {
    cvm::rmatrix result;
    result += project_force_to_C_from_dxdqi<0>(sum_dxdq[0]);
    result += project_force_to_C_from_dxdqi<1>(sum_dxdq[1]);
    result += project_force_to_C_from_dxdqi<2>(sum_dxdq[2]);
    result += project_force_to_C_from_dxdqi<3>(sum_dxdq[3]);
    return result;
  }

  /*! @brief Project the force on the correlation matrix \f$C\f$ to \f$\mathbf{r}_1\f$
   *
   *  Let \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, and the force on \f$C\f$ be
   *  \f$F_{C}\f$. This function returns the force on the i-th atom of \f$\mathbf{r}_1\f$.
   *
   *  @param[in] ia The atom index of the i-th atom in \f$\mathbf{r}_1\f$
   *  @param[in] dxdC The 3x3 matrix containing the forces on each element of \f$C\f$
   *
   *  @return The force on the i-th atom of \f$\mathbf{r}_1\f$
   */
  inline cvm::rvector project_force_to_group1(size_t ia, const cvm::rmatrix& dxdC) const {
    const cvm::real a2x = *(pos2x + ia);
    const cvm::real a2y = *(pos2y + ia);
    const cvm::real a2z = *(pos2z + ia);
    const cvm::rvector result{
      dxdC.xx * a2x + dxdC.xy * a2y + dxdC.xz * a2z,
      dxdC.yx * a2x + dxdC.yy * a2y + dxdC.yz * a2z,
      dxdC.zx * a2x + dxdC.zy * a2y + dxdC.zz * a2z};
    return result;
  }

  /*! @brief Project the force on the correlation matrix \f$C\f$ to \f$\mathbf{r}_2\f$
   *
   *  Let \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, and the force on \f$C\f$ be
   *  \f$F_{C}\f$. This function returns the force on the i-th atom of \f$\mathbf{r}_2\f$.
   *
   *  @param[in] ia The atom index of the i-th atom in \f$\mathbf{r}_2\f$
   *  @param[in] dxdC The 3x3 matrix containing the forces on each element of \f$C\f$
   *
   *  @return The force on the i-th atom of \f$\mathbf{r}_2\f$
   */
  inline cvm::rvector project_force_to_group2(size_t ia, const cvm::rmatrix& dxdC) const {
    const cvm::real a1x = *(pos1x + ia);
    const cvm::real a1y = *(pos1y + ia);
    const cvm::real a1z = *(pos1z + ia);
    const cvm::rvector result{
      dxdC.xx * a1x + dxdC.yx * a1y + dxdC.zx * a1z,
      dxdC.xy * a1x + dxdC.yy * a1y + dxdC.zy * a1z,
      dxdC.xz * a1x + dxdC.yz * a1y + dxdC.zz * a1z};
    return result;
  }
};

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
namespace colvars_gpu {
struct rotation_derivative_gpu {
  /// \brief Pointer to the parent colvarmodule
  colvarmodule* cvmodule;
  /// \brief Reference to the rotation
  const colvars_gpu::rotation_gpu *m_rot;
  /// \brief Reference to the atom positions of group 1
  const cvm::real* d_pos1x;
  const cvm::real* d_pos1y;
  const cvm::real* d_pos1z;
  /// \brief Reference to the atom positions of group 2
  const cvm::real* d_pos2x;
  const cvm::real* d_pos2y;
  const cvm::real* d_pos2z;
  /// \brief Number of atoms in group1 (used in SOA)
  size_t m_num_atoms_pos1;
  /// \brief Number of atoms in group1 (used in SOA)
  size_t m_num_atoms_pos2;
  /// \brief GPU temporary variable that will be updated if prepare_derivative called
  cvm::real* tmp_Q0Q0;
  cvm::real* tmp_Q0Q0_L;
  /*! @brief Constructor
   *
   *  The object of this class is expected to be constructed on host-pinned memory.
   */
  rotation_derivative_gpu(colvarmodule* cmodule_in);
  /*! @brief Initialization of the rotation_derivative_gpu class
   *  @param[in]  rot   The colvars_gpu::rotation_gpu object
   *  @param[in]  pos1  The atom positions of group 1
   *  @param[in]  pos2  The atom positions of group 2
   *  @param[in]  num_atoms_pos1 The number of atoms in group1
   *  @param[in]  num_atoms_pos2 The number of atoms in group2
   */
  int init(
    const colvars_gpu::rotation_gpu* rot,
    const cvm::real* d_pos1,
    const cvm::real* d_pos2,
    const size_t num_atoms_pos1,
    const size_t num_atoms_pos2);
  ~rotation_derivative_gpu();
  /*! @brief Add a "prepare_rotation_derivative" node to the CUDA graph
   *  @param[in] require_dl_dq Require the calculation of the derivatives of L or/and Q
   *                           with respect to atoms.
   *  @param[in] graph The CUDA graph where the node would be added
   *  @param[out] nodes_map A map from operation names to cudaGraphNode_t
   *  for dependency lookup. A node termed as "prepare_rotation_derivative" would
   *  be added into the nodes_map.
   *  @return COLVARS_OK if success and COLVARS_ERROR otherwise
   *
   * The corresponding graph must be executed before `calc_derivative_wrt_group1`
   * and `calc_derivative_wrt_group2` in order to prepare the tmp_Q0Q0
   * and tmp_Q0Q0_L.
   *
   */
  int add_prepare_derivative_nodes(
    rotation_derivative_dldq require_dl_dq,
    cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);

  /*! @brief Actual implementation of the derivative calculation
   *  @param[in]  ds  The derivative of matrix S with respect to an atom of
   *                  either group 1 or group 2
   *  @param[out] dl0_out The output of derivative of L
   *  @param[out] dq0_out The output of derivative of Q
   */
  template <bool use_dl, bool use_dq/*, bool use_ds*/>
  inline COLVARS_DEVICE void calc_derivative_impl(
    const cvm::rvector (&ds)[4][4],
    cvm::rvector* _noalias const dl0_out,
    cvm::rvector* _noalias const dq0_out) const {
    if (use_dl) {
      *dl0_out = tmp_Q0Q0[0*4+0] * ds[0][0] +
                 tmp_Q0Q0[0*4+1] * ds[0][1] +
                 tmp_Q0Q0[0*4+2] * ds[0][2] +
                 tmp_Q0Q0[0*4+3] * ds[0][3] +
                 tmp_Q0Q0[1*4+0] * ds[1][0] +
                 tmp_Q0Q0[1*4+1] * ds[1][1] +
                 tmp_Q0Q0[1*4+2] * ds[1][2] +
                 tmp_Q0Q0[1*4+3] * ds[1][3] +
                 tmp_Q0Q0[2*4+0] * ds[2][0] +
                 tmp_Q0Q0[2*4+1] * ds[2][1] +
                 tmp_Q0Q0[2*4+2] * ds[2][2] +
                 tmp_Q0Q0[2*4+3] * ds[2][3] +
                 tmp_Q0Q0[3*4+0] * ds[3][0] +
                 tmp_Q0Q0[3*4+1] * ds[3][1] +
                 tmp_Q0Q0[3*4+2] * ds[3][2] +
                 tmp_Q0Q0[3*4+3] * ds[3][3];
    }
    if (use_dq) {
      dq0_out[0] = tmp_Q0Q0_L[0*16+0*4+0] * ds[0][0] +
                   tmp_Q0Q0_L[0*16+0*4+1] * ds[0][1] +
                   tmp_Q0Q0_L[0*16+0*4+2] * ds[0][2] +
                   tmp_Q0Q0_L[0*16+0*4+3] * ds[0][3] +
                   tmp_Q0Q0_L[0*16+1*4+0] * ds[1][0] +
                   tmp_Q0Q0_L[0*16+1*4+1] * ds[1][1] +
                   tmp_Q0Q0_L[0*16+1*4+2] * ds[1][2] +
                   tmp_Q0Q0_L[0*16+1*4+3] * ds[1][3] +
                   tmp_Q0Q0_L[0*16+2*4+0] * ds[2][0] +
                   tmp_Q0Q0_L[0*16+2*4+1] * ds[2][1] +
                   tmp_Q0Q0_L[0*16+2*4+2] * ds[2][2] +
                   tmp_Q0Q0_L[0*16+2*4+3] * ds[2][3] +
                   tmp_Q0Q0_L[0*16+3*4+0] * ds[3][0] +
                   tmp_Q0Q0_L[0*16+3*4+1] * ds[3][1] +
                   tmp_Q0Q0_L[0*16+3*4+2] * ds[3][2] +
                   tmp_Q0Q0_L[0*16+3*4+3] * ds[3][3];

      dq0_out[1] = tmp_Q0Q0_L[1*16+0*4+0] * ds[0][0] +
                   tmp_Q0Q0_L[1*16+0*4+1] * ds[0][1] +
                   tmp_Q0Q0_L[1*16+0*4+2] * ds[0][2] +
                   tmp_Q0Q0_L[1*16+0*4+3] * ds[0][3] +
                   tmp_Q0Q0_L[1*16+1*4+0] * ds[1][0] +
                   tmp_Q0Q0_L[1*16+1*4+1] * ds[1][1] +
                   tmp_Q0Q0_L[1*16+1*4+2] * ds[1][2] +
                   tmp_Q0Q0_L[1*16+1*4+3] * ds[1][3] +
                   tmp_Q0Q0_L[1*16+2*4+0] * ds[2][0] +
                   tmp_Q0Q0_L[1*16+2*4+1] * ds[2][1] +
                   tmp_Q0Q0_L[1*16+2*4+2] * ds[2][2] +
                   tmp_Q0Q0_L[1*16+2*4+3] * ds[2][3] +
                   tmp_Q0Q0_L[1*16+3*4+0] * ds[3][0] +
                   tmp_Q0Q0_L[1*16+3*4+1] * ds[3][1] +
                   tmp_Q0Q0_L[1*16+3*4+2] * ds[3][2] +
                   tmp_Q0Q0_L[1*16+3*4+3] * ds[3][3];

      dq0_out[2] = tmp_Q0Q0_L[2*16+0*4+0] * ds[0][0] +
                   tmp_Q0Q0_L[2*16+0*4+1] * ds[0][1] +
                   tmp_Q0Q0_L[2*16+0*4+2] * ds[0][2] +
                   tmp_Q0Q0_L[2*16+0*4+3] * ds[0][3] +
                   tmp_Q0Q0_L[2*16+1*4+0] * ds[1][0] +
                   tmp_Q0Q0_L[2*16+1*4+1] * ds[1][1] +
                   tmp_Q0Q0_L[2*16+1*4+2] * ds[1][2] +
                   tmp_Q0Q0_L[2*16+1*4+3] * ds[1][3] +
                   tmp_Q0Q0_L[2*16+2*4+0] * ds[2][0] +
                   tmp_Q0Q0_L[2*16+2*4+1] * ds[2][1] +
                   tmp_Q0Q0_L[2*16+2*4+2] * ds[2][2] +
                   tmp_Q0Q0_L[2*16+2*4+3] * ds[2][3] +
                   tmp_Q0Q0_L[2*16+3*4+0] * ds[3][0] +
                   tmp_Q0Q0_L[2*16+3*4+1] * ds[3][1] +
                   tmp_Q0Q0_L[2*16+3*4+2] * ds[3][2] +
                   tmp_Q0Q0_L[2*16+3*4+3] * ds[3][3];

      dq0_out[3] = tmp_Q0Q0_L[3*16+0*4+0] * ds[0][0] +
                   tmp_Q0Q0_L[3*16+0*4+1] * ds[0][1] +
                   tmp_Q0Q0_L[3*16+0*4+2] * ds[0][2] +
                   tmp_Q0Q0_L[3*16+0*4+3] * ds[0][3] +
                   tmp_Q0Q0_L[3*16+1*4+0] * ds[1][0] +
                   tmp_Q0Q0_L[3*16+1*4+1] * ds[1][1] +
                   tmp_Q0Q0_L[3*16+1*4+2] * ds[1][2] +
                   tmp_Q0Q0_L[3*16+1*4+3] * ds[1][3] +
                   tmp_Q0Q0_L[3*16+2*4+0] * ds[2][0] +
                   tmp_Q0Q0_L[3*16+2*4+1] * ds[2][1] +
                   tmp_Q0Q0_L[3*16+2*4+2] * ds[2][2] +
                   tmp_Q0Q0_L[3*16+2*4+3] * ds[2][3] +
                   tmp_Q0Q0_L[3*16+3*4+0] * ds[3][0] +
                   tmp_Q0Q0_L[3*16+3*4+1] * ds[3][1] +
                   tmp_Q0Q0_L[3*16+3*4+2] * ds[3][2] +
                   tmp_Q0Q0_L[3*16+3*4+3] * ds[3][3];
    }
  }
  /*! @brief Calculate the derivatives of S, the leading eigenvalue L and
   *         the leading eigenvector Q with respect to `m_pos1`
   *  @param[in]  ia        The index the of atom
   *  @param[out] dl0_1_out The output of derivative of L with respect to
   *                        ia-th atom of group 1
   *  @param[out] dq0_1_out The output of derivative of Q with respect to
   *                        ia-th atom of group 1
   */
  template <bool use_dl, bool use_dq/*, bool use_ds*/>
  inline COLVARS_DEVICE void calc_derivative_wrt_group1(
    int ia,
    cvm::rvector* _noalias const dl0_1_out,
    cvm::rvector* _noalias const dq0_1_out) const {
    const cvm::real a2x = d_pos2x[ia];
    const cvm::real a2y = d_pos2y[ia];
    const cvm::real a2z = d_pos2z[ia];
    const cvm::rvector ds_1[4][4] = {
      {{ a2x,  a2y,  a2z}, { 0.0, a2z,  -a2y}, {-a2z,  0.0,  a2x}, { a2y, -a2x,  0.0}},
      {{ 0.0,  a2z, -a2y}, { a2x, -a2y, -a2z}, { a2y,  a2x,  0.0}, { a2z,  0.0,  a2x}},
      {{-a2z,  0.0,  a2x}, { a2y,  a2x,  0.0}, {-a2x,  a2y, -a2z}, { 0.0,  a2z,  a2y}},
      {{ a2y, -a2x,  0.0}, { a2z,  0.0,  a2x}, { 0.0,  a2z,  a2y}, {-a2x, -a2y,  a2z}}};
    calc_derivative_impl<use_dl, use_dq>(ds_1, dl0_1_out, dq0_1_out);
  }
  /*! @brief Calculate the derivatives of S, the leading eigenvalue L and
   *         the leading eigenvector Q with respect to `m_pos2`
   *  @param[in]  ia        The index the of atom
   *  @param[out] dl0_2_out The output of derivative of L with respect to
   *                        ia-th atom of group 2
   *  @param[out] dq0_2_out The output of derivative of Q with respect to
   *                        ia-th atom of group 2
   */
  template <bool use_dl, bool use_dq/*, bool use_ds*/>
  inline COLVARS_DEVICE void calc_derivative_wrt_group2(
    int ia,
    cvm::rvector* _noalias const dl0_2_out,
    cvm::rvector* _noalias const dq0_2_out) const {
    const cvm::real a1x = d_pos1x[ia];
    const cvm::real a1y = d_pos1y[ia];
    const cvm::real a1z = d_pos1z[ia];
    const cvm::rvector ds_2[4][4] = {
      {{ a1x,  a1y,  a1z}, { 0.0, -a1z,  a1y}, { a1z,  0.0, -a1x}, {-a1y,  a1x,  0.0}},
      {{ 0.0, -a1z,  a1y}, { a1x, -a1y, -a1z}, { a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}},
      {{ a1z,  0.0, -a1x}, { a1y,  a1x,  0.0}, {-a1x,  a1y, -a1z}, { 0.0,  a1z,  a1y}},
      {{-a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}, { 0.0,  a1z,  a1y}, {-a1x, -a1y,  a1z}}};
    calc_derivative_impl<use_dl, use_dq>(ds_2, dl0_2_out, dq0_2_out);
  }

  /*! @brief Project the force on \f$q_i\f$ (or the gradient of
   *  \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$) to the force on the correlation
   *  matrix \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, where \f$\mathbf{r}_1\f$ and
   *  \f$\mathbf{r}_2\f$ are \f$N\times 3\f$ matrices containing the atom
   *  positions of atom group 1 and atom group 2 after centering on origin.
   *
   *  This function can be only called after
   *  prepare_derivative(rotation_derivative_dldq::use_dq). It mulitplies the input with
   *  \f$\frac{\mathrm{d}q_i}{\mathrm{d}C}\f$.
   *
   *  @tparam i The i-th component of the quaternion, must be in the range [0, 3]
   *  @param[in] f_on_q The force on \f$q_i\f$ or the gradient of \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$
   *
   *  @return A 3x3 matrix. The matrix element at row \f$j\f$ and column \f$k\f$
   *  is \f$f_{q_i}\frac{\mathrm{d}q_i}{\mathrm{d}C_{jk}}\f$.
   */
  template <int i>
  inline COLVARS_DEVICE cvm::rmatrix project_force_to_C_from_dxdqi(cvm::real f_on_q) const {
    static_assert((i < 4) && (i >= 0), "i must be in [0, 3] in project_force_to_C_from_dxdqi.");
    return project_force_to_C_from_dxdqi(i, f_on_q);
  }

  inline COLVARS_DEVICE cvm::rmatrix project_force_to_C_from_dxdqi(int i, cvm::real f_on_q) const {
    cvm::rmatrix result;
    result.xx = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+0] + tmp_Q0Q0_L[i*16+1*4+1] - tmp_Q0Q0_L[i*16+2*4+2] - tmp_Q0Q0_L[i*16+3*4+3] );
    result.xy = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+3] + tmp_Q0Q0_L[i*16+1*4+2] + tmp_Q0Q0_L[i*16+2*4+1] + tmp_Q0Q0_L[i*16+3*4+0] );
    result.xz = f_on_q * (-tmp_Q0Q0_L[i*16+0*4+2] + tmp_Q0Q0_L[i*16+1*4+3] - tmp_Q0Q0_L[i*16+2*4+0] + tmp_Q0Q0_L[i*16+3*4+1] );
    result.yx = f_on_q * (-tmp_Q0Q0_L[i*16+0*4+3] + tmp_Q0Q0_L[i*16+1*4+2] + tmp_Q0Q0_L[i*16+2*4+1] - tmp_Q0Q0_L[i*16+3*4+0] );
    result.yy = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+0] - tmp_Q0Q0_L[i*16+1*4+1] + tmp_Q0Q0_L[i*16+2*4+2] - tmp_Q0Q0_L[i*16+3*4+3] );
    result.yz = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+1] + tmp_Q0Q0_L[i*16+1*4+0] + tmp_Q0Q0_L[i*16+2*4+3] + tmp_Q0Q0_L[i*16+3*4+2] );
    result.zx = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+2] + tmp_Q0Q0_L[i*16+1*4+3] + tmp_Q0Q0_L[i*16+2*4+0] + tmp_Q0Q0_L[i*16+3*4+1] );
    result.zy = f_on_q * (-tmp_Q0Q0_L[i*16+0*4+1] - tmp_Q0Q0_L[i*16+1*4+0] + tmp_Q0Q0_L[i*16+2*4+3] + tmp_Q0Q0_L[i*16+3*4+2] );
    result.zz = f_on_q * ( tmp_Q0Q0_L[i*16+0*4+0] - tmp_Q0Q0_L[i*16+1*4+1] - tmp_Q0Q0_L[i*16+2*4+2] + tmp_Q0Q0_L[i*16+3*4+3] );
    return result;
  }

  /*! @brief Project the force on \f$\mathbf{q}\f$ (or the gradient of
   *  \f$\frac{\mathrm{d}x}{\mathrm{d}q_i}\f$) to the force on the correlation
   *  matrix \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, where \f$\mathbf{r}_1\f$ and
   *  \f$\mathbf{r}_2\f$ are \f$N\times 3\f$ matrices containing the atom
   *  positions of atom group 1 and atom group 2 after centering on origin.
   *
   *  See also project_force_to_C_from_dxdqi().
   *
   *  @tparam dim4_array_t The type of force acting on \f$\mathbf{q}\f$.
   *
   *  @param[in] sum_dxdq The force on \f$\mathbf{q}\f$ or the gradient vector
   *  \f$(\frac{\mathrm{d}x}{\mathrm{d}q_0}, \frac{\mathrm{d}x}{\mathrm{d}q_1}, \frac{\mathrm{d}x}{\mathrm{d}q_2}, \frac{\mathrm{d}x}{\mathrm{d}q_3})\f$.
   *
   *  @return A 3x3 matrix. The matrix element at row \f$j\f$ and column \f$k\f$
   *  is \f$\sum_{i=0}^3 f_{q_i}\frac{\mathrm{d}q_i}{\mathrm{d}C_{jk}}\f$.
   *
   */
  template <typename dim4_array_t>
  inline COLVARS_DEVICE cvm::rmatrix project_force_to_C_from_dxdq(const dim4_array_t& sum_dxdq) const {
    cvm::rmatrix result;
    result += project_force_to_C_from_dxdqi<0>(sum_dxdq[0]);
    result += project_force_to_C_from_dxdqi<1>(sum_dxdq[1]);
    result += project_force_to_C_from_dxdqi<2>(sum_dxdq[2]);
    result += project_force_to_C_from_dxdqi<3>(sum_dxdq[3]);
    return result;
  }

  /*! @brief Project the force on the correlation matrix \f$C\f$ to \f$\mathbf{r}_1\f$
   *
   *  Let \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, and the force on \f$C\f$ be
   *  \f$F_{C}\f$. This function returns the force on the i-th atom of \f$\mathbf{r}_1\f$.
   *
   *  @param[in] ia The atom index of the i-th atom in \f$\mathbf{r}_1\f$
   *  @param[in] dxdC The 3x3 matrix containing the forces on each element of \f$C\f$
   *
   *  @return The force on the i-th atom of \f$\mathbf{r}_1\f$
   */
  inline COLVARS_DEVICE cvm::rvector project_force_to_group1(size_t ia, const cvm::rmatrix& dxdC) const {
    const cvm::real a2x = *(d_pos2x + ia);
    const cvm::real a2y = *(d_pos2y + ia);
    const cvm::real a2z = *(d_pos2z + ia);
    const cvm::rvector result{
      dxdC.xx * a2x + dxdC.xy * a2y + dxdC.xz * a2z,
      dxdC.yx * a2x + dxdC.yy * a2y + dxdC.yz * a2z,
      dxdC.zx * a2x + dxdC.zy * a2y + dxdC.zz * a2z};
    return result;
  }

  /*! @brief Project the force on the correlation matrix \f$C\f$ to \f$\mathbf{r}_2\f$
   *
   *  Let \f$C=\mathbf{r}_1^\intercal\mathbf{r}_2\f$, and the force on \f$C\f$ be
   *  \f$F_{C}\f$. This function returns the force on the i-th atom of \f$\mathbf{r}_2\f$.
   *
   *  @param[in] ia The atom index of the i-th atom in \f$\mathbf{r}_2\f$
   *  @param[in] dxdC The 3x3 matrix containing the forces on each element of \f$C\f$
   *
   *  @return The force on the i-th atom of \f$\mathbf{r}_2\f$
   */
  inline COLVARS_DEVICE cvm::rvector project_force_to_group2(size_t ia, const cvm::rmatrix& dxdC) const {
    const cvm::real a1x = *(d_pos1x + ia);
    const cvm::real a1y = *(d_pos1y + ia);
    const cvm::real a1z = *(d_pos1z + ia);
    const cvm::rvector result{
      dxdC.xx * a1x + dxdC.yx * a1y + dxdC.zx * a1z,
      dxdC.xy * a1x + dxdC.yy * a1y + dxdC.zy * a1z,
      dxdC.xz * a1x + dxdC.yz * a1y + dxdC.zz * a1z};
    return result;
  }
};
}
#elif defined(COLVARS_SYCL)
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)

#endif // COLVAR_ROTATION_DERIVATIVE
