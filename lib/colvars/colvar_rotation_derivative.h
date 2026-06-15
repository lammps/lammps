#ifndef COLVAR_ROTATION_DERIVATIVE
#define COLVAR_ROTATION_DERIVATIVE

#include "colvartypes.h"
#include <type_traits>
#include <cstring>

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

inline constexpr rotation_derivative_dldq operator|(rotation_derivative_dldq Lhs, rotation_derivative_dldq Rhs) {
  return static_cast<rotation_derivative_dldq>(
    static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Lhs) |
    static_cast<std::underlying_type<rotation_derivative_dldq>::type>(Rhs));
}

inline constexpr bool operator&(rotation_derivative_dldq Lhs, rotation_derivative_dldq Rhs)
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
  const std::vector<cvm::real> &m_pos1;
  /// \brief Reference to the atom positions of group 2
  const std::vector<cvm::real> &m_pos2;
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
    const std::vector<cvm::real> &pos1,
    const std::vector<cvm::real> &pos2,
    const size_t num_atoms_pos1,
    const size_t num_atoms_pos2):
      m_rot(rot), m_pos1(pos1), m_pos2(pos2),
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
      tmp_Q0Q0_L[1][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][0][0] = (Q1[0] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[0]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][0][1] = (Q1[0] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[1]) / (L0-L3) * Q3[3];


      tmp_Q0Q0_L[0][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][0][2] = (Q1[0] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[2]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][0][3] = (Q1[0] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[0] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[0] * Q0[3]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][1][0] = (Q1[1] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[0]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][1][1] = (Q1[1] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[1]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][1][2] = (Q1[1] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[2]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][1][3] = (Q1[1] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[1] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[1] * Q0[3]) / (L0-L3) * Q3[3];


      tmp_Q0Q0_L[0][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][2][0] = (Q1[2] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[0]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][2][1] = (Q1[2] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[1]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][2][2] = (Q1[2] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[2]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][2][3] = (Q1[2] * Q0[3]) / (L0-L1) * Q1[3] +
                            (Q2[2] * Q0[3]) / (L0-L2) * Q2[3] +
                            (Q3[2] * Q0[3]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][3][0] = (Q1[3] * Q0[0]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[0]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[0]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][3][1] = (Q1[3] * Q0[1]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[1]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[1]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[2];
      tmp_Q0Q0_L[3][3][2] = (Q1[3] * Q0[2]) / (L0-L1) * Q1[3] +
                            (Q2[3] * Q0[2]) / (L0-L2) * Q2[3] +
                            (Q3[3] * Q0[2]) / (L0-L3) * Q3[3];

      tmp_Q0Q0_L[0][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[0] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[0] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[0];
      tmp_Q0Q0_L[1][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[1] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[1] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[1];
      tmp_Q0Q0_L[2][3][3] = (Q1[3] * Q0[3]) / (L0-L1) * Q1[2] +
                            (Q2[3] * Q0[3]) / (L0-L2) * Q2[2] +
                            (Q3[3] * Q0[3]) / (L0-L3) * Q3[2];
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
    cvm::vector1d<cvm::rvector>* _noalias const dq0_out,
    cvm::matrix2d<cvm::rvector>* _noalias const ds_out) const {
    if (use_ds) {
      // this code path is for debug_gradients, so not necessary to unroll the loop
      *ds_out = cvm::matrix2d<cvm::rvector>(4, 4);
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
      if (dq0_out->size() != 4) dq0_out->resize(4);
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
    cvm::vector1d<cvm::rvector>* _noalias const dq0_1_out = nullptr,
    cvm::matrix2d<cvm::rvector>* _noalias const ds_1_out = nullptr) const {
      // if (dl0_1_out == nullptr && dq0_1_out == nullptr) return;
      const cvm::real a2x = m_pos2[ia];
      const cvm::real a2y = m_pos2[ia + m_num_atoms_pos2];
      const cvm::real a2z = m_pos2[ia + 2 * m_num_atoms_pos2];
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
    cvm::vector1d<cvm::rvector>* _noalias const dq0_2_out = nullptr,
    cvm::matrix2d<cvm::rvector>* _noalias const ds_2_out = nullptr) const {
    // if (dl0_2_out == nullptr && dq0_2_out == nullptr) return;
    const cvm::real a1x = m_pos1[ia];
    const cvm::real a1y = m_pos1[ia + m_num_atoms_pos1];
    const cvm::real a1z = m_pos1[ia + 2 * m_num_atoms_pos1];
    const cvm::rvector ds_2[4][4] = {
      {{ a1x,  a1y,  a1z}, { 0.0, -a1z,  a1y}, { a1z,  0.0, -a1x}, {-a1y,  a1x,  0.0}},
      {{ 0.0, -a1z,  a1y}, { a1x, -a1y, -a1z}, { a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}},
      {{ a1z,  0.0, -a1x}, { a1y,  a1x,  0.0}, {-a1x,  a1y, -a1z}, { 0.0,  a1z,  a1y}},
      {{-a1y,  a1x,  0.0}, { a1z,  0.0,  a1x}, { 0.0,  a1z,  a1y}, {-a1x, -a1y,  a1z}}};
    calc_derivative_impl<use_dl, use_dq, use_ds>(ds_2, dl0_2_out, dq0_2_out, ds_2_out);
  }
};

#endif // COLVAR_ROTATION_DERIVATIVE
