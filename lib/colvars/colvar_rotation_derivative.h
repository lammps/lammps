#ifndef COLVAR_ROTATION_DERIVATIVE
#define COLVAR_ROTATION_DERIVATIVE

#include "colvartypes.h"
#include <type_traits>
#include <cstring>

/// \brief Helper function for loading the ia-th atom in the vector pos to x, y and z (C++11 SFINAE is used)
template <typename T, typename std::enable_if<std::is_same<T, cvm::atom_pos>::value, bool>::type = true>
inline void read_atom_coord(
  size_t ia, const std::vector<T>& pos,
  cvm::real* x, cvm::real* y, cvm::real* z) {
  *x = pos[ia].x;
  *y = pos[ia].y;
  *z = pos[ia].z;
}

template <typename T, typename std::enable_if<std::is_same<T, cvm::atom>::value, bool>::type = true>
inline void read_atom_coord(
  size_t ia, const std::vector<T>& pos,
  cvm::real* x, cvm::real* y, cvm::real* z) {
  *x = pos[ia].pos.x;
  *y = pos[ia].pos.y;
  *z = pos[ia].pos.z;
}

/// \brief Helper enum class for specifying options in rotation_derivative::prepare_derivative
enum class rotation_derivative_dldq {
  /// Require the derivative of the leading eigenvalue with respect to the atom coordinats
  use_dl = 1 << 0,
  /// Require the derivative of the leading eigenvector with respect to the atom coordinats
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
template <typename T1, typename T2>
struct rotation_derivative {
  static_assert(std::is_same<T1, cvm::atom_pos>::value || std::is_same<T1, cvm::atom>::value,
                "class template rotation_derivative only supports cvm::atom_pos or cvm::atom types.");
  static_assert(std::is_same<T2, cvm::atom_pos>::value || std::is_same<T2, cvm::atom>::value,
                "class template rotation_derivative only supports cvm::atom_pos or cvm::atom types.");
  /// \brief Reference to the rotation
  const cvm::rotation &m_rot;
  /// \brief Reference to the atom positions of group 1
  const std::vector<T1> &m_pos1;
  /// \brief Reference to the atom positions of group 2
  const std::vector<T2> &m_pos2;
  /// \brief Temporary variable that will be updated if prepare_derivative called
  cvm::real tmp_Q0Q0[4][4];
  cvm::real tmp_Q0Q0_L[4][4][4];
  /*! @brief Constructor of the cvm::rotation::derivative class
    *  @param[in]  rot   The cvm::rotation object (must have called
    *                    `calc_optimal_rotation` before calling
    *                    `calc_derivative_wrt_group1` and
    *                    `calc_derivative_wrt_group2`)
    *  @param[in]  pos1  The atom positions of group 1
    *  @param[in]  pos2  The atom positions of group 2
    */
  rotation_derivative(
    const cvm::rotation &rot,
    const std::vector<T1> &pos1,
    const std::vector<T2> &pos2):
      m_rot(rot), m_pos1(pos1), m_pos2(pos2) {};
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
  void calc_derivative_impl(
    const cvm::rvector (&ds)[4][4],
    cvm::rvector* const dl0_out,
    cvm::vector1d<cvm::rvector>* const dq0_out,
    cvm::matrix2d<cvm::rvector>* const ds_out) const {
    if (ds_out != nullptr) {
      // this code path is for debug_gradients, so not necessary to unroll the loop
      *ds_out = cvm::matrix2d<cvm::rvector>(4, 4);
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          (*ds_out)[i][j] = ds[i][j];
        }
      }
    }
    if (dl0_out != nullptr) {
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
    if (dq0_out != nullptr) {
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
  void calc_derivative_wrt_group1(
    size_t ia, cvm::rvector* const dl0_1_out = nullptr,
    cvm::vector1d<cvm::rvector>* const dq0_1_out = nullptr,
    cvm::matrix2d<cvm::rvector>* const ds_1_out = nullptr) const {
      if (dl0_1_out == nullptr && dq0_1_out == nullptr) return;
      cvm::real a2x, a2y, a2z;
      // we can get rid of the helper function read_atom_coord if C++17 (constexpr) is available
      read_atom_coord(ia, m_pos2, &a2x, &a2y, &a2z);
      cvm::rvector ds_1[4][4];
      ds_1[0][0].set( a2x,  a2y,  a2z);
      ds_1[1][0].set( 0.0,  a2z, -a2y);
      ds_1[0][1] = ds_1[1][0];
      ds_1[2][0].set(-a2z,  0.0,  a2x);
      ds_1[0][2] = ds_1[2][0];
      ds_1[3][0].set( a2y, -a2x,  0.0);
      ds_1[0][3] = ds_1[3][0];
      ds_1[1][1].set( a2x, -a2y, -a2z);
      ds_1[2][1].set( a2y,  a2x,  0.0);
      ds_1[1][2] = ds_1[2][1];
      ds_1[3][1].set( a2z,  0.0,  a2x);
      ds_1[1][3] = ds_1[3][1];
      ds_1[2][2].set(-a2x,  a2y, -a2z);
      ds_1[3][2].set( 0.0,  a2z,  a2y);
      ds_1[2][3] = ds_1[3][2];
      ds_1[3][3].set(-a2x, -a2y,  a2z);
      calc_derivative_impl(ds_1, dl0_1_out, dq0_1_out, ds_1_out);
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
  void calc_derivative_wrt_group2(
    size_t ia, cvm::rvector* const dl0_2_out = nullptr,
    cvm::vector1d<cvm::rvector>* const dq0_2_out = nullptr,
    cvm::matrix2d<cvm::rvector>* const ds_2_out = nullptr) const {
    if (dl0_2_out == nullptr && dq0_2_out == nullptr) return;
    cvm::real a1x, a1y, a1z;
    // we can get rid of the helper function read_atom_coord if C++17 (constexpr) is available
    read_atom_coord(ia, m_pos1, &a1x, &a1y, &a1z);
    cvm::rvector ds_2[4][4];
    ds_2[0][0].set( a1x,  a1y,  a1z);
    ds_2[1][0].set( 0.0, -a1z,  a1y);
    ds_2[0][1] = ds_2[1][0];
    ds_2[2][0].set( a1z,  0.0, -a1x);
    ds_2[0][2] = ds_2[2][0];
    ds_2[3][0].set(-a1y,  a1x,  0.0);
    ds_2[0][3] = ds_2[3][0];
    ds_2[1][1].set( a1x, -a1y, -a1z);
    ds_2[2][1].set( a1y,  a1x,  0.0);
    ds_2[1][2] = ds_2[2][1];
    ds_2[3][1].set( a1z,  0.0,  a1x);
    ds_2[1][3] = ds_2[3][1];
    ds_2[2][2].set(-a1x,  a1y, -a1z);
    ds_2[3][2].set( 0.0,  a1z,  a1y);
    ds_2[2][3] = ds_2[3][2];
    ds_2[3][3].set(-a1x, -a1y,  a1z);
    calc_derivative_impl(ds_2, dl0_2_out, dq0_2_out, ds_2_out);
  }
};

/*! @brief  Function for debugging gradients (allow using either
 *          std::vector<cvm::atom_pos> or std::vector<cvm::atom> for
 *          pos1 and pos2)
 *  @param[in]  pos1  Atom positions of group 1
 *  @param[in]  pos2  Atom positions of group 2
 */
template<typename T1, typename T2>
void debug_gradients(
  cvm::rotation &rot,
  const std::vector<T1> &pos1,
  const std::vector<T2> &pos2) {
  static_assert(std::is_same<T1, cvm::atom_pos>::value || std::is_same<T1, cvm::atom>::value, "");
  static_assert(std::is_same<T2, cvm::atom_pos>::value || std::is_same<T2, cvm::atom>::value, "");
  // eigenvalues and eigenvectors
  cvm::real const L0 = rot.S_eigval[0];
  cvm::real const L1 = rot.S_eigval[1];
  cvm::real const L2 = rot.S_eigval[2];
  cvm::real const L3 = rot.S_eigval[3];
  cvm::quaternion const Q0(rot.S_eigvec[0]);
  cvm::quaternion const Q1(rot.S_eigvec[1]);
  cvm::quaternion const Q2(rot.S_eigvec[2]);
  cvm::quaternion const Q3(rot.S_eigvec[3]);

  cvm::log("L0 = "+cvm::to_str(L0, cvm::cv_width, cvm::cv_prec)+
            ", Q0 = "+cvm::to_str(Q0, cvm::cv_width, cvm::cv_prec)+
            ", Q0*Q0 = "+cvm::to_str(Q0.inner(Q0), cvm::cv_width, cvm::cv_prec)+
            "\n");
  cvm::log("L1 = "+cvm::to_str(L1, cvm::cv_width, cvm::cv_prec)+
            ", Q1 = "+cvm::to_str(Q1, cvm::cv_width, cvm::cv_prec)+
            ", Q0*Q1 = "+cvm::to_str(Q0.inner(Q1), cvm::cv_width, cvm::cv_prec)+
            "\n");
  cvm::log("L2 = "+cvm::to_str(L2, cvm::cv_width, cvm::cv_prec)+
            ", Q2 = "+cvm::to_str(Q2, cvm::cv_width, cvm::cv_prec)+
            ", Q0*Q2 = "+cvm::to_str(Q0.inner(Q2), cvm::cv_width, cvm::cv_prec)+
            "\n");
  cvm::log("L3 = "+cvm::to_str(L3, cvm::cv_width, cvm::cv_prec)+
            ", Q3 = "+cvm::to_str(Q3, cvm::cv_width, cvm::cv_prec)+
            ", Q0*Q3 = "+cvm::to_str(Q0.inner(Q3), cvm::cv_width, cvm::cv_prec)+
            "\n");

  rotation_derivative<T1, T2> deriv(rot, pos1, pos2);
  cvm::rvector dl0_2;
  cvm::vector1d<cvm::rvector> dq0_2(4);
  cvm::matrix2d<cvm::rvector> ds_2;
#ifdef COLVARS_LAMMPS
    MathEigen::Jacobi<cvm::real,
                      cvm::real[4],
                      cvm::real[4][4]> *ecalc =
        reinterpret_cast<MathEigen::Jacobi<cvm::real,
                                           cvm::real[4],
                                           cvm::real[4][4]> *>(rot.jacobi);
#endif
  deriv.prepare_derivative(rotation_derivative_dldq::use_dl | rotation_derivative_dldq::use_dq);
  cvm::real S_new[4][4];
  cvm::real S_new_eigval[4];
  cvm::real S_new_eigvec[4][4];
  for (size_t ia = 0; ia < pos2.size(); ++ia) {
    // cvm::real const &a1x = pos1[ia].x;
    // cvm::real const &a1y = pos1[ia].y;
    // cvm::real const &a1z = pos1[ia].z;
    deriv.calc_derivative_wrt_group2(ia, &dl0_2, &dq0_2, &ds_2);
    // make an infitesimal move along each cartesian coordinate of
    // this atom, and solve again the eigenvector problem
    for (size_t comp = 0; comp < 3; comp++) {
      std::memcpy(S_new, rot.S_backup, sizeof(cvm::real) * 4 * 4);
      std::memset(S_new_eigval, 0, sizeof(cvm::real) * 4);
      std::memset(S_new_eigvec, 0, sizeof(cvm::real) * 4 * 4);
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          S_new[i][j] +=
            colvarmodule::debug_gradients_step_size * ds_2[i][j][comp];
        }
      }
#ifdef COLVARS_LAMMPS
      ecalc->Diagonalize(S_new, S_new_eigval, S_new_eigvec);
#else
      NR::diagonalize_matrix(S_new, S_new_eigval, S_new_eigvec);
#endif
      cvm::real const &L0_new = S_new_eigval[0];
      cvm::quaternion const Q0_new(S_new_eigvec[0]);

      cvm::real const DL0 = (dl0_2[comp]) * colvarmodule::debug_gradients_step_size;
      cvm::quaternion const DQ0(dq0_2[0][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[1][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[2][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[3][comp] * colvarmodule::debug_gradients_step_size);

      cvm::log(  "|(l_0+dl_0) - l_0^new|/l_0 = "+
                cvm::to_str(cvm::fabs(L0+DL0 - L0_new)/L0, cvm::cv_width, cvm::cv_prec)+
                ", |(q_0+dq_0) - q_0^new| = "+
                cvm::to_str((Q0+DQ0 - Q0_new).norm(), cvm::cv_width, cvm::cv_prec)+
                "\n");
    }
  }
}

#endif // COLVAR_ROTATION_DERIVATIVE
