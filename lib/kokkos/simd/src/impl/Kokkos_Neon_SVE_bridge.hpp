// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEON_SVE_BRIDGE_HPP
#define KOKKOS_NEON_SVE_BRIDGE_HPP

#include <arm_neon.h>
#include <arm_sve.h>

#include <Kokkos_Macros.hpp>

/**
 * See https://clang.llvm.org/doxygen/arm__neon__sve__bridge_8h_source.html
 * for documentation on the Neon-SVE bridge
 */

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static uint64x2_t svget_neonq_u64(
    svuint64_t z) {
  uint64x2_t res;

  svbool_t pg0 = svpfirst(svptrue_b64(), svpfalse());
  svbool_t pg1 = svpnext_b64(pg0, pg0);

  res[0] = svlastb(pg0, z);
  res[1] = svlastb(pg1, z);

  return res;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static int32x4_t svget_neonq_s32(
    svint32_t z) {
  int32x4_t res;

  svbool_t pg0 = svpfirst(svptrue_b32(), svpfalse());
  svbool_t pg1 = svpnext_b32(pg0, pg0);
  svbool_t pg2 = svpnext_b32(pg1, pg1);
  svbool_t pg3 = svpnext_b32(pg2, pg2);

  res[0] = svlastb(pg0, z);
  res[1] = svlastb(pg1, z);
  res[2] = svlastb(pg2, z);
  res[3] = svlastb(pg3, z);

  return res;
}

#endif  // KOKKOS_NEON_SVE_BRIDGE_HPP
