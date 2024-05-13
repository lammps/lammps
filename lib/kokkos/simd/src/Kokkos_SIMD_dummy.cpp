//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

// This file is needed in order to get the linker language
// for the header only submodule.
// While we set the language properties in our normal cmake
// path it does not get set in the Trilinos environment.
// Furthermore, setting LINKER_LANGUAGE is only supported
// in CMAKE 3.19 and up.
void KOKKOS_SIMD_SRC_DUMMY_PREVENT_LINK_ERROR() {}
