// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include <impl/Kokkos_Tools.hpp>

int main(int argc, char* argv[]) {
  Kokkos::Tools::initialize(argc, argv);
  Kokkos::Tools::pushRegion(
      "The unanimous Declaration of the thirteen united States of America, "
      "When in the Course of human events, it becomes necessary for one people "
      "to dissolve the political bands which have connected them with another, "
      "and to assume among the powers of the earth, the separate and equal "
      "station to which the Laws of Nature and of Nature's God entitle them, a "
      "decent respect to the opinions of mankind requires that they should "
      "declare the causes which impel them to the separation.");
  Kokkos::Tools::popRegion();
  Kokkos::Tools::finalize();
}
