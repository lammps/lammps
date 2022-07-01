// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "meam_kokkos.h"

template<class DeviceType>
void MEAMKokkos<DeviceType>::meam_setup_done(double* cutmax)
{
  MEAM::meam_setup_done(cutmax);

  MemKK::realloc_kokkos(d_phir, "pair:phir", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar, "pair:phirar", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar1, "pair:phirar1", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar2, "pair:phirar2", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar3, "pair:phirar3", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar4, "pair:phirar4", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar5, "pair:phirar5", (neltypes * (neltypes + 1)) / 2, nr);
  MemKK::realloc_kokkos(d_phirar6, "pair:phirar6", (neltypes * (neltypes + 1)) / 2, nr);

  auto h_phir = Kokkos::create_mirror_view(d_phir);
  auto h_phirar = Kokkos::create_mirror_view(d_phirar);
  auto h_phirar1 = Kokkos::create_mirror_view(d_phirar1);
  auto h_phirar2 = Kokkos::create_mirror_view(d_phirar2);
  auto h_phirar3 = Kokkos::create_mirror_view(d_phirar3);
  auto h_phirar4 = Kokkos::create_mirror_view(d_phirar4);
  auto h_phirar5 = Kokkos::create_mirror_view(d_phirar5);
  auto h_phirar6 = Kokkos::create_mirror_view(d_phirar6);

  for (int i = 0; i <(neltypes * (neltypes + 1)) / 2; i++)
    for(int j = 0; j < nr; j++) {
      h_phir(i,j) = phir[i][j];
      h_phirar(i,j) = phirar[i][j];
      h_phirar1(i,j) = phirar1[i][j];
      h_phirar2(i,j) = phirar2[i][j];
      h_phirar3(i,j) = phirar3[i][j];
      h_phirar4(i,j) = phirar4[i][j];
      h_phirar5(i,j) = phirar5[i][j];
      h_phirar6(i,j) = phirar6[i][j];
    }

  Kokkos::deep_copy(d_phir,h_phir);
  Kokkos::deep_copy(d_phirar,h_phirar);
  Kokkos::deep_copy(d_phirar1,h_phirar1);
  Kokkos::deep_copy(d_phirar2,h_phirar2);
  Kokkos::deep_copy(d_phirar3,h_phirar3);
  Kokkos::deep_copy(d_phirar4,h_phirar4);
  Kokkos::deep_copy(d_phirar5,h_phirar5);
  Kokkos::deep_copy(d_phirar6,h_phirar6);
}
