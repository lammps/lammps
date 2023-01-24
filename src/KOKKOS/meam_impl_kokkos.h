// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naga Vydyanathan (NVIDIA), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "memory_kokkos.h"
#include "meam_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
MEAMKokkos<DeviceType>::MEAMKokkos(Memory *mem) : MEAM(mem)
{
  d_errorflag = typename AT::t_int_scalar("meam:errorflag");
}

template<class DeviceType>
MEAMKokkos<DeviceType>::~MEAMKokkos()
{
  if (copymode) return;

  MemoryKokkos *memoryKK = (MemoryKokkos *)memory;

  memoryKK->destroy_kokkos(k_rho,rho);
  memoryKK->destroy_kokkos(k_rho0,rho0);
  memoryKK->destroy_kokkos(k_rho1,rho1);
  memoryKK->destroy_kokkos(k_rho2,rho2);
  memoryKK->destroy_kokkos(k_rho3,rho3);
  memoryKK->destroy_kokkos(k_frhop,frhop);
  memoryKK->destroy_kokkos(k_gamma,gamma);
  memoryKK->destroy_kokkos(k_dgamma1,dgamma1);
  memoryKK->destroy_kokkos(k_dgamma2,dgamma2);
  memoryKK->destroy_kokkos(k_dgamma3,dgamma3);
  memoryKK->destroy_kokkos(k_arho2b,arho2b);

  memoryKK->destroy_kokkos(k_arho1,arho1);
  memoryKK->destroy_kokkos(k_arho2,arho2);
  memoryKK->destroy_kokkos(k_arho3,arho3);
  memoryKK->destroy_kokkos(k_arho3b,arho3b);
  memoryKK->destroy_kokkos(k_t_ave,t_ave);
  memoryKK->destroy_kokkos(k_tsq_ave,tsq_ave);

  memoryKK->destroy_kokkos(k_scrfcn,scrfcn);
  memoryKK->destroy_kokkos(k_dscrfcn,dscrfcn);
  memoryKK->destroy_kokkos(k_fcpair,fcpair);

  // msmeam

  memoryKK->destroy_kokkos(k_arho2mb, arho2mb);
  memoryKK->destroy_kokkos(k_arho1m, arho1m);
  memoryKK->destroy_kokkos(k_arho2m, arho2m);
  memoryKK->destroy_kokkos(k_arho3m, arho3m);
  memoryKK->destroy_kokkos(k_arho3mb, arho3mb);
}

#include "meam_setup_done_kokkos.h"
#include "meam_funcs_kokkos.h"
#include "meam_dens_init_kokkos.h"
#include "meam_dens_final_kokkos.h"
#include "meam_force_kokkos.h"

