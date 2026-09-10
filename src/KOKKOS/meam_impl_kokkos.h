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

}

#include "meam_setup_done_kokkos.h"
#include "meam_funcs_kokkos.h"
#include "meam_dens_init_kokkos.h"
#include "meam_dens_final_kokkos.h"
#include "meam_force_kokkos.h"

