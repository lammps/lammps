/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS

NPairStyle(copy/kk/device,
           NPairCopyKokkos<LMPDeviceType>,
           NP_COPY | NP_KOKKOS_DEVICE)

NPairStyle(copy/kk/host,
           NPairCopyKokkos<LMPHostType>,
           NP_COPY | NP_KOKKOS_HOST)

#else

#ifndef LMP_NPAIR_COPY_KOKKOS_H
#define LMP_NPAIR_COPY_KOKKOS_H

#include "npair.h"

namespace LAMMPS_NS {

template<class DeviceType>
class NPairCopyKokkos : public NPair {
 public:
  NPairCopyKokkos(class LAMMPS *);
  ~NPairCopyKokkos() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
