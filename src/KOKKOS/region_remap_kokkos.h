/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// clang-format off
#ifndef LMP_REGION_REMAP_KOKKOS_H
#define LMP_REGION_REMAP_KOKKOS_H

#include "domain.h"

#include "kokkos_type.h"

namespace LAMMPS_NS {

// device-side equivalent of the Domain::remap() call in Region::match()
// and Region::surface(): positions must be wrapped back into the periodic
// box before region matching.  The struct is a plain-data member of the
// region classes and is copied to the device with the region functor, so
// capture() must be called on the host (from prematch() and
// match_all_kokkos()) before launching any kernel that uses remap().

struct RegionRemapKokkos {
  int triclinic;
  int xperiodic, yperiodic, zperiodic;
  double boxlo[3], boxhi[3], prd[3];
  double boxlo_lamda[3], boxhi_lamda[3], prd_lamda[3];
  double h[6], h_inv[6];

  void capture(Domain *domain)
  {
    triclinic = domain->triclinic;
    xperiodic = domain->xperiodic;
    yperiodic = domain->yperiodic;
    zperiodic = domain->zperiodic;
    for (int i = 0; i < 3; i++) {
      boxlo[i] = domain->boxlo[i];
      boxhi[i] = domain->boxhi[i];
      prd[i] = domain->prd[i];
      boxlo_lamda[i] = domain->boxlo_lamda[i];
      boxhi_lamda[i] = domain->boxhi_lamda[i];
      prd_lamda[i] = domain->prd_lamda[i];
    }
    for (int i = 0; i < 6; i++) {
      h[i] = domain->h[i];
      h_inv[i] = domain->h_inv[i];
    }
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void remap(double &x, double &y, double &z) const
  {
    double coord[3];

    if (triclinic == 0) {
      coord[0] = x; coord[1] = y; coord[2] = z;
    } else {
      const double dx = x - boxlo[0];
      const double dy = y - boxlo[1];
      const double dz = z - boxlo[2];
      coord[0] = h_inv[0]*dx + h_inv[5]*dy + h_inv[4]*dz;
      coord[1] = h_inv[1]*dy + h_inv[3]*dz;
      coord[2] = h_inv[2]*dz;
    }

    const double *lo = (triclinic == 0) ? boxlo : boxlo_lamda;
    const double *hi = (triclinic == 0) ? boxhi : boxhi_lamda;
    const double *period = (triclinic == 0) ? prd : prd_lamda;

    if (xperiodic) {
      while (coord[0] < lo[0]) coord[0] += period[0];
      while (coord[0] >= hi[0]) coord[0] -= period[0];
      if (coord[0] < lo[0]) coord[0] = lo[0];
    }

    if (yperiodic) {
      while (coord[1] < lo[1]) coord[1] += period[1];
      while (coord[1] >= hi[1]) coord[1] -= period[1];
      if (coord[1] < lo[1]) coord[1] = lo[1];
    }

    if (zperiodic) {
      while (coord[2] < lo[2]) coord[2] += period[2];
      while (coord[2] >= hi[2]) coord[2] -= period[2];
      if (coord[2] < lo[2]) coord[2] = lo[2];
    }

    if (triclinic == 0) {
      x = coord[0]; y = coord[1]; z = coord[2];
    } else {
      x = h[0]*coord[0] + h[5]*coord[1] + h[4]*coord[2] + boxlo[0];
      y = h[1]*coord[1] + h[3]*coord[2] + boxlo[1];
      z = h[2]*coord[2] + boxlo[2];
    }
  }
};

}

#endif
