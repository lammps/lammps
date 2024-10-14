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

#ifndef LMP_GROUP_KOKKOS_H
#define LMP_GROUP_KOKKOS_H

#include "group.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class GroupKokkos : public Group {
 public:
  GroupKokkos(class LAMMPS *);
  double mass(int);               // total mass of atoms in group
  void xcm(int, double, double *);    // center-of-mass coords of group
};

}    // namespace LAMMPS_NS

#endif
