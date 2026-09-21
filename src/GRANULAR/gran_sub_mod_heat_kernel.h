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

/* ----------------------------------------------------------------------
   Stateless kernels for the granular heat conduction sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_HEAT_KERNEL_H
#define LMP_GRAN_SUB_MOD_HEAT_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

template <class T> struct GranHeatParams {
  int model;
  T coeff;    // conductivity (radius) or heat transfer coefficient (area)
};

// GranSubModHeatRadius::calculate_heat()

template <class T>
GRAN_KERNEL_FN T gran_heat_radius(const T conductivity, const T contact_radius, const T Ti,
                                  const T Tj)
{
  return (T) 2.0 * conductivity * contact_radius * (Tj - Ti);
}

// GranSubModHeatArea::calculate_heat()

template <class T>
GRAN_KERNEL_FN T gran_heat_area(const T heat_transfer_coeff, const T contact_radius, const T Ti,
                                const T Tj)
{
  return heat_transfer_coeff * (T) MathConst::MY_PI * contact_radius * contact_radius * (Tj - Ti);
}

/* ---------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN T gran_heat_flow(const GranHeatParams<T> &p, const T contact_radius, const T Ti,
                                const T Tj)
{
  switch (p.model) {
    case GRAN_HEAT_RADIUS:
      return gran_heat_radius(p.coeff, contact_radius, Ti, Tj);
    case GRAN_HEAT_AREA:
      return gran_heat_area(p.coeff, contact_radius, Ti, Tj);
    default:
      return (T) 0.0;
  }
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif
