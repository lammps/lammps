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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(rheo/thermal,AtomVecRHEOThermal);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_RHEO_THERMAL_H
#define LMP_ATOM_VEC_RHEO_THERMAL_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecRHEOThermal : virtual public AtomVec {
 public:
  AtomVecRHEOThermal(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  int property_atom(const std::string &) override;
  void pack_property_atom(int, double *, int, int) override;

 private:
  int *rheo_status;
  double *conductivity, *temperature, *heatflow, *esph;
  double *pressure, *rho, *drho, *viscosity;
};

}    // namespace LAMMPS_NS

#endif
#endif
