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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(rheo/property/atom,ComputeRHEOPropertyAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_PROPERTY_ATOM_H
#define LMP_COMPUTE_RHEO_PROPERTY_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOPropertyAtom : public Compute {
 public:
  ComputeRHEOPropertyAtom(class LAMMPS *, int, char **);
  ~ComputeRHEOPropertyAtom() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nvalues;
  int nmax;
  double *buf;

  typedef void (ComputeRHEOPropertyAtom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_rho(int);
  void pack_drho(int);
  void pack_temperature(int);
  void pack_heatflow(int);
  void pack_status(int);
  void pack_phase(int);
  void pack_surface(int);
  void pack_r_surface(int);
  void pack_divr_surface(int);
  void pack_nx_surface(int);
  void pack_ny_surface(int);
  void pack_nz_surface(int);
  void pack_coordination(int);
  void pack_viscosity(int);
  void pack_pressure(int);
  void pack_conductivity(int);
  void pack_cv(int);
  void pack_vx_shift(int);
  void pack_vy_shift(int);
  void pack_vz_shift(int);

};

}    // namespace LAMMPS_NS

#endif
#endif
