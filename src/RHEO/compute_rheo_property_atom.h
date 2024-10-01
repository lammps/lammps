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
  void init() override;
  void setup() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nvalues, nmax;
  int pressure_flag, thermal_flag, interface_flag;
  int surface_flag, shift_flag, shell_flag;
  int *avec_index;
  int *col_index, *col_t_index;
  double *buf;

  typedef void (ComputeRHEOPropertyAtom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_phase(int);
  void pack_status(int);
  void pack_chi(int);
  void pack_surface(int);
  void pack_surface_r(int);
  void pack_surface_divr(int);
  void pack_surface_n(int);
  void pack_coordination(int);
  void pack_cv(int);
  void pack_shift_v(int);
  void pack_gradv(int);
  void pack_pressure(int);
  void pack_viscous_stress(int);
  void pack_total_stress(int);
  void pack_nbond_shell(int);
  void pack_atom_style(int);

  int add_vector_component(char *, int, FnPtrPack);
  int add_tensor_component(char *, int, FnPtrPack);

  class FixRHEO *fix_rheo;
  class FixRHEOPressure *fix_pressure;
  class FixRHEOThermal *fix_thermal;
  class FixRHEOOxidation *fix_oxidation;
  class ComputeRHEOInterface *compute_interface;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOSurface *compute_surface;
  class ComputeRHEOVShift *compute_vshift;
  class ComputeRHEOGrad *compute_grad;
};

}    // namespace LAMMPS_NS

#endif
#endif
