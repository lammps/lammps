/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    Contributing authors:
             Rodrigo Freitas (UC Berkeley) - rodrigof@berkeley.edu
             Mark Asta (UC Berkeley) - mdasta@berkeley.edu
             Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ti/spring,FixTISpring);
// clang-format on
#else

#ifndef LMP_FIX_TI_SPRING_H
#define LMP_FIX_TI_SPRING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTISpring : public Fix {
 public:
  FixTISpring(class LAMMPS *, int, char **);
  ~FixTISpring() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  void initial_integrate(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  double switch_func(double);     // Switching function.
  double dswitch_func(double);    // Switching function derivative.

  double k;              // Spring constant.
  double espring;        // Springs energies.
  double **xoriginal;    // Original coords of atoms.
  double lambda;         // Coupling parameter.
  double dlambda;        // Lambda variation with t.
  double linfo[2];       // Current lambda status.
  bigint t_switch;       // Total switching steps.
  bigint t_equil;        // Equilibration time.
  bigint t0;             // Initial time.
  int sf;                // Switching function option.
  int nlevels_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
