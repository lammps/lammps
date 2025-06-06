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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_CHARGE_SOLVER_H
#define LMP_CHARGE_SOLVER_H

#include "lmptype.h"
#include <vector>

namespace LAMMPS_NS {
class ChargeSolver {
 public:
  ChargeSolver();
  virtual ~ChargeSolver();
  virtual void update_solver(std::vector<tagint>, std::vector<int>) = 0;
  virtual std::vector<double> solve(std::vector<double>) = 0;
  virtual void set_elyt_pot(double *) = 0;
  virtual double get_potential(int) = 0;
  virtual double get_sb_charges(int) = 0;
  virtual double get_macro_capacitance(int, int) = 0;
  virtual double get_macro_elastance(int, int) = 0;
  virtual void buffer_and_gather(double const *, double *) = 0;
  virtual double memory_use() = 0;
  double get_mult_time();

  // for electrode/thermo
  virtual std::vector<double> compute_potentials() = 0;
  virtual double vacuum_capacitance() = 0;

  // constraint setup
  void set_constraint(double);
  void set_constraint(std::vector<double>);

 protected:
  enum class ChargeConstraint { NONE, SINGLE, GROUP };
  ChargeConstraint constraint;
  double qtotal, mult_time;
  std::vector<double> qtotal_group;
};
}    // namespace LAMMPS_NS

#endif

