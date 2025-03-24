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

#ifndef BOSONIC_EXCHANGE_H
#define BOSONIC_EXCHANGE_H

#include "pointers.h"

namespace LAMMPS_NS {

class BosonicExchange : protected Pointers {
 public:
  BosonicExchange(class LAMMPS *, int nbosons, int np, int bead_num, bool mic,
                  bool beta_convention);
  ~BosonicExchange();

  void prepare_with_coordinates(const double *x, const double *x_prev, const double *x_next,
                                double beta, double spring_constant);

  double get_potential() const;
  double get_bead_spring_energy() const;

  void spring_force(double **f) const;

  double prim_estimator();

 private:
  void evaluate_cycle_energies();
  void diff_two_beads(const double *x1, int l1, const double *x2, int l2, double diff[3]) const;
  double get_interior_bead_spring_energy() const;
  double distance_squared_two_beads(const double *x1, int l1, const double *x2, int l2) const;
  double get_Enk(int m, int k) const;
  void set_Enk(int m, int k, double val);
  void evaluate_connection_probabilities();
  void spring_force_last_bead(double **f) const;
  void spring_force_first_bead(double **f) const;
  void spring_force_interior_bead(double **f) const;
  void Evaluate_VBn();
  void Evaluate_V_backwards();

  const int nbosons;
  const int np;
  const int bead_num;
  const bool apply_minimum_image;
  const bool physical_beta_convention;

  double spring_constant;
  double beta;
  const double *x;
  const double *x_prev;
  const double *x_next;

  double *E_kn;
  double *V;
  double *V_backwards;
  double *connection_probabilities;

  double *temp_nbosons_array;
};
}    // namespace LAMMPS_NS
#endif
