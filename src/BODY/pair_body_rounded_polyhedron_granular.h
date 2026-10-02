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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(body/rounded/polyhedron/granular,PairBodyRoundedPolyhedronGranular);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_GRANULAR_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_GRANULAR_H

#include "pair_body_rounded_polyhedron.h"

#include <vector>

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class PairBodyRoundedPolyhedronGranular : public PairBodyRoundedPolyhedron {
 public:
  PairBodyRoundedPolyhedronGranular(class LAMMPS *);
  ~PairBodyRoundedPolyhedronGranular() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void transfer_history(double *, double *, int, int) override;

  double normal_force(int ibody, int jbody, int kind, double R, int itype, int jtype, double &fe,
                      double &fc) override;

 protected:
  std::vector<Granular_NS::GranularModel *> models;    // all models, owned
  std::vector<double> energy_factor;    // elastic energy = energy_factor * Fne * delta
  int **model_index;                    // index of the model of each pair of types

  Granular_NS::GranularModel *type_model(int itype, int jtype) const
  {
    return models[model_index[itype][jtype]];
  }
  Granular_NS::GranularModel *add_model();
  void check_model(int m);
  static double effective_radius(int kind, double radi, double radj);
  void contact_geometry(Granular_NS::GranularModel *, int ibody, int jbody, int kind, double R);

  void contact_damping(int ibody, int jbody, int kind, double R, double w, double vnnr,
                       const double *n, const double *vn, const double *vt,
                       double *fdiss) override;
  void tangential_spring(int ibody, int jbody, const double *n, const double *vt, double fne,
                         Scratch &s, double *fs) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
