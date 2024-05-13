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

#include "Lepton.h"

#include <string>

// forward declarations

namespace LAMMPS_NS {
class LAMMPS;
}

// custom function zbl(zi,zj,r)

namespace Lepton {
class ZBLFunction : public CustomFunction {
 public:
  ZBLFunction(double _qqr2e, double _angstrom, double _qelectron) :
      qqr2e(_qqr2e), angstrom(_angstrom), qelectron(_qelectron){};
  ZBLFunction() = delete;

  int getNumArguments() const override { return 3; }
  CustomFunction *clone() const override { return new ZBLFunction(qqr2e, angstrom, qelectron); }
  double evaluate(const double *) const override;
  double evaluateDerivative(const double *, const int *) const override;

 private:
  double qqr2e, angstrom, qelectron;
};
}    // namespace Lepton

// utility functions and classes

namespace LeptonUtils {

/// remove whitespace and quotes from expression string
std::string condense(const std::string &);
/// substitute LAMMPS variable references with their value
std::string substitute(const std::string &, LAMMPS_NS::LAMMPS *);

}    // namespace LeptonUtils
