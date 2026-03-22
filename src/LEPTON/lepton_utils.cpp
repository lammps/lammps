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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "lepton_utils.h"

#include "error.h"
#include "input.h"
#include "lammps.h"
#include "pair_zbl_const.h"
#include "variable.h"

#include <cctype>
#include <cmath>
#include <exception>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace PairZBLConstants;

namespace Lepton {
class DerivativeException : public std::exception {
  std::string message;

 public:
  // remove unused default constructor
  DerivativeException() = delete;

  explicit DerivativeException(int deg, const std::string &fn, const std::string &vn)
  {
    message = fmt::format("Order {} derivative of function {} in {} is not supported", deg, fn, vn);
  }
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
};
}    // namespace Lepton

double Lepton::ZBLFunction::evaluate(const double *args) const
{
  const double zi = args[0];
  const double zj = args[1];
  const double r = args[2];

  const double rbya = r * (pow(zi, pzbl) + pow(zj, pzbl)) / (a0 * angstrom);
  return zi * zj * qqr2e * qelectron * qelectron / r *
      (c4 * exp(-d4 * rbya) + c3 * exp(-d3 * rbya) + c2 * exp(-d2 * rbya) + c1 * exp(-d1 * rbya));
}

double Lepton::ZBLFunction::evaluateDerivative(const double *args, const int *order) const
{
  if (order[0] > 0) throw DerivativeException(order[0], "zbl()", "'zi'");
  if (order[1] > 0) throw DerivativeException(order[0], "zbl()", "'zj'");
  if (order[2] > 1) throw DerivativeException(order[0], "zbl()", "'r'");

  if (order[2] == 1) {
    const double zi = args[0];
    const double zj = args[1];
    const double r = args[2];

    const double ainv = (pow(zi, pzbl) + pow(zj, pzbl)) / (a0 * angstrom);
    const double e1 = exp(-d1 * ainv * r);
    const double e2 = exp(-d2 * ainv * r);
    const double e3 = exp(-d3 * ainv * r);
    const double e4 = exp(-d4 * ainv * r);

    const double sum1 = c1 * e1 + c2 * e2 + c3 * e3 + c4 * e4;
    const double sum2 = ainv * (-c1 * d1 * e1 - c2 * d2 * e2 - c3 * d3 * e3 - c4 * d4 * e4);
    return (zi * zj * qqr2e * qelectron * qelectron) * (sum2 - sum1 / r) / r;
  }
  return 0.0;
}

namespace LeptonUtils {
class VariableException : public std::exception {
  std::string message;

 public:
  // remove unused default constructor
  VariableException() = delete;

  explicit VariableException(const std::string &var, const std::string &expr)
  {
    message = fmt::format("Variable {} in expression {} does not exist", var, expr);
  }
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }
};
}    // namespace LeptonUtils

/// remove whitespace and quotes from expression string
std::string LeptonUtils::condense(const std::string &in)
{
  std::string out;
  for (const auto &c : in)
    if (!isspace(c) && (c != '"') && (c != '\'')) out.push_back(c);
  return out;
}

/// substitute variable references with their values
std::string LeptonUtils::substitute(const std::string &in, LAMMPS_NS::LAMMPS *lmp)
{
  std::string output, name;
  bool in_var = false;
  char hold = ' ';
  auto *variable = lmp->input->variable;

  for (const auto &c : in) {
    if (in_var) {
      if (isalnum(c) || (c == '_')) {
        name.push_back(c);
      } else {
        in_var = false;
        const char *val = variable->retrieve(name.c_str());
        if (val)
          output.append(val);
        else
          lmp->error->all(FLERR, "Variable {} in expression {} does not exist", name, in);

        output.push_back(c);
      }
    } else {
      if (hold == 'v') {
        if (c == '_') {
          in_var = true;
          hold = ' ';
          name.clear();
        } else {
          output.push_back(hold);
          hold = ' ';
          output.push_back(c);
        }
      } else {
        if (c == 'v')
          hold = c;
        else
          output.push_back(c);
      }
    }
  }
  if (in_var) {
    const char *val = variable->retrieve(name.c_str());
    if (val)
      output.append(val);
    else
      lmp->error->all(FLERR, "Variable {} in expression {} does not exist", name, in);
  }

  return output;
}
