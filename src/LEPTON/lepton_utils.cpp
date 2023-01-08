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
#include "variable.h"

#include "fmt/args.h"

#include <cctype>
#include <exception>
#include <unordered_set>

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
  const char *what() const noexcept override { return message.c_str(); }
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
  std::string format, name;
  std::unordered_set<std::string> vars;
  bool in_var = false;
  char hold = ' ';

  for (const auto &c : in) {
    if (in_var) {
      if (isalnum(c) || (c == '_')) {
        format.push_back(c);
        name.push_back(c);
      } else {
        in_var = false;
        format.push_back('}');
        format.push_back(c);
        vars.insert(name);
      }
    } else {
      if (hold == 'v') {
        if (c == '_') {
          in_var = true;
          hold = ' ';
          name.clear();
          format.push_back('{');
        } else {
          format.push_back(hold);
          hold = ' ';
          format.push_back(c);
        }
      } else {
        if (c == 'v')
          hold = c;
        else
          format.push_back(c);
      }
    }
  }
  if (in_var) {
    format.push_back('}');
    vars.insert(name);
  }

  auto variable = lmp->input->variable;
  fmt::dynamic_format_arg_store<fmt::format_context> args;
  for (const auto &v : vars) {
    const char *val = variable->retrieve(v.c_str());
    if (val)
      args.push_back(fmt::arg(v.c_str(), val));
    else
      throw VariableException(v, in);
  }
  return fmt::vformat(format, args);
}
