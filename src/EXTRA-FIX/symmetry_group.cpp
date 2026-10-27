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

#include "symmetry_group.h"

#include "comm.h"
#include "error.h"
#include "json.h"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <unordered_map>

using namespace LAMMPS_NS;

namespace {
// tolerance for integer matrix elements and for comparing translations modulo 1
constexpr double INT_TOL = 1.0e-9;

// value of a number in the forms "1", "1/2", "0.25", or ".25"
double parse_number(const std::string &s)
{
  if (s.empty()) throw std::invalid_argument("missing number");
  const auto slash = s.find('/');
  std::size_t pos = 0;
  if (slash == std::string::npos) {
    double val = std::stod(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("invalid number '" + s + "'");
    return val;
  }
  const auto num = s.substr(0, slash);
  const auto den = s.substr(slash + 1);
  if (num.empty() || den.empty() || (num.find_first_not_of("0123456789") != std::string::npos) ||
      (den.find_first_not_of("0123456789") != std::string::npos))
    throw std::invalid_argument("invalid fraction '" + s + "'");
  const double dval = std::stod(den);
  if (dval == 0.0) throw std::invalid_argument("zero denominator in '" + s + "'");
  return std::stod(num) / dval;
}

// translation component given as a JSON number or a string like "1/2"
double get_t_value(const json &v)
{
  if (v.is_number()) return v.get<double>();
  if (v.is_string()) return parse_number(v.get<std::string>());
  throw std::invalid_argument("translation must be a number or a string like \"1/2\"");
}

double wrap01(double v)
{
  v -= std::floor(v);
  if (v >= 1.0) v -= 1.0;
  return v;
}

bool R_equal(const double A[3][3], const double B[3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (std::fabs(A[i][j] - B[i][j]) > INT_TOL) return false;
  return true;
}

bool t_equal_mod1(const double a[3], const double b[3])
{
  for (int k = 0; k < 3; ++k) {
    double d = wrap01(a[k] - b[k]);
    if (d > 0.5) d = 1.0 - d;
    if (d > INT_TOL) return false;
  }
  return true;
}

bool is_identity(const SymmOp &op)
{
  const double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const double zero[3] = {0, 0, 0};
  return R_equal(op.R, I) && t_equal_mod1(op.t, zero);
}
}    // namespace

/* ---------------------------------------------------------------------- */

SymmetryGroup::SymmetryGroup(LAMMPS *lmp) : Pointers(lmp), group_name("unknown")
{
  ncell[0] = ncell[1] = ncell[2] = 1;
}

/* ----------------------------------------------------------------------
   Convert an operator in the coordinate triplet notation of the
   International Tables for Crystallography and of CIF files, e.g.
   "-y,x-y,z+1/3", into a matrix R and a translation t.  The three comma
   separated components are the new x, y, and z.  Each is a sum of terms
   that are either an integer multiple of x, y, or z, or a number such as
   1/2 or 0.25.  Spaces are ignored and upper case X, Y, Z are accepted.
------------------------------------------------------------------------- */

SymmOp SymmetryGroup::parse_xyz(const std::string &str)
{
  SymmOp op{};
  std::string clean;
  for (char c : str)
    if (!isspace(c)) clean += static_cast<char>(tolower(c));

  std::vector<std::string> parts;
  std::size_t start = 0, comma;
  while ((comma = clean.find(',', start)) != std::string::npos) {
    parts.push_back(clean.substr(start, comma - start));
    start = comma + 1;
  }
  parts.push_back(clean.substr(start));
  if (parts.size() != 3)
    throw std::invalid_argument("operator '" + str + "' must have 3 comma separated components");

  for (int i = 0; i < 3; ++i) {
    const auto &part = parts[i];
    if (part.empty()) throw std::invalid_argument("operator '" + str + "' has an empty component");
    std::size_t pos = 0;
    bool first = true;
    while (pos < part.size()) {
      double sign = 1.0;
      if ((part[pos] == '+') || (part[pos] == '-')) {
        if (part[pos] == '-') sign = -1.0;
        ++pos;
      } else if (!first) {
        throw std::invalid_argument("operator '" + str + "': missing + or - in '" + part + "'");
      }
      first = false;

      // number: digits with an optional decimal point or fraction
      std::size_t end = pos;
      while ((end < part.size()) && (isdigit(part[end]) || (part[end] == '.') || (part[end] == '/')))
        ++end;
      const std::string num = part.substr(pos, end - pos);
      pos = end;
      if ((pos < part.size()) && (part[pos] == '*')) ++pos;

      if ((pos < part.size()) && (part[pos] >= 'x') && (part[pos] <= 'z')) {
        double coeff = 1.0;
        if (!num.empty()) {
          if (num.find_first_not_of("0123456789") != std::string::npos)
            throw std::invalid_argument("operator '" + str + "': factor '" + num +
                                        "' must be an integer");
          coeff = std::stod(num);
        }
        op.R[i][part[pos] - 'x'] += sign * coeff;
        ++pos;
      } else {
        if (num.empty() || (part[pos - 1] == '*'))
          throw std::invalid_argument("operator '" + str + "': cannot parse '" + part + "'");
        try {
          op.t[i] += sign * parse_number(num);
        } catch (std::invalid_argument &e) {
          throw std::invalid_argument("operator '" + str + "': " + e.what());
        }
      }
    }
  }
  return op;
}

/* ----------------------------------------------------------------------
   read + broadcast + parse + validate. Mirrors the molecule.cpp pattern:
   rank 0 parses the JSON file and re-serializes to UBJSON; the binary form
   is broadcast; every rank deserializes and processes its own copy.
------------------------------------------------------------------------- */

void SymmetryGroup::read(const char *filename)
{
  std::vector<std::uint8_t> ubjson;
  int ubjson_size = 0;

  if (comm->me == 0) {
    FILE *fp = std::fopen(filename, "r");
    if (!fp)
      error->one(FLERR, Error::NOLASTLINE, "Cannot open symmetry file {}{}", filename,
                 utils::getsyserror());
    try {
      json data = json::parse(fp);
      ubjson = json::to_ubjson(data);
      ubjson_size = static_cast<int>(ubjson.size());
    } catch (std::exception &e) {
      std::fclose(fp);
      error->one(FLERR, Error::NOLASTLINE, "Error parsing symmetry file {}: {}", filename,
                 e.what());
    }
    std::fclose(fp);
  }

  MPI_Bcast(&ubjson_size, 1, MPI_INT, 0, world);
  if (ubjson_size <= 0) error->all(FLERR, Error::NOLASTLINE, "Symmetry file produced no JSON data");
  if (comm->me != 0) ubjson.resize(ubjson_size);
  MPI_Bcast(ubjson.data(), ubjson_size, MPI_CHAR, 0, world);

  json data;
  try {
    data = json::from_ubjson(ubjson);
  } catch (std::exception &e) {
    error->all(FLERR, Error::NOLASTLINE, "Error deserializing broadcast symmetry data: {}",
               e.what());
  }
  ubjson.clear();

  populate(data);
  compute_inverses();
  validate_group();

  // the operators of a supercell form a group if those of the unit cell do
  if ((ncell[0] > 1) || (ncell[1] > 1) || (ncell[2] > 1)) {
    apply_supercell();
    compute_inverses();
  }
}

/* ----------------------------------------------------------------------
   Extract the operators from the JSON data.  Runs on every rank with
   identical data, so errors are raised via error->all().

   Format:
     {
       "group_name": "<name>",          // optional, default "unknown"
       "ops": [                         // required, must include the identity
         "x,y,z",                       // coordinate triplet notation, or
         { "R": [[..],[..],[..]],       // 3x3 integer matrix and
           "t": [tx,ty,tz] },           // translation, numbers or strings "1/2"
         ...
       ],
       "supercell": [nx, ny, nz]        // optional, default [1, 1, 1]
     }
------------------------------------------------------------------------- */

void SymmetryGroup::populate(const json &data)
{
  try {
    if (!data.is_object()) throw std::invalid_argument("the file must contain a JSON object");
    for (const auto &item : data.items()) {
      const auto &key = item.key();
      if ((key != "group_name") && (key != "ops") && (key != "supercell") && (comm->me == 0))
        error->warning(FLERR, "Fix symmetry: ignoring unknown entry \"{}\" in symmetry file", key);
    }

    if (data.contains("group_name")) group_name = data.at("group_name").get<std::string>();

    if (!data.contains("ops")) throw std::invalid_argument("the \"ops\" entry is missing");
    const auto &jops = data.at("ops");
    if (!jops.is_array() || jops.empty())
      throw std::invalid_argument("\"ops\" must be a non-empty list of operators");
    ops.assign(jops.size(), SymmOp{});
    for (std::size_t o = 0; o < jops.size(); ++o) {
      const auto &je = jops[o];
      const auto opname = "operator " + std::to_string(o + 1) + ": ";
      if (je.is_string()) {
        ops[o] = parse_xyz(je.get<std::string>());
        continue;
      }
      if (!je.is_object() || !je.contains("R") || !je.contains("t"))
        throw std::invalid_argument(opname + "must be a string like \"-x,-y,z\" or an object "
                                             "with \"R\" and \"t\" entries");
      const auto &jR = je.at("R");
      const auto &jt = je.at("t");
      if (!jR.is_array() || jR.size() != 3)
        throw std::invalid_argument(opname + "R must be a 3x3 array");
      if (!jt.is_array() || jt.size() != 3)
        throw std::invalid_argument(opname + "t must be an array with 3 elements");
      for (int i = 0; i < 3; ++i) {
        if (!jR[i].is_array() || jR[i].size() != 3)
          throw std::invalid_argument(opname + "R must be a 3x3 array");
        for (int j = 0; j < 3; ++j) ops[o].R[i][j] = jR[i][j].get<double>();
        ops[o].t[i] = get_t_value(jt[i]);
      }
    }

    for (std::size_t o = 0; o < ops.size(); ++o)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (std::fabs(ops[o].R[i][j] - std::round(ops[o].R[i][j])) > INT_TOL)
            throw std::invalid_argument(fmt::format(
                "operator {}: matrix elements must be integers, got {:.10g}", o + 1, ops[o].R[i][j]));

    if (data.contains("supercell")) {
      const auto &jsc = data.at("supercell");
      if (!jsc.is_array() || (jsc.size() != 3))
        throw std::invalid_argument("\"supercell\" must be a list of 3 positive integers");
      for (int i = 0; i < 3; ++i) {
        ncell[i] = jsc[i].get<int>();
        if (ncell[i] < 1)
          throw std::invalid_argument("\"supercell\" must be a list of 3 positive integers");
      }
    }
  } catch (std::exception &e) {
    error->all(FLERR, Error::NOLASTLINE, "Fix symmetry: error in symmetry file: {}", e.what());
  }

  // make the identity the first operator
  int ident = -1;
  for (int o = 0; o < static_cast<int>(ops.size()); ++o)
    if (is_identity(ops[o])) {
      ident = o;
      break;
    }
  if (ident < 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix symmetry: the operators in the symmetry file must include the identity x,y,z");
  std::swap(ops[0], ops[ident]);
}

/* ----------------------------------------------------------------------
   The operators are given for a unit cell and the simulation box consists
   of ncell[0] x ncell[1] x ncell[2] unit cells.  In box coordinates s = s_u/n
   operator (R,t) becomes R_b = N^-1 R N, t_b = N^-1 t with N = diag(ncell).
   Since every translation by a unit cell is a symmetry of the box, each
   operator is combined with all of them.
------------------------------------------------------------------------- */

void SymmetryGroup::apply_supercell()
{
  std::vector<SymmOp> cellops;
  for (int o = 0; o < static_cast<int>(ops.size()); ++o) {
    SymmOp op{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        op.R[i][j] = ops[o].R[i][j] * ncell[j] / ncell[i];
        if (std::fabs(op.R[i][j] - std::round(op.R[i][j])) > INT_TOL)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix symmetry: operator {} is not compatible with a supercell of {} x {} x {} "
                     "unit cells",
                     o + 1, ncell[0], ncell[1], ncell[2]);
        op.R[i][j] = std::round(op.R[i][j]);
      }
    for (int a = 0; a < ncell[0]; ++a)
      for (int b = 0; b < ncell[1]; ++b)
        for (int c = 0; c < ncell[2]; ++c) {
          const int shift[3] = {a, b, c};
          for (int i = 0; i < 3; ++i) op.t[i] = (ops[o].t[i] + shift[i]) / ncell[i];
          cellops.push_back(op);
        }
  }
  ops = cellops;
}

/* ----------------------------------------------------------------------
   invert each op.R into op.Rinv via the adjugate / det formula. R is an
   integer-valued lattice operator with |det R| = 1, so this is exact.
------------------------------------------------------------------------- */

void SymmetryGroup::compute_inverses()
{
  for (int o = 0; o < static_cast<int>(ops.size()); ++o) {
    auto &op = ops[o];
    double(&R)[3][3] = op.R;
    double det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
        R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
        R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
    if (std::fabs(std::fabs(det) - 1.0) > INT_TOL)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: operator {} is not a symmetry operator: the determinant of its "
                 "matrix is {:.6g} instead of 1 or -1",
                 o + 1, det);
    double inv_det = 1.0 / det;
    op.Rinv[0][0] = (R[1][1] * R[2][2] - R[1][2] * R[2][1]) * inv_det;
    op.Rinv[0][1] = (R[0][2] * R[2][1] - R[0][1] * R[2][2]) * inv_det;
    op.Rinv[0][2] = (R[0][1] * R[1][2] - R[0][2] * R[1][1]) * inv_det;
    op.Rinv[1][0] = (R[1][2] * R[2][0] - R[1][0] * R[2][2]) * inv_det;
    op.Rinv[1][1] = (R[0][0] * R[2][2] - R[0][2] * R[2][0]) * inv_det;
    op.Rinv[1][2] = (R[0][2] * R[1][0] - R[0][0] * R[1][2]) * inv_det;
    op.Rinv[2][0] = (R[1][0] * R[2][1] - R[1][1] * R[2][0]) * inv_det;
    op.Rinv[2][1] = (R[0][1] * R[2][0] - R[0][0] * R[2][1]) * inv_det;
    op.Rinv[2][2] = (R[0][0] * R[1][1] - R[0][1] * R[1][0]) * inv_det;
  }
}

/* ----------------------------------------------------------------------
   The operators must be distinct and form a group: the product of any
   two operators must be in the list, up to a lattice translation.
   Operators are looked up through a hash of the matrix and the translation
   rounded to a grid, with a search of the whole list as fallback for
   translations close to a grid boundary.
------------------------------------------------------------------------- */

namespace {
std::size_t op_hash(const double R[3][3], const double t[3])
{
  std::size_t h = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) h = h * 31 + static_cast<std::size_t>(std::lround(R[i][j]) + 7);
  for (int i = 0; i < 3; ++i)
    h = h * 65599 + static_cast<std::size_t>(std::llround(wrap01(t[i]) * 65536.0) & 0xffff);
  return h;
}
}    // namespace

void SymmetryGroup::validate_group()
{
  const int N = static_cast<int>(ops.size());
  std::unordered_multimap<std::size_t, int> index;

  auto find_op = [&](const double R[3][3], const double t[3]) {
    auto range = index.equal_range(op_hash(R, t));
    for (auto it = range.first; it != range.second; ++it)
      if (R_equal(R, ops[it->second].R) && t_equal_mod1(t, ops[it->second].t)) return it->second;
    for (int c = 0; c < static_cast<int>(index.size()); ++c)
      if (R_equal(R, ops[c].R) && t_equal_mod1(t, ops[c].t)) return c;
    return -1;
  };

  for (int a = 0; a < N; ++a) {
    const int b = find_op(ops[a].R, ops[a].t);
    if (b >= 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: operators {} and {} in the symmetry file are the same", b + 1, a + 1);
    index.emplace(op_hash(ops[a].R, ops[a].t), a);
  }

  for (int a = 0; a < N; ++a) {
    for (int b = 0; b < N; ++b) {
      double Rab[3][3];
      double tab[3];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          double s = 0.0;
          for (int k = 0; k < 3; ++k) s += ops[a].R[i][k] * ops[b].R[k][j];
          Rab[i][j] = s;
        }
      for (int i = 0; i < 3; ++i) {
        double s = ops[a].t[i];
        for (int k = 0; k < 3; ++k) s += ops[a].R[i][k] * ops[b].t[k];
        tab[i] = s;
      }
      if (find_op(Rab, tab) < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: the operators do not form a group: the combination of operators "
                   "{} and {} is missing",
                   a + 1, b + 1);
    }
  }
}
