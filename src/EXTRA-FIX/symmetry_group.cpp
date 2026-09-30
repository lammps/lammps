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
#include "math_eigen.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>

using namespace LAMMPS_NS;

namespace {
// Tolerances for validating that a parsed R matrix is an integer-valued
// lattice operator and that the group is closed modulo lattice translations.
constexpr double R_INT_TOL = 1.0e-9;
constexpr double T_MOD1_TOL = 1.0e-9;

// Parse "1/2", "-2/3", or a decimal string as a double. Used for the
// translation parts of operators (in the IT tables, always small rationals).
double parse_rational_str(const std::string &s)
{
  auto slash = s.find('/');
  if (slash == std::string::npos) return std::stod(s);
  double num = std::stod(s.substr(0, slash));
  double den = std::stod(s.substr(slash + 1));
  if (den == 0.0) throw std::runtime_error("SymmetryGroup: zero denominator in '" + s + "'");
  return num / den;
}

// Accept either a JSON number (passed through) or a JSON string (parsed as
// a rational like "1/2"). Used for op.t components.
double get_t_value(const json &v)
{
  if (v.is_number()) return v.get<double>();
  if (v.is_string()) return parse_rational_str(v.get<std::string>());
  throw std::runtime_error(
      "SymmetryGroup: translation component must be a number or rational string");
}

LatticeFamily parse_lattice(const std::string &s)
{
  if (s == "triclinic") return LATTICE_TRICLINIC;
  if (s == "monoclinic") return LATTICE_MONOCLINIC;
  if (s == "orthorhombic") return LATTICE_ORTHORHOMBIC;
  if (s == "tetragonal") return LATTICE_TETRAGONAL;
  if (s == "hexagonal") return LATTICE_HEXAGONAL;
  if (s == "trigonal") return LATTICE_TRIGONAL;
  if (s == "cubic") return LATTICE_CUBIC;
  throw std::runtime_error("SymmetryGroup: unknown lattice family '" + s + "'");
}
}    // namespace

/* ---------------------------------------------------------------------- */

SymmetryGroup::SymmetryGroup(LAMMPS *lmp) :
    Pointers(lmp), group_name("unknown"), lattice_family(LATTICE_TRICLINIC)
{
}

/* ----------------------------------------------------------------------
   read + broadcast + parse + validate. Mirrors the molecule.cpp pattern:
   rank 0 parses the JSON file and re-serializes to UBJSON; the binary form
   is broadcast; every rank deserializes and re-validates against its own
   in-memory copy.
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
    error->all(FLERR, "Error deserializing broadcast symmetry data: {}", e.what());
  }
  ubjson.clear();

  populate(data);
  compute_inverses();
  validate_group();
  compute_wyckoff_projectors();
}

/* ----------------------------------------------------------------------
   pull ops and orbits out of the JSON tree into the in-memory structs.
   Runs on every rank against an identical deserialized copy, so schema
   errors are raised via error->all().

   Schema:
     {
       "group_name": "<name>",        // optional, default "unknown"
       "lattice":    "<family>",      // required: triclinic / monoclinic /
                                      //   orthorhombic / tetragonal /
                                      //   hexagonal / trigonal / cubic
       "ops": [                       // required: >= 1 entry
         { "R": [[..],[..],[..]],     // 3x3 integer matrix
           "t": [tx,ty,tz] },         // numbers or rational strings "1/2"
         ...
       ],
       "orbits": [                    // required: >= 0 entries
         { "tags": [t0, t1, ..., tN-1] },   // parallel to "ops"; t0 is the
                                            // asymmetric representative
         ...
       ]
     }
------------------------------------------------------------------------- */

void SymmetryGroup::populate(const json &data)
{
  try {
    if (data.contains("group_name")) group_name = data.at("group_name").get<std::string>();
    lattice_family = parse_lattice(data.at("lattice").get<std::string>());

    const auto &jops = data.at("ops");
    if (!jops.is_array() || jops.empty())
      throw std::runtime_error("'ops' must be a non-empty array");
    ops.assign(jops.size(), SymmOp{});
    for (size_t o = 0; o < jops.size(); ++o) {
      const auto &je = jops[o];
      const auto &jR = je.at("R");
      const auto &jt = je.at("t");
      if (!jR.is_array() || jR.size() != 3)
        throw std::runtime_error("op " + std::to_string(o + 1) + ": R must be a 3x3 array");
      if (!jt.is_array() || jt.size() != 3)
        throw std::runtime_error("op " + std::to_string(o + 1) + ": t must be a length-3 array");
      for (int i = 0; i < 3; ++i) {
        if (!jR[i].is_array() || jR[i].size() != 3)
          throw std::runtime_error("op " + std::to_string(o + 1) + ": R row must have 3 entries");
        for (int j = 0; j < 3; ++j) {
          const auto rval = jR[i][j].get<double>();
          if (std::fabs(rval - std::round(rval)) > R_INT_TOL)
            throw std::runtime_error(fmt::format("op {}: R should be an integer but got {:.10g}", o + 1, rval));
          ops[o].R[i][j] = rval;
        }

        ops[o].t[i] = get_t_value(jt[i]);
      }
    }

    const auto &jorb = data.at("orbits");
    if (!jorb.is_array()) throw std::runtime_error("'orbits' must be an array");
    const int N = static_cast<int>(ops.size());
    orbits.assign(jorb.size(), Orbit{});
    for (size_t o = 0; o < jorb.size(); ++o) {
      const auto &jt = jorb[o].at("tags");
      if (!jt.is_array() || static_cast<int>(jt.size()) != N)
        throw std::runtime_error("orbit " + std::to_string(o + 1) + ": 'tags' must have exactly " +
                                 std::to_string(N) + " entries (one per op)");
      orbits[o].image_tag.assign(N, 0);
      for (int k = 0; k < N; ++k) orbits[o].image_tag[k] = jt[k].get<tagint>();
      orbits[o].asym_tag = orbits[o].image_tag[0];

      // optional site-symmetry: list of op indices (1-based in the file) that
      // map the asym atom onto itself. For each such op, tags[op-1] must
      // equal asym_tag -- otherwise the orbit is internally inconsistent.
      orbits[o].stabilizer.clear();
      if (jorb[o].contains("site_symmetry")) {
        const auto &jss = jorb[o].at("site_symmetry");
        if (!jss.is_array())
          throw std::runtime_error("orbit " + std::to_string(o + 1) +
                                   ": 'site_symmetry' must be an array of op indices");
        for (const auto &je : jss) {
          int k1 = je.get<int>();
          if (k1 < 2 || k1 > N)
            throw std::runtime_error("orbit " + std::to_string(o + 1) + ": site_symmetry op " +
                                     std::to_string(k1) + " out of range [2," + std::to_string(N) +
                                     "]");
          int k = k1 - 1;
          if (orbits[o].image_tag[k] != orbits[o].asym_tag)
            throw std::runtime_error("orbit " + std::to_string(o + 1) + ": site_symmetry op " +
                                     std::to_string(k1) + " requires tags[" + std::to_string(k1) +
                                     "] == asym tag, but they differ");
          orbits[o].stabilizer.push_back(k);
        }
      }
    }
  } catch (std::exception &e) {
    error->all(FLERR, Error::NOLASTLINE, "Symmetry file schema error: {}", e.what());
  }
}

/* ----------------------------------------------------------------------
   invert each op.R into op.Rinv via the adjugate / det formula. R is an
   integer-valued lattice operator with |det R| = 1, so this is exact.
------------------------------------------------------------------------- */

void SymmetryGroup::compute_inverses()
{
  for (auto &op : ops) {
    double (&R)[3][3] = op.R;
    double det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
        R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
        R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
    if (std::fabs(std::fabs(det) - 1.0) > R_INT_TOL)
      error->all(FLERR, Error::NOLASTLINE,
                 "Symmetry op R is not a unimodular lattice operator (det = {:.6g})", det);
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
   For each orbit on a special position, precompute the symmetric
   positive-semidefinite matrix B = sum_{k in stabilizer} (R_k - I)^T (R_k - I)
   and its Moore-Penrose pseudo-inverse. B depends only on the operator
   table, so it can be computed once at init -- the per-orbit shift c that
   depends on the asym atom's initial fractional position is filled in
   later by FixSymmetry.

   The pseudo-inverse is constructed via the eigendecomposition of B
   (B = Q Lambda Q^T from MathEigen::jacobi3) by zeroing out eigenvalues
   below a relative threshold:
       Binv = sum_i [ lambda_i > tol ? 1/lambda_i : 0 ] * v_i * v_i^T
   This is the standard SVD pseudo-inverse trick: A^+ = (A^T A)^+ A^T
   where (A^T A)^+ has truncated singular values.
------------------------------------------------------------------------- */

void SymmetryGroup::compute_wyckoff_projectors()
{
  for (auto &orb : orbits) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        orb.B[i][j] = 0.0;
        orb.Binv[i][j] = 0.0;
      }
    if (orb.stabilizer.empty()) continue;

    // B = sum_k (R_k - I)^T (R_k - I)
    for (int k : orb.stabilizer) {
      const auto &op = ops[k];
      double M[3][3];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) M[i][j] = op.R[i][j] - (i == j ? 1.0 : 0.0);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          double s = 0.0;
          for (int m = 0; m < 3; ++m) s += M[m][i] * M[m][j];
          orb.B[i][j] += s;
        }
    }

    // pseudo-inverse via 3x3 symmetric eigendecomposition
    double eval[3], evec[3][3];
    if (MathEigen::jacobi3(orb.B, eval, evec) != 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "SymmetryGroup: Jacobi eigensolver failed on Wyckoff matrix");
    const double max_eval = std::max({eval[0], eval[1], eval[2]});
    const double sv_tol = (max_eval > 0.0) ? max_eval * 1.0e-12 : 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double sum = 0.0;
        for (int m = 0; m < 3; ++m) {
          if (eval[m] > sv_tol) sum += evec[m][i] * evec[m][j] / eval[m];
        }
        orb.Binv[i][j] = sum;
      }
  }
}

/* ----------------------------------------------------------------------
   structural and group-theoretic validation
------------------------------------------------------------------------- */

static double wrap01(double v)
{
  v -= std::floor(v);
  if (v >= 1.0) v -= 1.0;
  return v;
}

static bool R_equal(const double A[3][3], const double B[3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (std::fabs(A[i][j] - B[i][j]) > R_INT_TOL) return false;
  return true;
}

static bool t_equal_mod1(const double a[3], const double b[3])
{
  for (int k = 0; k < 3; ++k) {
    double d = wrap01(a[k] - b[k]);
    if (d > 0.5) d = 1.0 - d;
    if (d > T_MOD1_TOL) return false;
  }
  return true;
}

void SymmetryGroup::validate_group()
{
  if (ops.empty()) error->all(FLERR, Error::NOLASTLINE, "Symmetry file declared zero operators");

  // op 0 must be the identity {R = I, t = 0}
  const double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const double zero[3] = {0, 0, 0};
  if (!R_equal(ops[0].R, I) || !t_equal_mod1(ops[0].t, zero))
    error->all(FLERR, Error::NOLASTLINE,
               "Symmetry file: first op must be the identity (R = I, t = 0)");

  // closure: for every a, b find c such that R_a * R_b == R_c and
  // t_a + R_a * t_b == t_c (mod lattice translations).
  const int N = static_cast<int>(ops.size());
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
      int found = -1;
      for (int c = 0; c < N; ++c) {
        if (R_equal(Rab, ops[c].R) && t_equal_mod1(tab, ops[c].t)) {
          found = c;
          break;
        }
      }
      if (found < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Symmetry group is not closed: op {} * op {} is not in the operator list", a + 1,
                   b + 1);
    }
  }

  // orbit structural sanity
  for (size_t o = 0; o < orbits.size(); ++o) {
    const Orbit &orb = orbits[o];
    if (static_cast<int>(orb.image_tag.size()) != N)
      error->all(FLERR, Error::NOLASTLINE, "Orbit {} has {} tag entries, expected {}", o + 1,
                 orb.image_tag.size(), N);
    for (int k = 0; k < N; ++k)
      if (orb.image_tag[k] <= 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Orbit {} has invalid (non-positive) tag {} at position {}", o + 1,
                   orb.image_tag[k], k + 1);
  }
}
