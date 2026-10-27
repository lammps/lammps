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
   find-symmetry: write symmetry files for fix symmetry with spglib

   find-symmetry [options] <data file>
       determine the space group of the atoms in a LAMMPS data file and
       write its operators in fractional coordinates of the simulation box

   find-symmetry [options] -g <space group number>
       write the operators of a space group in its standard setting
------------------------------------------------------------------------- */

#include "spglib.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

// tolerance for comparing translations, in fractional coordinates
constexpr double TTOL = 1.0e-5;

struct Op {
  int R[3][3];
  double t[3];
};

[[noreturn]] void fail(const std::string &mesg)
{
  std::cerr << "find-symmetry: " << mesg << "\n";
  exit(1);
}

void usage()
{
  std::cout
      << "Usage: find-symmetry [options] <data file>\n"
         "       find-symmetry [options] -g <space group number>\n\n"
         "Write a symmetry file for the LAMMPS fix symmetry command.\n\n"
         "With a LAMMPS data file, the space group of the atoms in the simulation box is\n"
         "determined with spglib and its operators are written in fractional coordinates\n"
         "of the box.  If the box consists of several unit cells, the operators of one\n"
         "unit cell are written together with a \"supercell\" entry.\n\n"
         "With -g, the operators of a space group in its standard setting are written.\n"
         "Then the simulation box must be the conventional unit cell of that setting.\n\n"
         "Options:\n"
         "  -g <number>       space group number from 1 to 230\n"
         "  -H <number>       Hall number from 1 to 530 to select a different setting\n"
         "  -s <nx> <ny> <nz> with -g or -H: add a supercell entry\n"
         "  -t <distance>     tolerance for finding the symmetry in distance units of\n"
         "                    the data file (default 1.0e-5)\n"
         "  -n                do not reduce the operators to a unit cell; this can\n"
         "                    produce a very large number of operators\n"
         "  -o <file>         write to file instead of the screen\n"
         "  -h                print this help\n";
}

double wrap01(double v)
{
  v -= std::floor(v);
  if (v >= 1.0) v -= 1.0;
  return v;
}

// true if a - b is an integer within TTOL
bool same_mod1(double a, double b)
{
  const double d = a - b;
  return std::fabs(d - std::round(d)) < TTOL;
}

bool same_op(const Op &a, const Op &b)
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      if (a.R[i][j] != b.R[i][j]) return false;
    if (!same_mod1(a.t[i], b.t[i])) return false;
  }
  return true;
}

bool is_unit(const int R[3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (R[i][j] != ((i == j) ? 1 : 0)) return false;
  return true;
}

// fraction p/q with the smallest q <= 48 for a translation, or a decimal number
std::string format_translation(double t)
{
  t = wrap01(t);
  for (int q = 1; q <= 48; ++q) {
    const double p = std::round(t * q);
    if (std::fabs(t - p / q) < TTOL) {
      if (p == 0.0 || p == q) return "";
      return "+" + std::to_string(static_cast<int>(p)) + "/" + std::to_string(q);
    }
  }
  char buf[32];
  snprintf(buf, sizeof(buf), "+%.10g", t);
  return buf;
}

// operator in coordinate triplet notation, e.g. "-y,x-y,z+1/3"
std::string format_op(const Op &op)
{
  std::string result;
  const char xyz[3] = {'x', 'y', 'z'};
  for (int i = 0; i < 3; ++i) {
    std::string comp;
    for (int j = 0; j < 3; ++j) {
      const int c = op.R[i][j];
      if (c == 0) continue;
      if (c < 0)
        comp += "-";
      else if (!comp.empty())
        comp += "+";
      if (std::abs(c) != 1) comp += std::to_string(std::abs(c));
      comp += xyz[j];
    }
    comp += format_translation(op.t[i]);
    if (i > 0) result += ",";
    result += comp;
  }
  return result;
}

/* ----------------------------------------------------------------------
   read the box and the atom positions and types from a LAMMPS data file
------------------------------------------------------------------------- */

struct Structure {
  double lattice[3][3];    // columns are the box vectors, as used by spglib
  std::vector<std::array<double, 3>> frac;
  std::vector<int> types;
};

// columns of the atom type and the x coordinate in the Atoms section
const std::map<std::string, std::pair<int, int>> style_columns = {
    {"angle", {2, 3}},     {"atomic", {1, 2}}, {"body", {1, 4}},       {"bond", {2, 3}},
    {"bpm/sphere", {2, 5}}, {"charge", {1, 3}}, {"dielectric", {1, 3}}, {"dipole", {1, 3}},
    {"dpd", {1, 3}},       {"electron", {1, 5}}, {"ellipsoid", {1, 4}}, {"full", {2, 4}},
    {"hybrid", {1, 2}},    {"molecular", {2, 3}}, {"peri", {1, 4}},     {"sphere", {1, 4}},
    {"spin", {1, 2}},      {"template", {1, 5}}};

std::vector<std::string> split_words(const std::string &line)
{
  std::vector<std::string> words;
  std::istringstream in(line.substr(0, line.find('#')));
  std::string w;
  while (in >> w) words.push_back(w);
  return words;
}

Structure read_data(const std::string &file)
{
  std::ifstream in(file);
  if (!in) fail("cannot open data file " + file);

  long natoms = -1;
  double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {0.0, 0.0, 0.0};
  double xy = 0.0, xz = 0.0, yz = 0.0;
  bool have_box[3] = {false, false, false};
  std::string line;
  std::getline(in, line);    // title line

  // header

  std::string style;
  bool found_atoms = false;
  while (std::getline(in, line)) {
    auto words = split_words(line);
    if (words.empty()) continue;
    if (words[0] == "Atoms") {
      found_atoms = true;
      auto hash = line.find('#');
      if (hash != std::string::npos) {
        std::istringstream sin(line.substr(hash + 1));
        sin >> style;
      }
      break;
    }
    const auto n = words.size();
    if ((n == 2) && (words[1] == "atoms")) {
      natoms = std::stol(words[0]);
    } else if ((n == 4) && (words[2] == "xlo") && (words[3] == "xhi")) {
      lo[0] = std::stod(words[0]);
      hi[0] = std::stod(words[1]);
      have_box[0] = true;
    } else if ((n == 4) && (words[2] == "ylo") && (words[3] == "yhi")) {
      lo[1] = std::stod(words[0]);
      hi[1] = std::stod(words[1]);
      have_box[1] = true;
    } else if ((n == 4) && (words[2] == "zlo") && (words[3] == "zhi")) {
      lo[2] = std::stod(words[0]);
      hi[2] = std::stod(words[1]);
      have_box[2] = true;
    } else if ((n == 6) && (words[3] == "xy") && (words[4] == "xz") && (words[5] == "yz")) {
      xy = std::stod(words[0]);
      xz = std::stod(words[1]);
      yz = std::stod(words[2]);
    } else if ((n > 3) && ((words[3] == "avec") || (words[3] == "bvec") || (words[3] == "cvec"))) {
      fail("data files for general triclinic boxes are not supported, use write_data without "
           "the triclinic/general keyword");
    }
  }
  if (!found_atoms) fail("no Atoms section in data file " + file);
  if (natoms <= 0) fail("no atoms in data file " + file);
  if (!have_box[0] || !have_box[1] || !have_box[2])
    fail("box dimensions are missing in data file " + file);

  int tcol = 1, xcol = 2;
  if (!style.empty()) {
    auto it = style_columns.find(style);
    if (it == style_columns.end()) {
      std::string known;
      for (const auto &kv : style_columns) known += " " + kv.first;
      fail("atom style " + style + " in data file is not supported, known styles:" + known);
    }
    tcol = it->second.first;
    xcol = it->second.second;
  }

  Structure st;
  const double prd[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
  const double lattice[3][3] = {{prd[0], xy, xz}, {0.0, prd[1], yz}, {0.0, 0.0, prd[2]}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) st.lattice[i][j] = lattice[i][j];

  // atoms; type labels are numbered in the order of their appearance

  std::map<std::string, int> labels;
  while ((static_cast<long>(st.types.size()) < natoms) && std::getline(in, line)) {
    auto words = split_words(line);
    if (words.empty()) continue;
    if (static_cast<int>(words.size()) < xcol + 3)
      fail("incorrect format of line in Atoms section: " + line);
    int type;
    const auto &tword = words[tcol];
    if (tword.find_first_not_of("0123456789") == std::string::npos) {
      type = std::stoi(tword);
    } else {
      auto it = labels.find(tword);
      if (it == labels.end()) it = labels.emplace(tword, -1 - static_cast<int>(labels.size())).first;
      type = it->second;
    }
    const double x = std::stod(words[xcol]) - lo[0];
    const double y = std::stod(words[xcol + 1]) - lo[1];
    const double z = std::stod(words[xcol + 2]) - lo[2];
    std::array<double, 3> s;
    s[2] = z / prd[2];
    s[1] = (y - yz * s[2]) / prd[1];
    s[0] = (x - xy * s[1] - xz * s[2]) / prd[0];
    for (auto &sval : s) sval = wrap01(sval);
    st.frac.push_back(s);
    st.types.push_back(type);
  }
  if (static_cast<long>(st.types.size()) != natoms)
    fail("found only " + std::to_string(st.types.size()) + " of " + std::to_string(natoms) +
         " atoms in data file " + file);
  return st;
}

/* ----------------------------------------------------------------------
   If the box consists of nx x ny x nz unit cells, i.e. the translations by
   1/nx, 1/ny, 1/nz along the box vectors are symmetries, convert the
   operators to fractional coordinates of the unit cell and keep one of each
   set that differs only by such a translation.  Returns false if the
   operators cannot be written that way.
------------------------------------------------------------------------- */

bool reduce_to_unitcell(const std::vector<Op> &ops, std::vector<Op> &cellops, int ncell[3])
{
  for (int d = 0; d < 3; ++d) {
    ncell[d] = 0;
    for (const auto &op : ops) {
      if (!is_unit(op.R)) continue;
      bool along = true;
      for (int e = 0; e < 3; ++e)
        if ((e != d) && !same_mod1(op.t[e], 0.0)) along = false;
      if (along) ++ncell[d];
    }
  }
  const int nprod = ncell[0] * ncell[1] * ncell[2];
  if ((nprod < 1) || (ops.size() % nprod != 0)) return false;

  cellops.clear();
  for (const auto &op : ops) {
    Op cop;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const int num = op.R[i][j] * ncell[j];
        if (num % ncell[i] != 0) return false;
        cop.R[i][j] = num / ncell[i];
      }
      cop.t[i] = wrap01(op.t[i] * ncell[i]);
    }
    bool known = false;
    for (const auto &other : cellops)
      if (same_op(cop, other)) known = true;
    if (!known) cellops.push_back(cop);
  }
  return cellops.size() * nprod == ops.size();
}

/* ----------------------------------------------------------------------
   largest deviation of the atom positions from exact symmetry: distance
   in fractional coordinates between the image of each atom under each
   operator and the closest atom of the same type.  Atoms are sorted into
   a grid of cells no smaller than the search range.
------------------------------------------------------------------------- */

double max_deviation(const Structure &st, const std::vector<Op> &ops, double range)
{
  const int ncell = static_cast<int>(std::max(1.0, std::min(1024.0, std::floor(0.5 / range))));
  auto cell_of = [ncell](double s) { return std::min(static_cast<int>(s * ncell), ncell - 1); };
  auto key = [ncell](int a, int b, int c) {
    return (static_cast<long long>(a) * ncell + b) * ncell + c;
  };
  std::map<long long, std::vector<int>> cells;
  const int natoms = static_cast<int>(st.types.size());
  for (int i = 0; i < natoms; ++i)
    cells[key(cell_of(st.frac[i][0]), cell_of(st.frac[i][1]), cell_of(st.frac[i][2]))].push_back(i);

  std::vector<int> offsets = {0};
  if (ncell > 1) offsets.push_back(1);
  if (ncell > 2) offsets.push_back(-1);

  double maxdev = 0.0;
  for (const auto &op : ops) {
    for (int i = 0; i < natoms; ++i) {
      double s[3];
      for (int d = 0; d < 3; ++d) {
        s[d] = op.t[d];
        for (int e = 0; e < 3; ++e) s[d] += op.R[d][e] * st.frac[i][e];
        s[d] = wrap01(s[d]);
      }
      double best = 1.0;
      const int c[3] = {cell_of(s[0]), cell_of(s[1]), cell_of(s[2])};
      for (int da : offsets)
        for (int db : offsets)
          for (int dc : offsets) {
            auto it = cells.find(key((c[0] + da + ncell) % ncell, (c[1] + db + ncell) % ncell,
                                     (c[2] + dc + ncell) % ncell));
            if (it == cells.end()) continue;
            for (int j : it->second) {
              if (st.types[j] != st.types[i]) continue;
              double dev = 0.0;
              for (int d = 0; d < 3; ++d) {
                const double delta = st.frac[j][d] - s[d];
                dev = std::max(dev, std::fabs(delta - std::round(delta)));
              }
              best = std::min(best, dev);
            }
          }
      maxdev = std::max(maxdev, best);
    }
  }
  return maxdev;
}

// number of the default setting (Hall number) of a space group
int default_hall_number(int number)
{
  for (int hall = 1; hall <= 530; ++hall)
    if (spg_get_spacegroup_type(hall).number == number) return hall;
  fail("space group number must be between 1 and 230");
}

// shape of the simulation box required by the standard setting of a space group
std::string box_shape(int number)
{
  if (number <= 2) return "any shape";
  if (number <= 15) return "triclinic with xy = yz = 0 (b is the unique axis)";
  if (number <= 74) return "orthogonal";
  if (number <= 142) return "orthogonal with equal lengths along x and y";
  if (number <= 194)
    return "hexagonal with equal lengths a and b at 120 degrees: xy = -a/2, ly = a*sqrt(3)/2";
  return "cubic";
}
}    // namespace

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  std::string datafile, outfile;
  int number = 0, hall = 0;
  int supercell[3] = {1, 1, 1};
  double symprec = 1.0e-5;
  bool reduce = true;

  for (int iarg = 1; iarg < argc; ++iarg) {
    const std::string arg = argv[iarg];
    auto next = [&]() -> std::string {
      if (iarg + 1 >= argc) fail("missing argument for option " + arg);
      return argv[++iarg];
    };
    try {
      if ((arg == "-h") || (arg == "--help")) {
        usage();
        return 0;
      } else if (arg == "-g") {
        number = std::stoi(next());
      } else if (arg == "-H") {
        hall = std::stoi(next());
      } else if (arg == "-s") {
        for (int &n : supercell) n = std::stoi(next());
      } else if (arg == "-t") {
        symprec = std::stod(next());
      } else if (arg == "-n") {
        reduce = false;
      } else if (arg == "-o") {
        outfile = next();
      } else if ((arg.size() > 1) && (arg[0] == '-')) {
        fail("unknown option " + arg + ", use -h for help");
      } else {
        if (!datafile.empty()) fail("more than one data file given");
        datafile = arg;
      }
    } catch (std::logic_error &) {
      fail("invalid argument for option " + arg);
    }
  }

  std::vector<Op> ops;
  std::string name;
  int ncell[3] = {1, 1, 1};

  if ((number > 0) || (hall > 0)) {
    // operators from the spglib database
    if (!datafile.empty()) fail("a data file cannot be combined with -g or -H");
    if (hall == 0) hall = default_hall_number(number);
    if ((hall < 1) || (hall > 530)) fail("Hall number must be between 1 and 530");
    const auto type = spg_get_spacegroup_type(hall);
    int rotations[192][3][3];
    double translations[192][3];
    const int nops = spg_get_symmetry_from_database(rotations, translations, hall);
    for (int k = 0; k < nops; ++k) {
      Op op;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) op.R[i][j] = rotations[k][i][j];
        op.t[i] = translations[k][i];
      }
      ops.push_back(op);
    }
    name = std::string(type.international_short) + " (" + std::to_string(type.number) + ")";
    if (type.choice[0] != '\0') name += " setting " + std::string(type.choice);
    for (int d = 0; d < 3; ++d) {
      if (supercell[d] < 1) fail("supercell sizes must be positive");
      ncell[d] = supercell[d];
    }
    std::cerr << "Space group " << name << " with " << ops.size()
              << " operators in the conventional unit cell\n"
              << "Shape of the simulation box for one unit cell: " << box_shape(type.number) << "\n";
  } else {
    // operators from the positions in a data file
    if (datafile.empty()) {
      usage();
      return 1;
    }
    auto st = read_data(datafile);
    const int natoms = static_cast<int>(st.types.size());
    std::vector<double> pos(3 * natoms);
    for (int i = 0; i < natoms; ++i)
      for (int d = 0; d < 3; ++d) pos[3 * i + d] = st.frac[i][d];

    SpglibDataset *data =
        spg_get_dataset(st.lattice, reinterpret_cast<double(*)[3]>(pos.data()), st.types.data(),
                        natoms, symprec);
    if (!data) fail(std::string("spglib failed: ") + spg_get_error_message(spg_get_error_code()));
    for (int k = 0; k < data->n_operations; ++k) {
      Op op;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) op.R[i][j] = data->rotations[k][i][j];
        op.t[i] = data->translations[k][i];
      }
      ops.push_back(op);
    }
    name = std::string(data->international_symbol) + " (" + std::to_string(data->spacegroup_number) +
        ")";
    std::cerr << "Space group " << name << " with " << ops.size()
              << " operators in the simulation box\n";
    spg_free_dataset(data);

    // fix symmetry needs a tolerance larger than the deviation from exact symmetry
    const double minlen = std::min({st.lattice[0][0], st.lattice[1][1], st.lattice[2][2]});
    const double maxdev = max_deviation(st, ops, 2.0 * symprec / minlen);
    if (maxdev > 0.5e-6) {
      char buf[64];
      snprintf(buf, sizeof(buf), "%.1e", 2.0 * maxdev);
      std::cerr << "The positions deviate from exact symmetry by up to " << maxdev
                << " in fractional coordinates;\nuse fix symmetry keyword tol " << buf << "\n";
    }

    std::vector<Op> cellops;
    if (reduce && reduce_to_unitcell(ops, cellops, ncell)) {
      if (ncell[0] * ncell[1] * ncell[2] > 1) {
        ops = cellops;
        std::cerr << "The box consists of " << ncell[0] << " x " << ncell[1] << " x " << ncell[2]
                  << " unit cells with " << ops.size() << " operators each\n";
      }
    } else {
      ncell[0] = ncell[1] = ncell[2] = 1;
    }
  }

  // the identity first
  for (std::size_t k = 0; k < ops.size(); ++k) {
    if (is_unit(ops[k].R) && same_mod1(ops[k].t[0], 0.0) && same_mod1(ops[k].t[1], 0.0) &&
        same_mod1(ops[k].t[2], 0.0)) {
      std::swap(ops[0], ops[k]);
      break;
    }
  }

  std::ostringstream out;
  out << "{\n  \"group_name\": \"" << name << "\",\n  \"ops\": [\n";
  for (std::size_t k = 0; k < ops.size(); ++k)
    out << "    \"" << format_op(ops[k]) << "\"" << ((k + 1 < ops.size()) ? ",\n" : "\n");
  out << "  ]";
  if (ncell[0] * ncell[1] * ncell[2] > 1)
    out << ",\n  \"supercell\": [" << ncell[0] << ", " << ncell[1] << ", " << ncell[2] << "]";
  out << "\n}\n";

  if (outfile.empty()) {
    std::cout << out.str();
  } else {
    std::ofstream fout(outfile);
    if (!fout) fail("cannot open output file " + outfile);
    fout << out.str();
  }
  return 0;
}
