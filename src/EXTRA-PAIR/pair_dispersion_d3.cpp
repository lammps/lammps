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
   Contributing authors:
      Sonia Salomoni^1,^2
      Arthur France-Lanord^1

      ^1: IMPMC, CNRS, Sorbonne Universite, Paris, France
      ^2: SCAI, Sorbonne Universite, Paris, France
------------------------------------------------------------------------- */

#include "pair_dispersion_d3.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;

// global ad hoc parameters
static constexpr double K1 = 16.0;
static constexpr double K3 = -4.0;

/*  reasonable choices for k3 are between 3 and 5 :
    this gives smoth curves with maxima around the integer values
    k3=3 give for CN=0 a slightly smaller value than computed
    for the free atom. This also yields to larger CN for atoms
    in larger molecules but with the same chemical environment
    which is physically not right.
    values >5 might lead to bumps in the potential.
*/

static constexpr int NUM_ELEMENTS = 94;      // maximum element number
static constexpr int N_PARS_COLS = 5;        // number of columns in C6 table
static constexpr int N_PARS_ROWS = 32385;    // number of rows in C6 table

static constexpr double autoang = 0.52917725;    // atomic units (Bohr) to Angstrom
static constexpr double autoev = 27.21140795;    // atomic units (Hartree) to eV

#include "d3_parameters.h"
/* ----------------------------------------------------------------------
   Constructor (Required)
------------------------------------------------------------------------- */

PairDispersionD3::PairDispersionD3(LAMMPS *lmp) :
    Pair(lmp), r2r4(nullptr), rcov(nullptr), mxci(nullptr), r0ab(nullptr), c6ab(nullptr),
    cn(nullptr), dc6(nullptr)
{
  nmax = 0;
  comm_forward = 2;
  comm_reverse = 2;

  restartinfo = 0;
  manybody_flag = 1;
  one_coeff = 1;
  single_enable = 0;

  s6 = s8 = s18 = rs6 = rs8 = rs18 = a1 = a2 = alpha = alpha6 = alpha8 = 0.0;
}

/* ----------------------------------------------------------------------
   Destructor (Required)
------------------------------------------------------------------------- */

PairDispersionD3::~PairDispersionD3()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(r0ab);
    memory->destroy(c6ab);
    memory->destroy(mxci);

    memory->destroy(r2r4);
    memory->destroy(rcov);

    memory->destroy(cn);
    memory->destroy(dc6);
  }
}

/* ----------------------------------------------------------------------
   Allocate all arrays (Required)
------------------------------------------------------------------------- */

void PairDispersionD3::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(mxci, n + 1, "pair:mxci");
  memory->create(r2r4, n + 1, "pair:r2r4");
  memory->create(rcov, n + 1, "pair:rcov");

  memory->create(r0ab, n + 1, n + 1, "pair:r0ab");
  memory->create(c6ab, n + 1, n + 1, 5, 5, 3, "pair:c6ab");
}

/* ----------------------------------------------------------------------
   Settings: read from pair_style (Required)
             pair_style   d3 rthr cn_thr damping_type
------------------------------------------------------------------------- */

void PairDispersionD3::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Pair_style dispersion/d3 needs 4 arguments");

  damping_type = arg[0];
  std::string functional_name = arg[1];

  std::transform(damping_type.begin(), damping_type.end(), damping_type.begin(), ::tolower);

  rthr = utils::numeric(FLERR, arg[2], false, lmp);
  cn_thr = utils::numeric(FLERR, arg[3], false, lmp);

  rthr = rthr * rthr;
  cn_thr = cn_thr * cn_thr;    // squared cutoff

  set_funcpar(functional_name);
}

/* ----------------------------------------------------------------------
   finds atomic number (used in PairDispersionD3::coeff)
------------------------------------------------------------------------- */

int PairDispersionD3::find_atomic_number(std::string &key)
{
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);
  if (key.length() == 1) key += " ";
  key.resize(2);

  std::vector<std::string> element_table = {
      "h ", "he", "li", "be", "b ", "c ", "n ", "o ", "f ", "ne", "na", "mg", "al", "si",
      "p ", "s ", "cl", "ar", "k ", "ca", "sc", "ti", "v ", "cr", "mn", "fe", "co", "ni",
      "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y ", "zr", "nb", "mo",
      "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i ", "xe", "cs", "ba",
      "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb",
      "lu", "hf", "ta", "w ", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po",
      "at", "rn", "fr", "ra", "ac", "th", "pa", "u ", "np", "pu"};

  for (size_t i = 0; i < element_table.size(); ++i) {
    if (element_table[i] == key) {
      int atomic_number = i + 1;
      return atomic_number;
    }
  }

  // if not the case
  return -1;
}

/* ----------------------------------------------------------------------
   Check whether an integer value in an integer array (used in PairDispersionD3::coeff)
------------------------------------------------------------------------- */

std::vector<int> PairDispersionD3::is_int_in_array(int array[], int size, int value)
{
  std::vector<int> indices;
  for (int i = 0; i < size; ++i) {
    if (array[i] == value) indices.push_back(i + 1);
  }
  return indices;
}

/* ----------------------------------------------------------------------
   Read r0ab values from r0ab.csv (used in PairDispersionD3::coeff)
------------------------------------------------------------------------- */

void PairDispersionD3::read_r0ab(int *atomic_numbers, int ntypes)
{
  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= ntypes; j++) {
      r0ab[i][j] = r0ab_table[atomic_numbers[i - 1] - 1][atomic_numbers[j - 1] - 1];
    }
  }
}

/* ----------------------------------------------------------------------
   Get atom pair indices and grid indices (used in PairDispersionD3::read_c6ab)
------------------------------------------------------------------------- */

void PairDispersionD3::set_limit_in_pars_array(int &idx_atom_1, int &idx_atom_2, int &idx_i,
                                               int &idx_j)
{
  idx_i = 0;
  idx_j = 0;
  int shift = 100;

  while (idx_atom_1 > shift) {
    idx_atom_1 -= shift;
    idx_i++;
  }

  while (idx_atom_2 > shift) {
    idx_atom_2 -= shift;
    idx_j++;
  }
}

/* ----------------------------------------------------------------------
   Read c6ab values from c6ab.csv (used in PairDispersionD3::coeff)
------------------------------------------------------------------------- */

void PairDispersionD3::read_c6ab(int *atomic_numbers, int ntypes)
{
  for (int i = 0; i <= ntypes; i++) mxci[i] = 0;
  int grid_i = 0, grid_j = 0;

  for (int i = 0; i < N_PARS_ROWS; i++) {
    const double ref_c6 = c6ab_table[i][0];

    int atom_number_1 = std::round(c6ab_table[i][1]);
    int atom_number_2 = std::round(c6ab_table[i][2]);

    set_limit_in_pars_array(atom_number_1, atom_number_2, grid_i, grid_j);

    std::vector<int> idx_atoms_1 = is_int_in_array(atomic_numbers, ntypes, atom_number_1);
    if (idx_atoms_1.empty()) continue;

    std::vector<int> idx_atoms_2 = is_int_in_array(atomic_numbers, ntypes, atom_number_2);
    if (idx_atoms_2.empty()) continue;

    const double ref_cn1 = c6ab_table[i][3];
    const double ref_cn2 = c6ab_table[i][4];

    for (int idx_atom_1 : idx_atoms_1) {
      for (int idx_atom_2 : idx_atoms_2) {
        mxci[idx_atom_1] = std::max(mxci[idx_atom_1], grid_i);
        mxci[idx_atom_2] = std::max(mxci[idx_atom_2], grid_j);

        c6ab[idx_atom_1][idx_atom_2][grid_i][grid_j][0] = ref_c6;
        c6ab[idx_atom_1][idx_atom_2][grid_i][grid_j][1] = ref_cn1;
        c6ab[idx_atom_1][idx_atom_2][grid_i][grid_j][2] = ref_cn2;
        c6ab[idx_atom_2][idx_atom_1][grid_j][grid_i][0] = ref_c6;
        c6ab[idx_atom_2][idx_atom_1][grid_j][grid_i][1] = ref_cn2;
        c6ab[idx_atom_2][idx_atom_1][grid_j][grid_i][2] = ref_cn1;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Coeff: read from pair_coeff (Required)
          pair_coeff * * path_r0ab.csv path_c6ab.csv functional element1 element2 ...
------------------------------------------------------------------------- */

void PairDispersionD3::coeff(int narg, char **arg)
{
  int ntypes = atom->ntypes;
  if (narg != ntypes + 2) error->all(FLERR, "Pair_coeff * * needs: element1 element2 ...");

  if (!allocated) allocate();

  std::string element;
  int *atomic_numbers = (int *) malloc(sizeof(int) * ntypes);
  for (int i = 0; i < ntypes; i++) {
    element = arg[i + 2];
    atomic_numbers[i] = find_atomic_number(element);
  }

  int count = 0;
  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= ntypes; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

  for (int i = 1; i <= ntypes; i++) {
    r2r4[i] = r2r4_ref[atomic_numbers[i - 1]];
    rcov[i] = rcov_ref[atomic_numbers[i - 1]];
  }

  // set r0ab
  read_r0ab(atomic_numbers, ntypes);

  // read c6ab
  read_c6ab(atomic_numbers, ntypes);

  free(atomic_numbers);
}

/* ----------------------------------------------------------------------
   Calculate coordination number of atoms
------------------------------------------------------------------------- */

void PairDispersionD3::calc_coordination_number()
{

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(cn, nmax, "pair:cn");
    memory->grow(dc6, nmax, "pair:dc6");
  }

  // zero out coordination number
  memset(cn, 0, sizeof(double) * (newton_pair ? nall : nlocal));
  memset(dc6, 0, sizeof(double) * (newton_pair ? nall : nlocal));

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {

    int i = ilist[ii];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {

      int j = jlist[jj];
      j &= NEIGHMASK;
      int jtype = type[j];

      double delrj[3];
      delrj[0] = x[i][0] - x[j][0];
      delrj[1] = x[i][1] - x[j][1];
      delrj[2] = x[i][2] - x[j][2];

      double rsq = delrj[0] * delrj[0] + delrj[1] * delrj[1] + delrj[2] * delrj[2];

      // if the atoms are too far away don't consider the contribution
      if (rsq > cn_thr) continue;

      double rr = sqrt(rsq);
      double rcov_ij = (rcov[itype] + rcov[jtype]) * autoang;
      double cn_ij = 1.0f / (1.0f + expf(-K1 * ((rcov_ij / rr) - 1.0f)));

      // update coordination number
      cn[i] += cn_ij;
      if (newton_pair || j < nlocal) { cn[j] += cn_ij; }
    }
  }

  // communicate coordination number
  communicationStage = 1;
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
   Get derivative of C6
------------------------------------------------------------------------- */

double *PairDispersionD3::get_dC6(int iat, int jat, double cni, double cnj)
{

  static double c6_res[3] = {};
  double c6_ref, cni_ref, cnj_ref;
  double c6mem, r_save, r;
  double expterm, term;
  double num, den, d_num_i, d_num_j, d_den_i, d_den_j;

  c6mem = -1.0e20f, r_save = 1.0e20f;
  num = 0;
  den = 0;
  d_num_i = 0;
  d_num_j = 0;
  d_den_i = 0;
  d_den_j = 0;

  for (int ci = 0; ci <= mxci[iat]; ci++) {
    for (int cj = 0; cj <= mxci[jat]; cj++) {

      c6_ref = c6ab[iat][jat][ci][cj][0];
      c6_ref *= autoev * pow(autoang, 6);

      if (c6_ref > 0) {
        cni_ref = c6ab[iat][jat][ci][cj][1];
        cnj_ref = c6ab[iat][jat][ci][cj][2];

        r = (cni - cni_ref) * (cni - cni_ref) + (cnj - cnj_ref) * (cnj - cnj_ref);

        if (r < r_save) {
          r_save = r;
          c6mem = c6_ref;
        }

        expterm = exp(static_cast<double>(K3) * static_cast<double>(r));

        num += c6_ref * expterm;
        den += expterm;

        expterm = expterm * 2.0 * K3;

        term = expterm * (cni - cni_ref);
        d_num_i += c6_ref * term;
        d_den_i += term;

        term = expterm * (cnj - cnj_ref);
        d_num_j += c6_ref * term;
        d_den_j += term;
      }
    }
  }

  if (den > 1.0E-99) {
    c6_res[0] = num / den;
    c6_res[1] = ((d_num_i * den) - (d_den_i * num)) / (den * den);
    c6_res[2] = ((d_num_j * den) - (d_den_j * num)) / (den * den);
  } else {
    c6_res[0] = c6mem;
    c6_res[1] = 0;
    c6_res[2] = 0;
  }
  return c6_res;
}

/* ----------------------------------------------------------------------
   Compute : energy, force, and stress (Required)
------------------------------------------------------------------------- */

void PairDispersionD3::compute(int eflag, int vflag)
{

  std::unordered_map<std::string, int> dampingMap = {
      {"zero", 1}, {"zerom", 2}, {"bj", 3}, {"bjm", 4}};
  int dampingCode = dampingMap[damping_type];

  double evdwl = 0.0;
  ev_init(eflag, vflag);

  calc_coordination_number();

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // start computing quantities
  for (int ii = 0; ii < inum; ii++) {

    int i = ilist[ii];

    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    int jnum = numneigh[i];
    int *jlist = firstneigh[i];

    for (int jj = 0; jj < jnum; jj++) {

      int j = jlist[jj];
      double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];

      double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[type[i]][type[j]]) {

        double r = sqrt(rsq);
        double r2inv = 1.0f / rsq;
        double r6inv = r2inv * r2inv * r2inv;
        double r8inv = r2inv * r2inv * r2inv * r2inv;
        double r10inv = r2inv * r2inv * r2inv * r2inv * r2inv;

        double *c6_res = get_dC6(type[i], type[j], cn[i], cn[j]);

        double C6 = c6_res[0];
        double C8 = 3.0 * C6 * r2r4[type[i]] * r2r4[type[j]] * autoang * autoang;

        double alpha6 = alpha;
        double alpha8 = alpha + 2;

        double t6, t8, damp6, damp8, e6, e8;
        double tmp6, tmp8, fpair1, fpair2, fpair;
        t6 = t8 = e6 = e8 = evdwl = fpair = fpair1 = fpair2 = 0.0;

        switch (dampingCode) {
          case 1: {    // zero

            double r0 = r / r0ab[type[i]][type[j]];

            t6 = pow(rs6 / r0, alpha6);
            damp6 = 1.0f / (1.0f + 6.0f * t6);
            t8 = pow(rs8 / r0, alpha8);
            damp8 = 1.0f / (1.0f + 6.0f * t8);

            e6 = C6 * damp6 * r6inv;
            e8 = C8 * damp8 * r8inv;

            tmp6 = 6 * s6 * C6 * r8inv * damp6;
            tmp8 = 8 * s8 * C8 * r10inv * damp8;

            fpair1 = -tmp6 - tmp8;
            fpair2 = tmp6 * alpha6 * t6 * damp6 + (3.0f / 4) * tmp8 * alpha8 * t8 * damp8;

            fpair = fpair1 + fpair2;
            fpair *= factor_lj;
          } break;
          case 2: {    // zerom

            double r0 = r0ab[type[i]][type[j]];

            t6 = pow((r / (rs6 * r0)) + rs8 * r0, -alpha6);
            damp6 = 1.0f / (1.0f + 6.0f * t6);
            t8 = pow((r / r0) + rs8 * r0, -alpha8);
            damp8 = 1.0f / (1.0f + 6.0f * t8);

            e6 = C6 * damp6 * r6inv;
            e8 = C8 * damp8 * r8inv;

            tmp6 = 6 * s6 * C6 * r8inv * damp6;
            tmp8 = 8 * s8 * C8 * r10inv * damp8;

            fpair1 = -tmp6 - tmp8;

            double fp26 = tmp6 * alpha6 * t6 * damp6 * r / (r + rs6 * rs8 * r0 * r0);
            double fp28 = tmp8 * alpha8 * t8 * damp8 * r / (r + rs8 * r0 * r0);

            fpair2 = fp26 + (3.0f / 4) * fp28;

            fpair = fpair1 + fpair2;
            fpair *= factor_lj;
          } break;
          case 3: {    // bj

            double r0 = sqrt(C8 / C6);

            double r4 = rsq * rsq;
            double r6 = rsq * rsq * rsq;
            double r8 = rsq * rsq * rsq * rsq;

            t6 = r6 + pow((a1 * r0 + a2), 6);
            t8 = r8 + pow((a1 * r0 + a2), 8);

            e6 = C6 / t6;
            e8 = C8 / t8;

            tmp6 = 6.0 * s6 * C6 * r4 / (t6 * t6);
            tmp8 = 8.0 * s8 * C8 * r6 / (t8 * t8);

            fpair = -(tmp6 + tmp8);
            fpair *= factor_lj;
          } break;
          case 4: {    // bjm

            double r0 = sqrt(C8 / C6);

            double r4 = rsq * rsq;
            double r6 = rsq * rsq * rsq;
            double r8 = rsq * rsq * rsq * rsq;

            t6 = r6 + pow((a1 * r0 + a2), 6);
            t8 = r8 + pow((a1 * r0 + a2), 8);

            e6 = C6 / t6;
            e8 = C8 / t8;

            tmp6 = 6.0 * s6 * C6 * r4 / (t6 * t6);
            tmp8 = 8.0 * s8 * C8 * r6 / (t8 * t8);

            fpair = -(tmp6 + tmp8);
            fpair *= factor_lj;
          }
        }

        if (eflag) { evdwl = -(s6 * e6 + s8 * e8) * factor_lj; }

        double rest = (s6 * e6 + s8 * e8) / C6;

        dc6[i] += rest * c6_res[1];
        if (newton_pair || j < nlocal) { dc6[j] += rest * c6_res[2]; }

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  // communicate derivatives of C6
  communicationStage = 2;
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);

  // After calculating all derivatives dE/dr_ij w.r.t. distances,
  // the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
  for (int ii = 0; ii < inum; ii++) {

    int i = ilist[ii];

    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    int jnum = numneigh[i];
    int *jlist = firstneigh[i];

    for (int jj = 0; jj < jnum; jj++) {

      int j = jlist[jj];
      double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];

      double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[type[i]][type[j]]) {

        double r = sqrt(rsq);
        double dcn;

        // here we calculate dcn = dCNi/dr = dCNj/dr
        if (rsq < cn_thr) {
          double rcovij = (rcov[type[i]] + rcov[type[j]]) * autoang;
          double expterm = exp(-K1 * (rcovij / r - 1.0));
          dcn = -K1 * rcovij * expterm / (rsq * (expterm + 1.0) * (expterm + 1.0));
        } else {
          dcn = 0.0;
        };

        double fpair = dcn * (dc6[i] + dc6[j]) / r;
        fpair *= factor_lj;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }
        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

void PairDispersionD3::set_funcpar(std::string &functional_name)
{

  std::unordered_map<std::string, int> dampingMap = {
      {"zero", 1}, {"zerom", 2}, {"bj", 3}, {"bjm", 4}};

  int dampingCode = dampingMap[damping_type];

  switch (dampingCode) {

    case 1: {    // zero

      s6 = 1.0;
      alpha = 14.0;
      rs8 = 1.0;

      // default def2-QZVP (almost basis set limit)
      std::unordered_map<std::string, int> functionalMap = {{"slater-dirac-exchange", 1},
                                                            {"b-lyp", 2},
                                                            {"b-p", 3},
                                                            {"b97-d", 4},
                                                            {"revpbe", 5},
                                                            {"pbe", 6},
                                                            {"pbesol", 7},
                                                            {"rpw86-pbe", 8},
                                                            {"rpbe", 9},
                                                            {"tpss", 10},
                                                            {"b3-lyp", 11},
                                                            {"pbe0", 12},
                                                            {"hse06", 13},
                                                            {"revpbe38", 14},
                                                            {"pw6b95", 15},
                                                            {"tpss0", 16},
                                                            {"b2-plyp", 17},
                                                            {"pwpb95", 18},
                                                            {"b2gp-plyp", 19},
                                                            {"ptpss", 20},
                                                            {"hf", 21},
                                                            {"mpwlyp", 22},
                                                            {"bpbe", 23},
                                                            {"bh-lyp", 24},
                                                            {"tpssh", 25},
                                                            {"pwb6k", 26},
                                                            {"b1b95", 27},
                                                            {"bop", 28},
                                                            {"o-lyp", 29},
                                                            {"o-pbe", 30},
                                                            {"ssb", 31},
                                                            {"revssb", 32},
                                                            {"otpss", 33},
                                                            {"b3pw91", 34},
                                                            {"revpbe0", 35},
                                                            {"pbe38", 36},
                                                            {"mpw1b95", 37},
                                                            {"mpwb1k", 38},
                                                            {"bmk", 39},
                                                            {"cam-b3lyp", 40},
                                                            {"lc-wpbe", 41},
                                                            {"m05", 42},
                                                            {"m052x", 43},
                                                            {"m06l", 44},
                                                            {"m06", 45},
                                                            {"m062x", 46},
                                                            {"m06hf", 47},
                                                            {"hcth120", 48}};

      std::transform(functional_name.begin(), functional_name.end(), functional_name.begin(),
                     ::tolower);

      int functionalCode = functionalMap[functional_name];

      switch (functionalCode) {
        case 1:
          rs6 = 0.999;
          s8 = -1.957;
          rs8 = 0.697;
          break;
        case 2:
          rs6 = 1.094;
          s8 = 1.682;
          break;
        case 3:
          rs6 = 1.139;
          s8 = 1.683;
          break;
        case 4:
          rs6 = 0.892;
          s8 = 0.909;
          break;
        case 5:
          rs6 = 0.923;
          s8 = 1.010;
          break;
        case 6:
          rs6 = 1.217;
          s8 = 0.722;
          break;
        case 7:
          rs6 = 1.345;
          s8 = 0.612;
          break;
        case 8:
          rs6 = 1.224;
          s8 = 0.901;
          break;
        case 9:
          rs6 = 0.872;
          s8 = 0.514;
          break;
        case 10:
          rs6 = 1.166;
          s8 = 1.105;
          break;
        case 11:
          rs6 = 1.261;
          s8 = 1.703;
          break;
        case 12:
          rs6 = 1.287;
          s8 = 0.928;
          break;
        case 13:
          rs6 = 1.129;
          s8 = 0.109;
          break;
        case 14:
          rs6 = 1.021;
          s8 = 0.862;
          break;
        case 15:
          rs6 = 1.532;
          s8 = 0.862;
          break;
        case 16:
          rs6 = 1.252;
          s8 = 1.242;
          break;
        case 17:
          rs6 = 1.427;
          s8 = 1.022;
          s6 = 0.64;
          break;
        case 18:
          rs6 = 1.557;
          s8 = 0.705;
          s6 = 0.82;
          break;
        case 19:
          rs6 = 1.586;
          s8 = 0.760;
          s6 = 0.56;
          break;
        case 20:
          rs6 = 1.541;
          s8 = 0.879;
          s6 = 0.75;
          break;
        case 21:
          rs6 = 1.158;
          s8 = 1.746;
          break;
        case 22:
          rs6 = 1.239;
          s8 = 1.098;
          break;
        case 23:
          rs6 = 1.087;
          s8 = 2.033;
          break;
        case 24:
          rs6 = 1.370;
          s8 = 1.442;
          break;
        case 25:
          rs6 = 1.223;
          s8 = 1.219;
          break;
        case 26:
          rs6 = 1.660;
          s8 = 0.550;
          break;
        case 27:
          rs6 = 1.613;
          s8 = 1.868;
          break;
        case 28:
          rs6 = 0.929;
          s8 = 1.975;
          break;
        case 29:
          rs6 = 0.806;
          s8 = 1.764;
          break;
        case 30:
          rs6 = 0.837;
          s8 = 2.055;
          break;
        case 31:
          rs6 = 1.215;
          s8 = 0.663;
          break;
        case 32:
          rs6 = 1.221;
          s8 = 0.560;
          break;
        case 33:
          rs6 = 1.128;
          s8 = 1.494;
          break;
        case 34:
          rs6 = 1.176;
          s8 = 1.775;
          break;
        case 35:
          rs6 = 0.949;
          s8 = 0.792;
          break;
        case 36:
          rs6 = 1.333;
          s8 = 0.998;
          break;
        case 37:
          rs6 = 1.605;
          s8 = 1.118;
          break;
        case 38:
          rs6 = 1.671;
          s8 = 1.061;
          break;
        case 39:
          rs6 = 1.931;
          s8 = 2.168;
          break;
        case 40:
          rs6 = 1.378;
          s8 = 1.217;
          break;
        case 41:
          rs6 = 1.355;
          s8 = 1.279;
          break;
        case 42:
          rs6 = 1.373;
          s8 = 0.595;
          break;
        case 43:
          rs6 = 1.417;
          s8 = 0.000;
          break;
        case 44:
          rs6 = 1.581;
          s8 = 0.000;
          break;
        case 45:
          rs6 = 1.325;
          s8 = 0.000;
          break;
        case 46:
          rs6 = 1.619;
          s8 = 0.000;
          break;
        case 47:
          rs6 = 1.446;
          s8 = 0.000;
          break;
        /* DFTB3(zeta = 4.0), old deprecated parameters; case ("dftb3"); rs6 = 1.235; s8 = 0.673; */
        case 48:
          rs6 = 1.221;
          s8 = 1.206;
          break;
        default:
          error->all(FLERR, "Functional name unknown");
          break;
      }
      //fprintf(stderr,"s6    : %f\n", s6);
      //fprintf(stderr,"s8    : %f\n", s8);
      //fprintf(stderr,"rs6   : %f\n", rs6);
      //fprintf(stderr,"rs8   : %f\n", rs8);
      //fprintf(stderr,"alpha : %f\n", alpha);
    } break;

    case 2: {    // zerom
      s6 = 1.0;
      alpha = 14.0;

      std::unordered_map<std::string, int> functionalMap = {
          {"b2-plyp", 1}, {"b3-lyp", 2}, {"b97-d", 3}, {"b-lyp", 4},
          {"b-p", 5},     {"pbe", 6},    {"pbe0", 7},  {"lc-wpbe", 8}};

      int functionalCode = functionalMap[functional_name];
      switch (functionalCode) {
        case 1:
          rs6 = 1.313134;
          s8 = 0.717543;
          rs8 = 0.016035;
          s6 = 0.640000;
          break;
        case 2:
          rs6 = 1.338153;
          s8 = 1.532981;
          rs8 = 0.013988;
          break;
        case 3:
          rs6 = 1.151808;
          s8 = 1.020078;
          rs8 = 0.035964;
          break;
        case 4:
          rs6 = 1.279637;
          s8 = 1.841686;
          rs8 = 0.014370;
          break;
        case 5:
          rs6 = 1.233460;
          s8 = 1.945174;
          rs8 = 0.000000;
          break;
        case 6:
          rs6 = 2.340218;
          s8 = 0.000000;
          rs8 = 0.129434;
          break;
        case 7:
          rs6 = 2.077949;
          s8 = 0.000081;
          rs8 = 0.116755;
          break;
        case 8:
          rs6 = 1.366361;
          s8 = 1.280619;
          rs8 = 0.003160;
          break;
        default:
          error->all(FLERR, "Functional name unknown");
          break;
      }
      //fprintf(stderr,"s6    : %f\n", s6);
      //fprintf(stderr,"s8    : %f\n", s8);
      //fprintf(stderr,"rs6   : %f\n", rs6);
      //fprintf(stderr,"rs8   : %f\n", rs8);
      //fprintf(stderr,"alpha : %f\n", alpha);

      rs8 = rs8 / autoang;
    } break;

    case 3: {    // bj

      s6 = 1.0;
      alpha = 14.0;

      std::unordered_map<std::string, int> functionalMap = {
          {"b-p", 1},       {"b-lyp", 2},        {"revpbe", 3},    {"rpbe", 4},
          {"b97-d", 5},     {"pbe", 6},          {"rpw86-pbe", 7}, {"b3-lyp", 8},
          {"tpss", 9},      {"hf", 10},          {"tpss0", 11},    {"pbe0", 12},
          {"hse06", 13},    {"revpbe38", 14},    {"pw6b95", 15},   {"b2-plyp", 16},
          {"dsd-blyp", 17}, {"dsd-blyp-fc", 18}, {"bop", 19},      {"mpwlyp", 20},
          {"o-lyp", 21},    {"pbesol", 22},      {"bpbe", 23},     {"opbe", 24},
          {"ssb", 25},      {"revssb", 26},      {"otpss", 27},    {"b3pw91", 28},
          {"bh-lyp", 29},   {"revpbe0", 30},     {"tpssh", 31},    {"mpw1b95", 32},
          {"pwb6k", 33},    {"b1b95", 34},       {"bmk", 35},      {"cam-b3lyp", 36},
          {"lc-wpbe", 37},  {"b2gp-plyp", 38},   {"ptpss", 39},    {"pwpb95", 40},
          {"hf/mixed", 41}, {"hf/sv", 42},       {"hf/minis", 43}, {"b3-lyp/6-31gd", 44},
          {"hcth120", 45},  {"pw1pw", 46},       {"pwgga", 47},    {"hsesol", 48},
          {"hf3c", 49},     {"hf3cv", 50},       {"pbeh3c", 51},   {"pbeh-3c", 52}};

      int functionalCode = functionalMap[functional_name];
      switch (functionalCode) {
        case 1:
          a1 = 0.3946;
          s8 = 3.2822;
          a2 = 4.8516;
          break;
        case 2:
          a1 = 0.4298;
          s8 = 2.6996;
          a2 = 4.2359;
          break;
        case 3:
          a1 = 0.5238;
          s8 = 2.3550;
          a2 = 3.5016;
          break;
        case 4:
          a1 = 0.1820;
          s8 = 0.8318;
          a2 = 4.0094;
          break;
        case 5:
          a1 = 0.5545;
          s8 = 2.2609;
          a2 = 3.2297;
          break;
        case 6:
          a1 = 0.4289;
          s8 = 0.7875;
          a2 = 4.4407;
          break;
        case 7:
          a1 = 0.4613;
          s8 = 1.3845;
          a2 = 4.5062;
          break;
        case 8:
          a1 = 0.3981;
          s8 = 1.9889;
          a2 = 4.4211;
          break;
        case 9:
          a1 = 0.4535;
          s8 = 1.9435;
          a2 = 4.4752;
          break;
        case 10:
          a1 = 0.3385;
          s8 = 0.9171;
          a2 = 2.8830;
          break;
        case 11:
          a1 = 0.3768;
          s8 = 1.2576;
          a2 = 4.5865;
          break;
        case 12:
          a1 = 0.4145;
          s8 = 1.2177;
          a2 = 4.8593;
          break;
        case 13:
          a1 = 0.383;
          s8 = 2.310;
          a2 = 5.685;
          break;
        case 14:
          a1 = 0.4309;
          s8 = 1.4760;
          a2 = 3.9446;
          break;
        case 15:
          a1 = 0.2076;
          s8 = 0.7257;
          a2 = 6.3750;
          break;
        case 16:
          a1 = 0.3065;
          s8 = 0.9147;
          a2 = 5.0570;
          s6 = 0.64;
          break;
        case 17:
          a1 = 0.0000;
          s8 = 0.2130;
          a2 = 6.0519;
          s6 = 0.50;
          break;
        case 18:
          a1 = 0.0009;
          s8 = 0.2112;
          a2 = 5.9807;
          s6 = 0.50;
          break;
        case 19:
          a1 = 0.4870;
          s8 = 3.2950;
          a2 = 3.5043;
          break;
        case 20:
          a1 = 0.4831;
          s8 = 2.0077;
          a2 = 4.5323;
          break;
        case 21:
          a1 = 0.5299;
          s8 = 2.6205;
          a2 = 2.8065;
          break;
        case 22:
          a1 = 0.4466;
          s8 = 2.9491;
          a2 = 6.1742;
          break;
        case 23:
          a1 = 0.4567;
          s8 = 4.0728;
          a2 = 4.3908;
          break;
        case 24:
          a1 = 0.5512;
          s8 = 3.3816;
          a2 = 2.9444;
          break;
        case 25:
          a1 = -0.0952;
          s8 = -0.1744;
          a2 = 5.2170;
          break;
        case 26:
          a1 = 0.4720;
          s8 = 0.4389;
          a2 = 4.0986;
          break;
        case 27:
          a1 = 0.4634;
          s8 = 2.7495;
          a2 = 4.3153;
          break;
        case 28:
          a1 = 0.4312;
          s8 = 2.8524;
          a2 = 4.4693;
          break;
        case 29:
          a1 = 0.2793;
          s8 = 1.0354;
          a2 = 4.9615;
          break;
        case 30:
          a1 = 0.4679;
          s8 = 1.7588;
          a2 = 3.7619;
          break;
        case 31:
          a1 = 0.4529;
          s8 = 2.2382;
          a2 = 4.6550;
          break;
        case 32:
          a1 = 0.1955;
          s8 = 1.0508;
          a2 = 6.4177;
          break;
        case 33:
          a1 = 0.1805;
          s8 = 0.9383;
          a2 = 7.7627;
          break;
        case 34:
          a1 = 0.2092;
          s8 = 1.4507;
          a2 = 5.5545;
          break;
        case 35:
          a1 = 0.1940;
          s8 = 2.0860;
          a2 = 5.9197;
          break;
        case 36:
          a1 = 0.3708;
          s8 = 2.0674;
          a2 = 5.4743;
          break;
        case 37:
          a1 = 0.3919;
          s8 = 1.8541;
          a2 = 5.0897;
          break;
        case 38:
          a1 = 0.0000;
          s8 = 0.2597;
          a2 = 6.3332;
          s6 = 0.560;
          break;
        case 39:
          a1 = 0.0000;
          s8 = 0.2804;
          a2 = 6.5745;
          s6 = 0.750;
          break;
        case 40:
          a1 = 0.0000;
          s8 = 0.2904;
          a2 = 7.3141;
          s6 = 0.820;
          break;
        // special HF / DFT with eBSSE correction;
        case 41:
          a1 = 0.5607;
          s8 = 3.9027;
          a2 = 4.5622;
          break;
        case 42:
          a1 = 0.4249;
          s8 = 2.1849;
          a2 = 4.2783;
          break;
        case 43:
          a1 = 0.1702;
          s8 = 0.9841;
          a2 = 3.8506;
          break;
        case 44:
          a1 = 0.5014;
          s8 = 4.0672;
          a2 = 4.8409;
          break;
        case 45:
          a1 = 0.3563;
          s8 = 1.0821;
          a2 = 4.3359;
          break;
        /*     DFTB3 old, deprecated parameters : ;
            *     case ("dftb3"); a1 = 0.7461; s8 = 3.209; a2 = 4.1906;
            *     special SCC - DFTB parametrization;
            *     full third order DFTB, self consistent charges, hydrogen pair damping with; exponent 4.2;
        */
        case 46:
          a1 = 0.3807;
          s8 = 2.3363;
          a2 = 5.8844;
          break;
        case 47:
          a1 = 0.2211;
          s8 = 2.6910;
          a2 = 6.7278;
          break;
        case 48:
          a1 = 0.4650;
          s8 = 2.9215;
          a2 = 6.2003;
          break;
        // special HF - D3 - gCP - SRB / MINIX parametrization;
        case 49:
          a1 = 0.4171;
          s8 = 0.8777;
          a2 = 2.9149;
          break;
        // special HF - D3 - gCP - SRB2 / ECP - 2G parametrization;
        case 50:
          a1 = 0.3063;
          s8 = 0.5022;
          a2 = 3.9856;
          break;
        // special PBEh - D3 - gCP / def2 - mSVP parametrization;
        case 51:
          a1 = 0.4860;
          s8 = 0.0000;
          a2 = 4.5000;
          break;
        case 52:
          a1 = 0.4860;
          s8 = 0.0000;
          a2 = 4.5000;
          break;
        default:
          error->all(FLERR, "Functional name unknown");
          break;
      }

      //fprintf(stderr,"s6    : %f\n", s6);
      //fprintf(stderr,"s8    : %f\n", s8);
      //fprintf(stderr,"a1    : %f\n", a1);
      //fprintf(stderr,"a2    : %f\n", a2);
      //fprintf(stderr,"alpha : %f\n", alpha);

      a2 = a2 * autoang;
    } break;

    case 4: {    // bjm

      s6 = 1.0;
      alpha = 14.0;

      std::unordered_map<std::string, int> functionalMap = {
          {"b2-plyp", 1}, {"b3-lyp", 2}, {"b97-d", 3}, {"b-lyp", 4},
          {"b-p", 5},     {"pbe", 6},    {"pbe0", 7},  {"lc-wpbe", 8}};

      int functionalCode = functionalMap[functional_name];
      switch (functionalCode) {
        case 1:
          a1 = 0.486434;
          s8 = 0.672820;
          a2 = 3.656466;
          s6 = 0.640000;
          break;
        case 2:
          a1 = 0.278672;
          s8 = 1.466677;
          a2 = 4.606311;
          break;
        case 3:
          a1 = 0.240184;
          s8 = 1.206988;
          a2 = 3.864426;
          break;
        case 4:
          a1 = 0.448486;
          s8 = 1.875007;
          a2 = 3.610679;
          break;
        case 5:
          a1 = 0.821850;
          s8 = 3.140281;
          a2 = 2.728151;
          break;
        case 6:
          a1 = 0.012092;
          s8 = 0.358940;
          a2 = 5.938951;
          break;
        case 7:
          a1 = 0.007912;
          s8 = 0.528823;
          a2 = 6.162326;
          break;
        case 8:
          a1 = 0.563761;
          s8 = 0.906564;
          a2 = 3.593680;
          break;
        default:
          error->all(FLERR, "Functional name unknown");
          break;
      }

      //fprintf(stderr,"s6    : %f\n", s6);
      //fprintf(stderr,"s8    : %f\n", s8);
      //fprintf(stderr,"a1    : %f\n", a1);
      //fprintf(stderr,"a2    : %f\n", a2);
      //fprintf(stderr,"alpha : %f\n", alpha);

      a2 = a2 * autoang;

    } break;
    default:
      error->all(FLERR, "Damping type unknown");
      break;
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDispersionD3::init_one(int i, int j)
{

  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  r0ab[j][i] = r0ab[i][j];

  return std::sqrt(rthr);
}

void PairDispersionD3::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style D3 requires atom IDs");
  //if (force->newton_pair == 0)
  //  error->all(FLERR,"Pair style D3 requires newton pair on");

  // need an half neighbor list
  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   Communication
------------------------------------------------------------------------- */

int PairDispersionD3::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  if (communicationStage == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = cn[j];
    }
  }
  if (communicationStage == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = dc6[j];
    }
  }

  return m;
}

void PairDispersionD3::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (communicationStage == 1) {
    for (i = first; i < last; i++) { cn[i] = buf[m++]; }
  }
  if (communicationStage == 2) {
    for (i = first; i < last; i++) { dc6[i] = buf[m++]; }
  }
}

int PairDispersionD3::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (communicationStage == 1) {
    for (i = first; i < last; i++) { buf[m++] = cn[i]; }
  }
  if (communicationStage == 2) {
    for (i = first; i < last; i++) { buf[m++] = dc6[i]; }
  }
  return m;
}

void PairDispersionD3::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  if (communicationStage == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      cn[j] += buf[m++];
    }
  }
  if (communicationStage == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      dc6[j] += buf[m++];
    }
  }
}
