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

#include "compute_uf3.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "uf3_bspline_basis3.h"
#include "update.h"

#include <algorithm>
#include <vector>
#include <climits>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

UF3Reader::UF3Reader(LAMMPS *lmp) : PairUF3(lmp) {}

void UF3Reader::load_settings_2body()
{
  char *deg = (char *) "2";
  settings(1, &deg);
}

void UF3Reader::load_coeffs(int narg, char **argv) { coeff(narg, argv); }

void UF3Reader::ensure_bsplines()
{
  if (!bsplines_created) create_bsplines();
}

double UF3Reader::get_cutsq(int i, int j) const { return cutsq[i][j]; }

int UF3Reader::knot_start_2b(int it, int jt, double r)
{
  return (this->*get_starting_index_2b)(it, jt, r);
}

int UF3Reader::ncoeff_2b(int i, int j) const { return n2b_coeff_array_size[i][j]; }

double UF3Reader::basis_2b_phi(int it, int jt, int lcoeff, double rsq, double rij)
{
  int ks = (this->*get_starting_index_2b)(it, jt, rij);
  if (ks < 0) return 0.0;
  if (lcoeff < ks - 3 || lcoeff > ks) return 0.0;
  if (lcoeff < 0 || lcoeff >= n2b_coeff_array_size[it][jt]) return 0.0;
  uf3_bspline_basis3 bs(lmp, &n2b_knots_array[it][jt][lcoeff], 1.0);
  double rth = rsq * rij;
  int slot = lcoeff - (ks - 3);
  if (slot == 0) return bs.eval0(rth, rsq, rij);
  if (slot == 1) return bs.eval1(rth, rsq, rij);
  if (slot == 2) return bs.eval2(rth, rsq, rij);
  return bs.eval3(rth, rsq, rij);
}

double UF3Reader::basis_2b_dphi_dr(int it, int jt, int lcoeff, double rsq, double rij)
{
  int ks = (this->*get_starting_index_2b)(it, jt, rij);
  if (ks < 0) return 0.0;
  if (lcoeff < ks - 3 || lcoeff > ks) return 0.0;
  if (lcoeff < 0 || lcoeff >= n2b_coeff_array_size[it][jt]) return 0.0;
  uf3_bspline_basis3 bs(lmp, &n2b_knots_array[it][jt][lcoeff], 1.0);
  double rth = rsq * rij;
  int slot = lcoeff - (ks - 3);
  if (slot == 0)
    return (3.0 * rsq) * bs.constants[3] + (2.0 * rij) * bs.constants[2] + bs.constants[1];
  if (slot == 1)
    return (3.0 * rsq) * bs.constants[7] + (2.0 * rij) * bs.constants[6] + bs.constants[5];
  if (slot == 2)
    return (3.0 * rsq) * bs.constants[11] + (2.0 * rij) * bs.constants[10] + bs.constants[9];
  return (3.0 * rsq) * bs.constants[15] + (2.0 * rij) * bs.constants[14] + bs.constants[13];
}

/* ---------------------------------------------------------------------- */

ComputeUF3::ComputeUF3(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), uf3_reader(nullptr), cutmax(0.0), lastcol(0), ncoeff(0),
    col_off_2b(nullptr), list(nullptr), c_pe(nullptr), c_virial(nullptr), uf3local(nullptr),
    uf3all(nullptr)
{
  if (narg < 4 + atom->ntypes)
    error->all(FLERR,
               "Illegal compute uf3 command: expected 'compute ID GROUP uf3 FILE TYPE1 TYPE2 ...' "
               "({} element symbols after file)",
               atom->ntypes);

  array_flag = 1;
  extarray = 0;

  auto *reader = new UF3Reader(lmp);
  uf3_reader = reader;
  reader->load_settings_2body();

  int ntypes = atom->ntypes;
  std::vector<char *> coeff_arg(3 + ntypes + 1, nullptr);
  coeff_arg[0] = (char *) "*";
  coeff_arg[1] = (char *) "*";
  coeff_arg[2] = arg[3];
  for (int t = 0; t < ntypes; t++) coeff_arg[3 + t] = arg[4 + t];

  reader->load_coeffs(3 + ntypes, coeff_arg.data());

  if (reader->has_3b_pot())
    error->all(FLERR, "compute uf3 currently supports only 2-body UF3 potentials (pair_style uf3 2)");

  cutmax = 0.0;
  for (int i = 1; i <= ntypes; i++)
    for (int j = 1; j <= ntypes; j++) cutmax = std::max(cutmax, reader->cut_ij(i, j));

  ncoeff = 0;
  memory->create(col_off_2b, ntypes + 1, ntypes + 1, "compute_uf3:col_off_2b");
  for (int i = 1; i <= ntypes; i++)
    for (int j = 1; j <= ntypes; j++) col_off_2b[i][j] = -1;

  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) {
      if (reader->setflag[i][j]) {
        col_off_2b[i][j] = ncoeff;
        col_off_2b[j][i] = ncoeff;
        ncoeff += reader->n2b_coeff_count(i, j);
      }
    }
  }

  const bigint natoms = atom->natoms;
  if (natoms > static_cast<bigint>(INT_MAX)) error->all(FLERR, "Too many atoms for compute uf3");
  int n = static_cast<int>(natoms);

  int ndims_force = 3;
  int ndims_virial = 6;
  int bik_rows = 1;
  size_array_rows = bik_rows + ndims_force * n + ndims_virial;
  size_array_cols = ncoeff + 1;
  lastcol = ncoeff;
}

/* ---------------------------------------------------------------------- */

ComputeUF3::~ComputeUF3()
{
  if (copymode) return;
  free_local();
  if (modify->find_compute(id_virial) != -1) modify->delete_compute(id_virial);
  delete uf3_reader;
  uf3_reader = nullptr;
  memory->destroy(col_off_2b);
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::free_local()
{
  memory->destroy(uf3local);
  memory->destroy(uf3all);
  uf3local = uf3all = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::build_column_offsets() {}

/* ---------------------------------------------------------------------- */

void ComputeUF3::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute uf3 requires a pair style be defined (e.g. pair_style zero)");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR, "Compute uf3 cutoff {} is longer than pairwise cutoff {}", cutmax,
               force->pair->cutforce);

  // Perpetual full list (same pattern as pair_style uf3). Occasional + build_one() can skip the
  // actual npair::build() when last_build >= lastcall, leaving inum/jlist stale and zero descriptors.
  neighbor->add_request(this, NeighConst::REQ_FULL);

  memory->create(uf3local, size_array_rows, size_array_cols, "compute_uf3:uf3local");
  memory->create(uf3all, size_array_rows, size_array_cols, "compute_uf3:uf3all");
  array = uf3all;

  uf3_reader->ensure_bsplines();

  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR, "Compute thermo_pe does not exist (required by compute uf3)");

  id_virial = id + std::string("_press");
  c_virial = modify->add_compute(id_virial + " all pressure NULL virial");
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::init_list(int /*id*/, NeighList *ptr) { list = ptr; }

/* ---------------------------------------------------------------------- */

void ComputeUF3::compute_array()
{
  invoked_array = update->ntimestep;

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int ic = 0; ic < size_array_cols; ic++) uf3local[irow][ic] = 0.0;

  compute_array_core();
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::compute_array_core()
{
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  const bigint natoms = atom->natoms;
  int n = static_cast<int>(natoms);

  for (int i = 0; i < nlocal; i++) {
    const tagint tg = atom->tag[i];
    if (tg < 1 || tg > n)
      error->all(FLERR, Error::NOLASTLINE,
                 "Compute uf3: atom tag {} is out of range [1,{}] (Natoms). "
                 "FitSNAP expects tags 1..N matching row layout.",
                 tg, n);
  }

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  const int bik_rows = 1;
  const int ndims_virial = 6;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    const int itype = type[i];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const int row_offset_i = bik_rows + 3 * (atom->tag[i] - 1);

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];
      if (rsq >= uf3_reader->get_cutsq(itype, jtype)) continue;

      double rij = sqrt(rsq);
      int row_offset_j = bik_rows + 3 * (atom->tag[j] - 1);
      int c0 = col_off_2b[itype][jtype];
      if (c0 < 0) continue;

      int ks = uf3_reader->knot_start_2b(itype, jtype, rij);
      if (ks < 0) continue;

      for (int lcoeff = ks - 3; lcoeff <= ks; lcoeff++) {
        if (lcoeff < 0 || lcoeff >= uf3_reader->ncoeff_2b(itype, jtype)) continue;
        double phi = uf3_reader->basis_2b_phi(itype, jtype, lcoeff, rsq, rij);
        double dphi_dr = uf3_reader->basis_2b_dphi_dr(itype, jtype, lcoeff, rsq, rij);
        int gcol = c0 + lcoeff;

        uf3local[0][gcol] += phi;

        double fac = dphi_dr / rij;
        double fx = delx * fac;
        double fy = dely * fac;
        double fz = delz * fac;

        uf3local[row_offset_i][gcol] += fx;
        uf3local[row_offset_i + 1][gcol] += fy;
        uf3local[row_offset_i + 2][gcol] += fz;

        uf3local[row_offset_j][gcol] -= fx;
        uf3local[row_offset_j + 1][gcol] -= fy;
        uf3local[row_offset_j + 2][gcol] -= fz;

        int vbase = bik_rows + 3 * n;
        uf3local[vbase + 0][gcol] += (fx * x[i][0] - fx * x[j][0]);
        uf3local[vbase + 1][gcol] += (fy * x[i][1] - fy * x[j][1]);
        uf3local[vbase + 2][gcol] += (fz * x[i][2] - fz * x[j][2]);
        uf3local[vbase + 3][gcol] += (fz * x[i][1] - fz * x[j][1]);
        uf3local[vbase + 4][gcol] += (fz * x[i][0] - fz * x[j][0]);
        uf3local[vbase + 5][gcol] += (fy * x[i][0] - fy * x[j][0]);
      }
    }
  }

  MPI_Allreduce(&uf3local[0][0], &uf3all[0][0], size_array_rows * size_array_cols, MPI_DOUBLE,
                MPI_SUM, world);

  for (int i = 0; i < bik_rows; i++) uf3all[i][lastcol] = 0.0;
  uf3all[0][lastcol] = c_pe->compute_scalar();

  for (int i = 0; i < nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = bik_rows + 3 * (iglobal - 1);
    uf3all[irow++][lastcol] = atom->f[i][0];
    uf3all[irow++][lastcol] = atom->f[i][1];
    uf3all[irow][lastcol] = atom->f[i][2];
  }

  c_virial->compute_vector();
  int iv = bik_rows + 3 * n;
  uf3all[iv++][lastcol] = c_virial->vector[0];
  uf3all[iv++][lastcol] = c_virial->vector[1];
  uf3all[iv++][lastcol] = c_virial->vector[2];
  uf3all[iv++][lastcol] = c_virial->vector[5];
  uf3all[iv++][lastcol] = c_virial->vector[4];
  uf3all[iv][lastcol] = c_virial->vector[3];
}

/* ---------------------------------------------------------------------- */

double ComputeUF3::memory_usage()
{
  double bytes = (double) size_array_rows * size_array_cols * sizeof(double);
  bytes += (double) size_array_rows * size_array_cols * sizeof(double);
  return bytes;
}
