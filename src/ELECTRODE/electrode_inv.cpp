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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#include "electrode_inv.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

#include <cassert>
#include <unordered_map>

using namespace LAMMPS_NS;

static constexpr double SMALL = 1e-16;

extern "C" {
void dgetrf_(const int *M, const int *N, double *A, const int *lda, int *ipiv, int *info);
void dgetri_(const int *N, double *A, const int *lda, const int *ipiv, double *work,
             const int *lwork, int *info);
}

ElectrodeInv::ElectrodeInv(LAMMPS *lmp) : Pointers(lmp), ChargeSolver()
{
  setup = cap_set = vac_cap_computed = false;
  nmax = 0;
  memory->create(potential_i, nmax, "ElectrodeInv:potential_i");
  int const nprocs = comm->nprocs;
  recvcounts = new int[nprocs];
  displs = new int[nprocs];
  elyt_step = -1;
}

/* ---------------------------------------------------------------------- */

ElectrodeInv::~ElectrodeInv() noexcept
{
  memory->destroy(potential_i);
  if (setup) {
    memory->destroy(iele_gathered);
    memory->destroy(buf_gathered);
    memory->destroy(potential_iele);
  }
  delete[] recvcounts;
  delete[] displs;
}


/* ---------------------------------------------------------------------- */

double ElectrodeInv::memory_use()
{
  double bytes = 0.;
  if (setup) bytes += 3 * nele_world * sizeof(double);
  bytes += qvec.capacity() * sizeof(double);
  bytes += iele_to_group.capacity() * sizeof(int);
  bytes += taglist_local.capacity() * sizeof(tagint);
  bytes += iele_local.capacity() * sizeof(int);
  bytes += buf_iele.capacity() * sizeof(double);
  bytes += tag_to_iele.bucket_count() * (sizeof(tagint) + sizeof(int));    // TODO check
  return bytes;
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::set_elastance(int nele_world, double **elastance)
{
  cap_set = true;
  this->nele_world = nele_world;
  this->capacitance = elastance;
  // invert elastance to obtain capacitance
  MPI_Barrier(world);
  // TODO timer flag
  //double invert_time = MPI_Wtime();
  //if (timer_flag && (comm->me == 0)) utils::logmesg(lmp, "CONP inverting matrix\n");
  int m = nele_world, n = nele_world, lda = nele_world;
  std::vector<int> ipiv(nele_world);
  int const lwork = nele_world * nele_world;
  std::vector<double> work(lwork);

  int info_rf, info_ri;
  dgetrf_(&m, &n, &capacitance[0][0], &lda, ipiv.data(), &info_rf);
  dgetri_(&n, &capacitance[0][0], &lda, ipiv.data(), work.data(), &lwork, &info_ri);
  if (info_rf != 0 || info_ri != 0) error->all(FLERR, "CONP matrix inversion failed!");
  MPI_Barrier(world);
  //if (timer_flag && (comm->me == 0))
  //utils::logmesg(lmp, "Invert time: {:.4g} s\n", MPI_Wtime() - invert_time);
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::set_capacitance(int nele_world, double **capacitance)
{
  cap_set = true;
  this->nele_world = nele_world;
  this->capacitance = capacitance;
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::setup_solver(int groupbit, std::unordered_map<tagint, int> tag_to_iele,
                                std::vector<int> group_bits,  bool ffield)
{
  assert(cap_set);
  setup = true;
  evscale = force->qe2f / force->qqrd2e;
  this->groupbit = groupbit;
  if (ffield) symmetrize();
  this->tag_to_iele = tag_to_iele;
  memory->create(iele_gathered, nele_world, "ElectrodeInv:iele_gathered");
  memory->create(buf_gathered, nele_world, "ElectrodeInv:buf_gathered");
  memory->create(potential_iele, nele_world, "ElectrodeInv:potential_iele");
  int const nlocal = atom->nlocal;
  ngroups = group_bits.size();
  group_pot = std::vector<double>(ngroups, 0.);
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  iele_to_group = std::vector<int>(nele_world, -1);
  for (int i = 0; i < nlocal; i++) {
    for (int g = 0; g < ngroups; g++) {
      if (mask[i] & group_bits[g]) { iele_to_group[tag_to_iele[tag[i]]] = g; }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, iele_to_group.data(), nele_world, MPI_INT, MPI_MAX, world);
  sb_charges = std::vector<double>(ngroups);
  if (ffield) {
    compute_sd_vectors_ffield(group_bits);
  } else
    compute_sd_vectors();
  compute_macro_matrices(ffield);
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::update_solver(std::vector<tagint> taglist_local,
                                 std::vector<int> /*iele_to_group_local*/)
{
  assert(setup);
  nlocalele = taglist_local.size();
  this->taglist_local = taglist_local;
  int const nprocs = comm->nprocs;
  qvec = std::vector<double>(nlocalele);
  delete[] recvcounts;
  delete[] displs;
  recvcounts = new int[nprocs];
  displs = new int[nprocs];
  MPI_Allgather(&nlocalele, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];

  qvec.resize(nlocalele);
  iele_local.clear();
  iele_local.reserve(nlocalele);
  for (tagint t : taglist_local) iele_local.push_back(tag_to_iele[t]);
  MPI_Allgatherv(iele_local.data(), nlocalele, MPI_INT, iele_gathered, recvcounts, displs, MPI_INT,
                 world);
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::set_elyt_pot(double *b_nall)
{
  elyt_step = update->ntimestep;
  buffer_and_gather(b_nall, potential_iele);
  // calculate charges due to electrolyte
  std::fill(sb_charges.begin(), sb_charges.end(), 0.);
  MPI_Barrier(world);
  double mult_start = MPI_Wtime();
  for (int i = 0; i < nlocalele; i++) {
    double q_tmp = 0.;
    int const iele = iele_local[i];
    double *_noalias caprow = capacitance[iele];
    for (int j = 0; j < nele_world; j++) { q_tmp -= caprow[j] * potential_iele[j]; }
    sb_charges[iele_to_group[iele]] += q_tmp;
    qvec[i] = q_tmp;
  }
  MPI_Allreduce(MPI_IN_PLACE, sb_charges.data(), ngroups, MPI_DOUBLE, MPI_SUM, world);
  MPI_Barrier(world);
  mult_time += MPI_Wtime() - mult_start;
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeInv::solve(std::vector<double> v)
{
  assert(setup);
  assert(update->ntimestep == elyt_step);    // qvec already has charges due to electrolyte
  MPI_Barrier(world);
  double mult_start = MPI_Wtime();
  group_pot = apply_constraint(v);
  // calculate final charges
  for (int g = 0; g < ngroups; g++)
    for (int j = 0; j < nlocalele; j++) qvec[j] += sd_vectors[g][iele_local[j]] * group_pot[g];
  MPI_Barrier(world);
  mult_time += MPI_Wtime() - mult_start;
  return qvec;
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeInv::apply_constraint(std::vector<double> v)
{
  switch (constraint) {
    case ChargeConstraint::NONE:
      break;
    case ChargeConstraint::SINGLE: {
      double q_current = 0.;
      for (int i = 0; i < ngroups; i++) {
        q_current += sb_charges[i];
        for (int j = 0; j < ngroups; j++) q_current += macro_capacitance[i][j] * v[j];
      }
      double add_psi = (qtotal - q_current) / macro_capacitance_sum;
      for (int i = 0; i < ngroups; i++) v[i] += add_psi;
      break;
    }
    case ChargeConstraint::GROUP: {
      std::vector<double> group_remainder_q(ngroups);
      for (int g = 0; g < ngroups; g++) group_remainder_q[g] = qtotal_group[g] - sb_charges[g];
      for (int g = 0; g < ngroups; g++) {
        double vtmp = 0;
        for (int h = 0; h < ngroups; h++) { vtmp += macro_elastance[g][h] * group_remainder_q[h]; }
        v[g] = vtmp;
      }
      break;
    }
    default:
      error->all(FLERR, "Constraint not implemented");
  }
  return v;
}

/* ----------------------------------------------------------------------
    calculate potentials that would give current charges
------------------------------------------------------------------------- */

std::vector<double> ElectrodeInv::compute_potentials()
{
  assert(setup);
  assert(update->ntimestep == elyt_step);    // assert sb_charges up to date
  // sum charges for each group
  int *tag = atom->tag;
  double *q = atom->q;
  auto group_q = std::vector<double>(ngroups, 0.);
  for (int i = 0; i < atom->nlocal; i++) group_q[iele_to_group[tag_to_iele[tag[i]]]] += q[i];
  MPI_Allreduce(MPI_IN_PLACE, group_q.data(), ngroups, MPI_DOUBLE, MPI_SUM, world);

  // compute potentials
  for (int g = 0; g < ngroups; g++) group_q[g] -= sb_charges[g];
  for (int g = 0; g < ngroups; g++) {
    double vtmp = 0.;
    for (int h = 0; h < ngroups; h++) vtmp += macro_elastance[g][h] * group_q[h];
    group_pot[g] = vtmp;
  }
  return group_pot;
}

/* ---------------------------------------------------------------------- */

double ElectrodeInv::get_sb_charges(int igroup)
{
  assert(setup);
  assert(igroup < ngroups);
  assert(update->ntimestep == elyt_step);
  return sb_charges[igroup];
}

/* ---------------------------------------------------------------------- */

double ElectrodeInv::get_macro_capacitance(int igroup, int jgroup)
{
  assert(setup);
  assert(igroup < ngroups && jgroup < ngroups);
  return macro_capacitance[igroup][jgroup];
}

/* ---------------------------------------------------------------------- */

double ElectrodeInv::get_macro_elastance(int igroup, int jgroup)
{
  assert(setup);
  assert(igroup < ngroups && jgroup < ngroups);
  return macro_elastance[igroup][jgroup];
}

/* ---------------------------------------------------------------------- */

double ElectrodeInv::get_potential(int igroup)
{
  assert(setup);
  assert(igroup < ngroups);
  return group_pot[igroup];
}

/* ---------------------------------------------------------------------- */

double ElectrodeInv::vacuum_capacitance()
{
  assert(ngroups == 2);
  if (!vac_cap_computed) {
    vac_cap_computed = true;
    vac_cap = (macro_capacitance[0][0] * macro_capacitance[1][1] -
               macro_capacitance[0][1] * macro_capacitance[0][1]) /
        (macro_capacitance[0][0] + macro_capacitance[1][1] + 2 * macro_capacitance[0][1]);
  }
  return vac_cap;
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::buffer_and_gather(double const *ivec, double *elevec)
{
  buf_iele.resize(nlocalele);
  for (int i_iele = 0; i_iele < nlocalele; i_iele++) {
    buf_iele[i_iele] = ivec[atom->map(taglist_local[i_iele])];
  }
  MPI_Allgatherv(buf_iele.data(), nlocalele, MPI_DOUBLE, buf_gathered, recvcounts, displs,
                 MPI_DOUBLE, world);

  for (int i = 0; i < nele_world; i++) elevec[iele_gathered[i]] = buf_gathered[i];
}

/* ----------------------------------------------------------------------
    S matrix to enforce charge neutrality constraint
------------------------------------------------------------------------- */

void ElectrodeInv::symmetrize()
{

  std::vector<double> AinvE(nele_world, 0.);
  double EAinvE = 0.0;
  for (int i = 0; i < nele_world; i++) {
    double AinvEtmp = 0.0;
    for (int j = 0; j < nele_world; j++) AinvEtmp += capacitance[i][j];
    AinvE[i] = AinvEtmp;    // use temp accumulator to enable vectorization
    EAinvE += AinvEtmp;
  }
  for (int i = 0; i < nele_world; i++) {
    double iAinvE = AinvE[i];
    for (int j = 0; j < nele_world; j++) capacitance[i][j] -= AinvE[j] * iAinvE / EAinvE;
  }
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::compute_sd_vectors()
{
  sd_vectors = std::vector<std::vector<double>>(ngroups, std::vector<double>(nele_world, 0.));
  for (int g = 0; g < ngroups; g++) {
    for (int j = 0; j < nele_world; j++) {
      if (iele_to_group[j] == g) {
        for (int k = 0; k < nele_world; k++) { sd_vectors[g][k] += capacitance[k][j] * evscale; }
      }
    }
  }
}
/* ---------------------------------------------------------------------- */

void ElectrodeInv::compute_sd_vectors_ffield(std::vector<int> group_bits)
{
  int top_group = get_top_group(group_bits);
  sd_vectors = std::vector<std::vector<double>>(ngroups, std::vector<double>(nele_world, 0.));
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  double zprd = domain->prd[2];
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      int const i_iele = tag_to_iele[tag[i]];
      double const zprd_offset = (mask[i] & group_bits[top_group]) ? 0.0 : 1.0;
      double const evscale_elez = evscale * (x[i][2] / zprd + zprd_offset);
      for (int g = 0; g < ngroups; g++) {
        double gmult = (g == top_group) ? -1.0 : 1.0;
        for (int k = 0; k < nele_world; k++) {
          sd_vectors[g][k] += gmult * capacitance[k][i_iele] * evscale_elez;
        }
      }
    }
  }
  for (int g = 0; g < ngroups; g++) {
    MPI_Allreduce(MPI_IN_PLACE, sd_vectors[g].data(), nele_world, MPI_DOUBLE, MPI_SUM, world);
  }
}

/* ---------------------------------------------------------------------- */

int ElectrodeInv::get_top_group(std::vector<int> group_bits)
{
  double *zmax = new double[ngroups];
  double **x = atom->x;
  for (int g = 0; g < ngroups; g++) { zmax[g] = domain->boxlo[2]; }
  int *mask = atom->mask;
  for (int i = 0; i < atom->nlocal; i++) {
    for (int g = 0; g < ngroups; g++) {
      if (mask[i] & group_bits[g]) {
        if (x[i][2] > zmax[g]) zmax[g] = x[i][2];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, zmax, ngroups, MPI_DOUBLE, MPI_MAX, world);
  int gmax = 0;
  for (int g = 0; g < ngroups; g++) { gmax = (zmax[g] > zmax[gmax]) ? g : gmax; }
  delete[] zmax;
  return gmax;
}

/* ---------------------------------------------------------------------- */

void ElectrodeInv::compute_macro_matrices(bool symm)
{
  // capacitance
  macro_capacitance = std::vector<std::vector<double>>(ngroups, std::vector<double>(ngroups, 0.));
  for (int g = 0; g < ngroups; g++) {
    for (int k = 0; k < nele_world; k++) {
      macro_capacitance[iele_to_group[k]][g] += sd_vectors[g][k];
    }
  }
  if (symm) {    // scaling with C[0][0] improves numerical stability
    double scalar = macro_capacitance[0][0];
    macro_capacitance.back() = std::vector<double>(ngroups, scalar);
  }

  macro_capacitance_sum = 0.;
  for (int i = 0; i < ngroups; i++)
    for (int j = 0; j < ngroups; j++) macro_capacitance_sum += macro_capacitance[i][j];

  // elastance
  macro_elastance = std::vector<std::vector<double>>(ngroups, std::vector<double>(ngroups));
  switch (ngroups) {
    case 1: {
      macro_elastance[0][0] = 1. / macro_capacitance[0][0];
      break;
    }
    case 2: {
      double const det = macro_capacitance[0][0] * macro_capacitance[1][1] -
          macro_capacitance[0][1] * macro_capacitance[1][0];
      if (fabs(det) < SMALL) error->all(FLERR, "ELECTRODE macro matrix inversion failed!");
      double const detinv = 1 / det;
      macro_elastance[0][0] = macro_capacitance[1][1] * detinv;
      macro_elastance[1][1] = macro_capacitance[0][0] * detinv;
      macro_elastance[0][1] = -macro_capacitance[0][1] * detinv;
      macro_elastance[1][0] = -macro_capacitance[1][0] * detinv;
      break;
    }
    default:
      int m = ngroups;
      int n = m, lda = m;
      std::vector<int> ipiv(m);
      int const lwork = m * m;
      std::vector<double> work(lwork);
      std::vector<double> tmp(lwork);
      for (int i = 0; i < ngroups; i++) {
        for (int j = 0; j < ngroups; j++) {
          int idx = i * ngroups + j;
          tmp[idx] = macro_capacitance[i][j];
        }
      }
      int info_rf, info_ri;
      dgetrf_(&m, &n, tmp.data(), &lda, ipiv.data(), &info_rf);
      dgetri_(&n, tmp.data(), &lda, ipiv.data(), work.data(), &lwork, &info_ri);
      if (info_rf != 0 || info_ri != 0)
        error->all(FLERR, "ELECTRODE macro matrix inversion failed!");
      for (int i = 0; i < ngroups; i++) {
        for (int j = 0; j < ngroups; j++) {
          int idx = i * ngroups + j;
          macro_elastance[i][j] = tmp[idx];
        }
      }
      break;
  }
}
