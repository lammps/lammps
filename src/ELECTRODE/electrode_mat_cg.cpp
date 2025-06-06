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

#include "electrode_mat_cg.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"

#include <cassert>

using namespace LAMMPS_NS;

ElectrodeMatCG::ElectrodeMatCG(LAMMPS *lmp) : ElectrodeCG(lmp)
{
  matrix_set = false;
}

/* ---------------------------------------------------------------------- */

ElectrodeMatCG::~ElectrodeMatCG() noexcept {}

/* ---------------------------------------------------------------------- */

double ElectrodeMatCG::memory_use()
{
  double bytes = ElectrodeCG::memory_use();
  bytes += qele_world.capacity() * sizeof(double);
  bytes += iele_local.capacity() * sizeof(double);
  bytes += tag_to_iele.bucket_count() * (sizeof(tagint) + sizeof(int));    // TODO check
  return bytes;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatCG::set_elastance(int nele_world, double **elastance)
{
  matrix_set = true;
  this->nele_world = nele_world;
  n_mat = nele_world;
  qele_world = std::vector<double>(nele_world);
  this->elastance = elastance;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatCG::setup_solver(double cg_threshold, std::unordered_map<tagint, int> tag_to_iele)
{
  assert(matrix_set);
  ElectrodeCG::setup_cg(cg_threshold);
  this->tag_to_iele = tag_to_iele;
}

/* ---------------------------------------------------------------------- */

void ElectrodeMatCG::update_solver(std::vector<tagint> taglist_local,
                                   std::vector<int> iele_to_group)
{
  ElectrodeCG::update_solver(taglist_local, iele_to_group);
  if (n_mat != nele_world) error->all(FLERR, "Number of electrode atoms has changed");
  iele_local.clear();
  iele_local.reserve(nele);
  for (tagint t : taglist_local) iele_local.push_back(tag_to_iele[t]);
}

/* ----------------------------------------------------------------------
   Calculate ele ele interaction by multiplying elastance matrix with charges
------------------------------------------------------------------------- */

std::vector<double> ElectrodeMatCG::ele_ele_interaction(const std::vector<double> &q_vec)
{
  MPI_Barrier(world);
  double mult_start = MPI_Wtime();
  auto a = std::vector<double>(nele, 0.);
  std::fill(qele_world.begin(), qele_world.end(), 0.);
  for (int i = 0; i < nele; i++) qele_world[iele_local[i]] = q_vec[i];
  MPI_Allreduce(MPI_IN_PLACE, qele_world.data(), nele_world, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < nele; i++) {
    double a_tmp = 0.;
    double *_noalias row = elastance[iele_local[i]];
    for (int j = 0; j < nele_world; j++) a_tmp += row[j] * qele_world[j];
    a[i] = a_tmp;
  }
  MPI_Barrier(world);
  mult_time += MPI_Wtime() - mult_start;
  return a;
}

