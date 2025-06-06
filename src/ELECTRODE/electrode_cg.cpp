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

#include "electrode_cg.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

#include <cassert>

using namespace LAMMPS_NS;

ElectrodeCG::ElectrodeCG(LAMMPS *lmp) :
    Fix(lmp, 0,
        std::vector<char *>{(char *) "fix_electrode_cg", (char *) "all", (char *) "electrode/cg"}
            .data()),
    ChargeSolver()
{
  setup = a_cached_flag = false;
  nstep = ncall = 0;
  comm_forward = 1;
  nmax = 0;
  memory->create(potential_i, nmax, "ElectrodeCG:potential_i");
  elyt_step = -1;
}

/* ---------------------------------------------------------------------- */

ElectrodeCG::~ElectrodeCG() noexcept
{
  memory->destroy(potential_i);
  if (comm->me == 0)
    utils::logmesg(lmp, "Average conjugate gradient steps: {:.3g}\n", nstep * 1. / ncall);
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::memory_use()
{
  double bytes = 0.;
  bytes += q_ele.capacity() * sizeof(double);
  bytes += taglist.capacity() * sizeof(tagint);
  bytes += iele_to_group.capacity() * sizeof(int);
  bytes += bvec.capacity() * sizeof(double);
  bytes += a_cached.capacity() * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int ElectrodeCG::setmask()
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::setup_solver(double cg_threshold, ElectrodeVector *vec)
{
  setup_cg(cg_threshold);
  elec_vec = vec;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::setup_cg(double cg_threshold)
{
  setup = true;
  evscale = force->qe2f / force->qqrd2e;
  threshold = cg_threshold;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::update_solver(std::vector<tagint> taglist_local, std::vector<int> iele_to_group)
{
  assert(setup);
  taglist = taglist_local;
  nele = taglist.size();    // local number of electrode atoms
  q_ele.resize(nele);
  MPI_Allreduce(&nele, &nele_world, 1, MPI_INT, MPI_SUM, world);
  this->iele_to_group = iele_to_group;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::set_elyt_pot(double *b_nall)
{
  elyt_step = update->ntimestep;
  bvec = pot_to_vector(b_nall);
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::solve(std::vector<double> v)
{
  assert(setup);
  assert(update->ntimestep == elyt_step);    // bvec is up to date
  a_cached_flag = false;
  ncall++;
  auto b = std::vector<double>(nele);
  for (int i = 0; i < nele; i++) b[i] = bvec[i] - evscale * v[iele_to_group[i]];
  predict_q();
  q_ele = constraint_projection(q_ele, true);
  auto r = add(b, ele_ele_interaction(q_ele));
  auto d = constraint_projection(r, false);
  double dot_old = dot_product(r, d);
  double delta = dot_old;
  for (int k = 0; k < nele_world && delta > threshold; k++, nstep++) {
    auto y = ele_ele_interaction(d);
    double alpha = dot_old / -dot_product(d, y);
    q_ele = add(q_ele, scale_vector(alpha, d));
    // prepare next step
    if ((k + 1) % 20 == 0) {
      // avoid shifting residual. This rarely happens.
      q_ele = constraint_projection(q_ele, true);
      r = add(b, ele_ele_interaction(q_ele));
    } else {
      r = add(r, scale_vector(alpha, std::move(y)));
    }
    auto p = constraint_projection(r, false);
    double dot_new = dot_product(r, p);
    d = add(std::move(p), scale_vector(dot_new / dot_old, d));
    delta = dot_product(r, d);
    dot_old = dot_new;
  }
  if ((delta > threshold) && (comm->me == 0)) error->warning(FLERR, "CG threshold not reached");
  return q_ele;
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::get_potential(int igroup)
{
  assert(update->ntimestep == elyt_step);    // bvec is up to date
  if (!a_cached_flag) {
    a_cached_flag = true;
    a_cached = ele_ele_interaction(q_ele);
  }
  double pot = 0.;
  int count = 0;
  for (int i = 0; i < nele; i++) {
    if (iele_to_group[i] == igroup) {
      pot += (a_cached[i] + bvec[i]) / evscale;
      count++;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &pot, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, world);
  return pot / count;
}

/* ----------------------------------------------------------------------
    predict with current charges
------------------------------------------------------------------------- */

void ElectrodeCG::predict_q()
{
  double *q = atom->q;
  for (int i = 0; i < nele; i++) q_ele[i] = q[atom->map(taglist[i])];
}

/* ----------------------------------------------------------------------
   Calculate ele ele interaction by calling ElectrodeVector and translate
   from nlocal indices to local electrode indices
------------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::ele_ele_interaction(const std::vector<double> &q_vec)
{
  MPI_Barrier(world);
  double mult_start = MPI_Wtime();
  set_charges(q_vec);
  if (atom->nmax > nmax) {
    memory->destroy(potential_i);
    nmax = atom->nmax;
    memory->create(potential_i, nmax, "ElectrodeCG:potential_i");
  }
  memset(potential_i, 0, atom->nmax * sizeof(double));
  elec_vec->compute_pot(potential_i);
  auto a = pot_to_vector(potential_i);
  MPI_Barrier(world);
  mult_time += MPI_Wtime() - mult_start;
  return a;
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::pot_to_vector(double *pot)
{
  auto vec = std::vector<double>(nele, 0.);
  for (int i = 0; i < nele; i++) vec[i] = pot[atom->map(taglist[i])];
  return vec;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::set_charges(std::vector<double> q_vec)
{
  double *q = atom->q;
  for (int i = 0; i < nele; i++) q[atom->map(taglist[i])] = q_vec[i];
  comm->forward_comm(this);
  //intel_pack_buffers(); // TODO
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::scale_vector(double alpha, std::vector<double> x)
{
  for (double &xi : x) xi *= alpha;
  return x;
}
/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::add(std::vector<double> a, std::vector<double> b)
{
  assert(((int) a.size() == nele) && ((int) b.size() == nele));
  for (int i = 0; i < nele; i++) a[i] += b[i];
  return a;
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::dot_product(std::vector<double> a, std::vector<double> b)
{
  assert(((int) a.size() == nele) && ((int) b.size() == nele));
  double out = 0.;
  for (int i = 0; i < nele; i++) out += a[i] * b[i];
  MPI_Allreduce(MPI_IN_PLACE, &out, 1, MPI_DOUBLE, MPI_SUM, world);
  return out;
}

/* ----------------------------------------------------------------------
   project into direction that conserves total charge (cf. Gingrich master thesis)
   or correct total electrode charge to qtotal if correction
------------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::constraint_projection(std::vector<double> x, bool correction)
{
  switch (constraint) {
    case ChargeConstraint::NONE:
      return x;
    case ChargeConstraint::SINGLE: {
      double sum = 0.;
      for (double xi : x) sum += xi;
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
      if (correction) sum -= qtotal;
      sum /= nele_world;
      for (double &xi : x) xi -= sum;
      return x;
    }
    case ChargeConstraint::GROUP: {
      int const n = x.size();
      int const ngroups = qtotal_group.size();
      auto counts = std::vector<int>(ngroups, 0);
      auto sums = std::vector<double>(ngroups, 0);
      for (int i = 0; i < n; i++) {
        int const g = iele_to_group[i];
        sums[g] += x[i];
        counts[g]++;
      }
      MPI_Allreduce(MPI_IN_PLACE, sums.data(), ngroups, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, counts.data(), ngroups, MPI_INT, MPI_SUM, world);
      for (int g = 0; g < ngroups; g++) {
        if (correction) sums[g] -= qtotal_group[g];
        sums[g] /= counts[g];
      }
      for (int i = 0; i < n; i++) x[i] -= sums[iele_to_group[i]];
      return x;
    }
    default:
      error->all(FLERR, "Constraint not implemented");
  }
}

/* ---------------------------------------------------------------------- */

int ElectrodeCG::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int const j = list[i];
    buf[m++] = atom->q[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::unpack_forward_comm(int n, int first, double *buf)
{
  int const last = first + n;
  for (int i = first, m = 0; i < last; i++) atom->q[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

std::vector<double> ElectrodeCG::compute_potentials()
{
  error->all(FLERR, "Method not implemented");
  return std::vector<double>(0, 0.);
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::get_sb_charges(int)
{
  error->all(FLERR, "Method not implemented");
  return 0.;
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::get_macro_capacitance(int, int)
{
  error->all(FLERR, "Method not implemented");
  return 0.;
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::get_macro_elastance(int, int)
{
  error->all(FLERR, "Method not implemented");
  return 0.;
}

/* ---------------------------------------------------------------------- */

double ElectrodeCG::vacuum_capacitance()
{
  error->all(FLERR, "Method not implemented");
  return 0.;
}

/* ---------------------------------------------------------------------- */

void ElectrodeCG::buffer_and_gather(double const * /*ivec*/, double * /*elevec*/)
{
  error->all(FLERR, "Method not implemented");
}
