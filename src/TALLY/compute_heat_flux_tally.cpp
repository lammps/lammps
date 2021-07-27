/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_heat_flux_tally.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHeatFluxTally::ComputeHeatFluxTally(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal compute heat/flux/tally command");

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1) error->all(FLERR, "Could not find compute heat/flux/tally second group ID");
  groupbit2 = group->bitmask[igroup2];

  vector_flag = 1;
  timeflag = 1;
  dynamic_group_allow = 0;

  comm_reverse = 7;
  extvector = 1;
  size_vector = 6;
  peflag = 1;    // we need Pair::ev_tally() to be run

  did_setup = 0;
  invoked_peratom = invoked_scalar = -1;
  nmax = -1;
  stress = nullptr;
  eatom = nullptr;
  vector = new double[size_vector];
  heatj = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeHeatFluxTally::~ComputeHeatFluxTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(stress);
  memory->destroy(eatom);
  delete[] heatj;
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxTally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute heat/flux/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR, "Compute heat/flux/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute heat/flux/tally only called from pair style");
  }
  did_setup = -1;
}

/* ---------------------------------------------------------------------- */
void ComputeHeatFluxTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(stress);
    memory->destroy(eatom);
    nmax = atom->nmax;
    memory->create(stress, nmax, 6, "heat/flux/tally:stress");
    memory->create(eatom, nmax, "heat/flux/tally:eatom");
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i) {
    eatom[i] = 0.0;
    stress[i][0] = 0.0;
    stress[i][1] = 0.0;
    stress[i][2] = 0.0;
    stress[i][3] = 0.0;
    stress[i][4] = 0.0;
    stress[i][5] = 0.0;
  }

  for (int i = 0; i < size_vector; ++i) vector[i] = heatj[i] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputeHeatFluxTally::pair_tally_callback(int i, int j, int nlocal, int newton, double evdwl,
                                               double ecoul, double fpair, double dx, double dy,
                                               double dz)
{
  const int *const mask = atom->mask;

  if (((mask[i] & groupbit) && (mask[j] & groupbit2)) ||
      ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    const double epairhalf = 0.5 * (evdwl + ecoul);
    fpair *= 0.5;
    const double v0 = dx * dx * fpair;    // dx*fpair = Fij_x
    const double v1 = dy * dy * fpair;
    const double v2 = dz * dz * fpair;
    const double v3 = dx * dy * fpair;
    const double v4 = dx * dz * fpair;
    const double v5 = dy * dz * fpair;

    if (newton || i < nlocal) {
      eatom[i] += epairhalf;
      stress[i][0] += v0;
      stress[i][1] += v1;
      stress[i][2] += v2;
      stress[i][3] += v3;
      stress[i][4] += v4;
      stress[i][5] += v5;
    }
    if (newton || j < nlocal) {
      eatom[j] += epairhalf;
      stress[j][0] += v0;
      stress[j][1] += v1;
      stress[j][2] += v2;
      stress[j][3] += v3;
      stress[j][4] += v4;
      stress[j][5] += v5;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeHeatFluxTally::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = eatom[i];
    buf[m++] = stress[i][0];
    buf[m++] = stress[i][1];
    buf[m++] = stress[i][2];
    buf[m++] = stress[i][3];
    buf[m++] = stress[i][4];
    buf[m++] = stress[i][5];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxTally::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    eatom[j] += buf[m++];
    stress[j][0] += buf[m++];
    stress[j][1] += buf[m++];
    stress[j][2] += buf[m++];
    stress[j][3] += buf[m++];
    stress[j][4] += buf[m++];
    stress[j][5] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxTally::compute_vector()
{
  invoked_vector = update->ntimestep;
  if ((did_setup != invoked_vector) || (update->eflag_global != invoked_vector))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // collect contributions from ghost atoms

  if (force->newton_pair) {
    comm->reverse_comm_compute(this);

    const int nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; ++i) {
      eatom[i] = 0.0;
      stress[i][0] = 0.0;
      stress[i][1] = 0.0;
      stress[i][2] = 0.0;
      stress[i][3] = 0.0;
      stress[i][4] = 0.0;
      stress[i][5] = 0.0;
    }
  }

  // compute heat currents
  // heat flux vector = jc[3] + jv[3]
  // jc[3] = convective portion of heat flux = sum_i (ke_i + pe_i) v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included
  // J = sum_i( (0.5*m*v_i^2 + 0.5*(evdwl_i+ecoul_i))*v_i +
  //              + (F_ij . v_i)*dR_ij/2 )

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  const double pfactor = 0.5 * force->mvv2e;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;

  double jc[3] = {0.0, 0.0, 0.0};
  double jv[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double *const vi = v[i];
      const double *const si = stress[i];
      double ke_i;

      if (rmass)
        ke_i = pfactor * rmass[i];
      else
        ke_i = pfactor * mass[type[i]];
      ke_i *= (vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]);
      ke_i += eatom[i];

      jc[0] += ke_i * vi[0];
      jc[1] += ke_i * vi[1];
      jc[2] += ke_i * vi[2];
      jv[0] += si[0] * vi[0] + si[3] * vi[1] + si[4] * vi[2];
      jv[1] += si[3] * vi[0] + si[1] * vi[1] + si[5] * vi[2];
      jv[2] += si[4] * vi[0] + si[5] * vi[1] + si[2] * vi[2];
    }
  }

  // sum accumulated heatj across procs
  heatj[0] = jc[0] + jv[0];
  heatj[1] = jc[1] + jv[1];
  heatj[2] = jc[2] + jv[2];
  heatj[3] = jc[0];
  heatj[4] = jc[1];
  heatj[5] = jc[2];
  MPI_Allreduce(heatj, vector, size_vector, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeHeatFluxTally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double)comm_reverse * sizeof(double);
  return bytes;
}
