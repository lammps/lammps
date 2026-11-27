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

#include "compute_heat_flux_virial_tally.h"

#include "atom.h"
#include "comm.h"
#include "citeme.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

static const char cite_tally_pair_manybody[] =
    "compute tallying of general many-body forces: "
    "doi:10.1103/mtkk-kyyy\n\n"
    "@article{Poulos2026,\n"
    " title = {Exact formula and spectral decomposition of the heat flux in molecular dynamics \n"
    "          for arbitrary many-body potentials},\n"
    " author = {Poulos, Markos and Surblys, Donatas and Termentzidis, Konstantinos},\n"
    " journal = {Phys. Rev. B},\n"
    " volume = {113},\n"
    " issue = {4},\n"
    " pages = {045414},\n"
    " numpages = {10},\n"
    " year = {2026},\n"
    " month = {Jan},\n"
    " publisher = {American Physical Society},\n"
    " doi = {10.1103/mtkk-kyyy},\n"
    " url = {https://link.aps.org/doi/10.1103/mtkk-kyyy}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVirialTally::ComputeHeatFluxVirialTally(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute heat/flux/virial/tally", error);

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1)
    error->all(FLERR, "Could not find compute heat/flux/virial/tally second group ID {}", arg[3]);
  groupbit2 = group->bitmask[igroup2];

  scalar_flag = 1;
  vector_flag = 0;
  peratom_flag = 1;
  timeflag = 1;

  comm_reverse = size_peratom_cols = 3;
  extscalar = 1;
  peflag = 1;    // we need Pair::ev_tally() to be run

  did_setup = invoked_peratom = invoked_scalar = -1;
  nmax = -1;
  fatom = nullptr;

  // process optional args (n-body contribution flags)

  if (narg == 4) {
    two_bdflag = three_bdflag = four_bdflag = 1;
  } else {
    two_bdflag = three_bdflag = four_bdflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "two_body") == 0)
        two_bdflag = 1;
      else if (strcmp(arg[iarg], "three_body") == 0)
        three_bdflag = 1;
      else if (strcmp(arg[iarg], "four_body") == 0)
        four_bdflag = 1;
      else
        error->all(FLERR, "Illegal compute heat/flux/virial/tally command: unknown keyword {}", arg[iarg]);
      iarg++;
    }
  }
  if (lmp->citeme && force && force->pair && force->pair->manybody_flag)
    lmp->citeme->add(cite_tally_pair_manybody);
}

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVirialTally::~ComputeHeatFluxVirialTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(fatom);
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute heat/flux/virial/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->centroidstressflag != CENTROID_AVAIL &&
        (force->pair->single_enable == 0 || force->pair->manybody_flag))
      error->warning(FLERR, "Compute heat/flux/virial/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute heat/flux/virial/tally only called from pair style");
  }

  // error out if any atoms are in both groups
  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit) && (atom->mask[i] & groupbit2)) {
      if (atom->tag_enable) {
        error->all(FLERR,
                   "Atom {} belonging to both groups is not allowed"
                   " with compute heat/flux/virial/tally",
                   atom->tag[i]);
      } else {
        error->all(FLERR,
                   "Atom belonging to both groups is not allowed"
                   " with compute heat/flux/virial/tally");
      }
    }
  }

  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(fatom);
    nmax = atom->nmax;
    memory->create(fatom, nmax, size_peratom_cols, "heat/flux/virial/tally:fatom");
    array_atom = fatom;
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i)
    for (int j = 0; j < size_peratom_cols; ++j) fatom[i][j] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputeHeatFluxVirialTally::pair_tally_callback(int i, int j, int nlocal, int newton, double,
                                                     double, double fpair, double dx, double dy,
                                                     double dz)
{
  if (two_bdflag == 0) return;

  const int *const mask = atom->mask;

  const bool ingroup1 = (mask[i] & groupbit);
  if ((ingroup1 && (mask[j] & groupbit2)) || ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    // signs set to calculate heat flux from group2 into group1
    const double fx = (ingroup1 ? 0.5 : -0.5) * dx * fpair;
    const double fy = (ingroup1 ? 0.5 : -0.5) * dy * fpair;
    const double fz = (ingroup1 ? 0.5 : -0.5) * dz * fpair;

    if (newton || i < nlocal) {
      fatom[i][0] += fx;
      fatom[i][1] += fy;
      fatom[i][2] += fz;
    }
    if (newton || j < nlocal) {
      fatom[j][0] += fx;
      fatom[j][1] += fy;
      fatom[j][2] += fz;
    }
  }
}

/* -----------------------------------  CENTROID  ADDITION  ----------------------------------------- */
/*
The H_ij is equal to (H(x_i-x_cs) - H(x_j-x_cs)) where H is the Heaviside step function and x_cs is the control surface coordinate. (Torii 2008)

- H_ij = 1 if i in group1 and j in group2
- H_ij = -1 if i in group2 and j in group1 (even if i and j belong to both groups; this is a problem of incorrect simulation setup, not of our implementation)
- H_ij = 0 otherwise (i.e if i and j are both in group1 or both in group2)

These are true even if i and j belong to both groups; this is then a problem of incorrect simulation setup, not of our implementation ;)
*/
void ComputeHeatFluxVirialTally::pair_cv_tally3_callback(int i, int j, int k, double *fi, double *fj, double *fk, double Uijk, double pi, double pj, double pk)
{
  if (three_bdflag == 0) return;

  const int *const mask = atom->mask;
  int H_ij = H_ab(i, j, mask, groupbit, groupbit2);
  int H_ik = H_ab(i, k, mask, groupbit, groupbit2);
  int H_jk = H_ab(j, k, mask, groupbit, groupbit2);

  for (int m=0; m<3; m++) {
    fatom[i][m] += (pj*H_ij + pk*H_ik)*fi[m];
    fatom[j][m] += (-pi*H_ij + pk*H_jk)*fj[m];
    fatom[k][m] += (-pi*H_ik - pj*H_jk)*fk[m];
  }

}

void ComputeHeatFluxVirialTally::pair_cv_tally4_callback(int i, int j, int k, int l, double *fi, double *fj, double *fk, double *fl, double Uijkl, double pi, double pj, double pk, double pl)
{
  if (four_bdflag == 0) return;

  const int *const mask = atom->mask;
  int H_ij = H_ab(i, j, mask, groupbit, groupbit2);
  int H_ik = H_ab(i, k, mask, groupbit, groupbit2);
  int H_il = H_ab(i, l, mask, groupbit, groupbit2);
  int H_jk = H_ab(j, k, mask, groupbit, groupbit2);
  int H_jl = H_ab(j, l, mask, groupbit, groupbit2);
  int H_kl = H_ab(k, l, mask, groupbit, groupbit2);

  for (int m=0; m<3; m++) {
    fatom[i][m] += (pj*H_ij + pk*H_ik + pl*H_il)*fi[m];
    fatom[j][m] += (-pi*H_ij + pk*H_jk + pl*H_jl)*fj[m];
    fatom[k][m] += (-pi*H_ik - pj*H_jk + pl*H_kl)*fk[m];
    fatom[l][m] += (-pi*H_il - pj*H_jl - pk*H_kl)*fl[m];
  }
}

/* -----------------------------------  CENTROID  ADDITION  ----------------------------------------- */



int ComputeHeatFluxVirialTally::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = fatom[i][0];
    buf[m++] = fatom[i][1];
    buf[m++] = fatom[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    fatom[j][0] += buf[m++];
    fatom[j][1] += buf[m++];
    fatom[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeHeatFluxVirialTally::compute_scalar()
{
  // need to collect contributions from ghost atoms
  // because we need to use local velocities to compute heat flux
  if (invoked_peratom != update->ntimestep) compute_peratom();

  invoked_scalar = update->ntimestep;
  if ((did_setup != invoked_scalar) || (update->eflag_global != invoked_scalar))
    error->all(FLERR, Error::NOLASTLINE,
               "Stress was not tallied on needed timestep" + utils::errorurl(22));

  if ((comm->me == 0) && !force->pair->did_tally_callback())
    error->warning(FLERR, "Stress was not tallied by pair style" + utils::errorurl(11));

  // sum heat flux across procs
  double hflux = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    hflux +=
        fatom[i][0] * atom->v[i][0] + fatom[i][1] * atom->v[i][1] + fatom[i][2] * atom->v[i][2];

  MPI_Allreduce(&hflux, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if ((did_setup != invoked_peratom) || (update->eflag_global != invoked_peratom))
    error->all(FLERR, Error::NOLASTLINE,
               "Stress was not tallied on needed timestep" + utils::errorurl(22));

  if ((comm->me == 0) && !force->pair->did_tally_callback())
    error->warning(FLERR, "Stress was not tallied by pair style" + utils::errorurl(11));

  // collect contributions from ghost atoms

  if (force->newton_pair) {
    comm->reverse_comm(this);

    // clear out ghost atom data after it has been collected to local atoms
    const int nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; ++i)
      for (int j = 0; j < size_peratom_cols; ++j) fatom[i][j] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeHeatFluxVirialTally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double) size_peratom_cols * sizeof(double);
  return bytes;
}
