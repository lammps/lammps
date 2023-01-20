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

#include "compute_nbond_atom.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNBondAtom::ComputeNBondAtom(LAMMPS *_lmp, int narg, char **arg) :
    Compute(_lmp, narg, arg), nbond(nullptr)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "compute nbond/atom", error);

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNBondAtom::~ComputeNBondAtom()
{
  memory->destroy(nbond);
}

/* ---------------------------------------------------------------------- */

void ComputeNBondAtom::compute_peratom()
{
  // grow local nbond array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(nbond);
    nmax = atom->nmax;
    memory->create(nbond, nmax, "nbond/atom:nbond");
    vector_atom = nbond;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  int nlocal = atom->nlocal;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;

  int ntotal = nlocal;
  if (force->newton) ntotal += atom->nghost;

  // set local nbond array
  int i, j, k;
  int *num_bond = atom->num_bond;
  int newton_bond = force->newton_bond;

  for (i = 0; i < ntotal; i++) nbond[i] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] <= 0) continue;

      k = atom->map(bond_atom[i][j]);
      if (k < 0) continue;

      nbond[i] += 1;
      if (newton_bond) nbond[k] += 1;
    }
  }

  // communicate ghost nbond between neighbor procs
  if (force->newton) comm->reverse_comm(this);

  // zero nbond of atoms not in group
  // only do this after comm since ghost contributions must be included
  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) nbond[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputeNBondAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = nbond[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeNBondAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    nbond[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeNBondAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
