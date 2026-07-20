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
   Contributing author: Claude Opus 4.8 (Anthropic), under the direction of
   Joel Clemmer (SNL) and Axel Kohlmeyer (Temple U). Derived from the PERI
   package compute_damage_atom.cpp (Mike Parks, SNL).
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Volume-weighted bond damage for the BPM peridynamics model:
       damage[i] = 1 - (sum of nodal volume over surviving bonds) / vinter[i]
   where vinter[i] is the reference interaction volume stored by bond_style
   bpm/peri. Damage is 0 for a fully intact node and approaches 1 as a node
   loses bonds. Companion diagnostic to pair/bond style bpm/peri.
------------------------------------------------------------------------- */

#include "compute_bpm_peri_damage_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBPMPeriDamageAtom::ComputeBPMPeriDamageAtom(LAMMPS *_lmp, int narg, char **arg) :
    Compute(_lmp, narg, arg), damage(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute bpm/peri/damage/atom command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR, "Compute bpm/peri/damage/atom used in system without bonds");

  peratom_flag = 1;
  size_peratom_cols = 0;

  index_vfrac = index_vinter = -1;
  nmax = 0;
}

/* ----------------------------------------------------------------------
   resolve the per-atom properties here (not in the constructor): vinter is
   created by bond_style bpm/peri during force->init(), which runs before the
   compute's init() in the setup sequence.
------------------------------------------------------------------------- */

void ComputeBPMPeriDamageAtom::init()
{
  // the nodal volume vfrac is user-supplied; the reference interaction volume
  // vinter is created and filled by bond_style bpm/peri
  int flag, cols;
  index_vfrac = atom->find_custom("vfrac", flag, cols);
  if ((index_vfrac < 0) || (flag != 1) || (cols != 0))
    error->all(FLERR, "Compute bpm/peri/damage/atom requires a per-atom vfrac property");
  index_vinter = atom->find_custom("vinter", flag, cols);
  if ((index_vinter < 0) || (flag != 1) || (cols != 0))
    error->all(FLERR,
               "Compute bpm/peri/damage/atom requires bond style bpm/peri (missing vinter)");
}

/* ---------------------------------------------------------------------- */

ComputeBPMPeriDamageAtom::~ComputeBPMPeriDamageAtom()
{
  memory->destroy(damage);
}

/* ---------------------------------------------------------------------- */

void ComputeBPMPeriDamageAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow damage array if necessary
  if (atom->nmax > nmax) {
    memory->destroy(damage);
    nmax = atom->nmax;
    memory->create(damage, nmax, "bpm/peri/damage/atom:damage");
    vector_atom = damage;
  }

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  double *vfrac = atom->dvector[index_vfrac];
  double *vinter = atom->dvector[index_vinter];

  // bond_style bpm/peri mandates newton bond off, so every bond touching atom i
  // is stored at i; the surviving interaction volume is a purely local sum.
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) {
      damage[i] = 0.0;
      continue;
    }

    double surviving = 0.0;
    for (int m = 0; m < num_bond[i]; m++) {
      if (bond_type[i][m] <= 0) continue;    // broken/turned-off bond
      int j = atom->map(bond_atom[i][m]);
      if (j < 0) continue;
      surviving += vfrac[j];
    }

    damage[i] = (vinter[i] > 0.0) ? 1.0 - surviving / vinter[i] : 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeBPMPeriDamageAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
