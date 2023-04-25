// clang-format off
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
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "compute_rheo_property_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::ComputeRHEOPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  index(nullptr), colindex(nullptr), pack_choice(nullptr)
{
  if (narg < 4)  utils::missing_cmd_args(FLERR, "compute property/atom", error);

  peratom_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_id;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Compute property/atom {} is not available", arg[iarg]);
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_molecule;
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_proc;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_type;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[i] = &ComputeRHEOPropertyAtom::pack_mass;


    } else {
      error->all(FLERR,"Invalid keyword {} for compute rheo/property/atom command ", arg[iarg]);
    }
  }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOPropertyAtom::~ComputeRHEOPropertyAtom()
{
  delete[] pack_choice;
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow vector or array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    if (nvalues == 1) {
      memory->destroy(vector_atom);
      memory->create(vector_atom,nmax,"rheo/property/atom:vector");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom,nmax,nvalues,"rheo/property/atom:array");
    }
  }

  // fill vector or array with per-atom values

  if (nvalues == 1) {
    buf = vector_atom;
    (this->*pack_choice[0])(0);
  } else {
    if (nmax) buf = &array_atom[0][0];
    else buf = nullptr;
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOPropertyAtom::memory_usage()
{
  double bytes = (double)nmax * nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute rheo/property/atom can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_rho(int n)
{
  double *rho = atom->rho;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = rho[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOPropertyAtom::pack_status(int n)
{
  int *status = atom->status;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = status[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}
