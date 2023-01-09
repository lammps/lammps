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

#include "compute_dihedral.h"

#include "dihedral.h"
#include "dihedral_hybrid.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDihedral::ComputeDihedral(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), emine(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute dihedral command");

  vector_flag = 1;
  extvector = 1;
  peflag = 1;
  timeflag = 1;

  // check if dihedral style hybrid exists

  dihedral = dynamic_cast<DihedralHybrid *>(force->dihedral_match("hybrid"));
  if (!dihedral) error->all(FLERR, "Dihedral style for compute dihedral command is not hybrid");
  size_vector = nsub = dihedral->nstyles;

  emine = new double[nsub];
  vector = new double[nsub];
}

/* ---------------------------------------------------------------------- */

ComputeDihedral::~ComputeDihedral()
{
  delete[] emine;
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeDihedral::init()
{
  // recheck dihedral style in case it has been changed

  dihedral = dynamic_cast<DihedralHybrid *>(force->dihedral_match("hybrid"));
  if (!dihedral) error->all(FLERR, "Dihedral style for compute dihedral command is not hybrid");
  if (dihedral->nstyles != nsub)
    error->all(FLERR, "Dihedral style for compute dihedral command has changed");
}

/* ---------------------------------------------------------------------- */

void ComputeDihedral::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->eflag_global != invoked_vector)
    error->all(FLERR, "Energy was not tallied on needed timestep");

  for (int i = 0; i < nsub; i++) emine[i] = dihedral->styles[i]->energy;

  MPI_Allreduce(emine, vector, nsub, MPI_DOUBLE, MPI_SUM, world);
}
