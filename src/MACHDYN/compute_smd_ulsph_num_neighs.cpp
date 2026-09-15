/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#include "compute_smd_ulsph_num_neighs.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMDULSPHNumNeighs::ComputeSMDULSPHNumNeighs(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, 2, "Illegal compute smd/ulsph_num_neighs command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  numNeighsOutput = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMDULSPHNumNeighs::~ComputeSMDULSPHNumNeighs()
{
  memory->destroy(numNeighsOutput);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDULSPHNumNeighs::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("^smd/ulsph/num/neighs").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDULSPHNumNeighs::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (atom->nmax > nmax) {
    memory->destroy(numNeighsOutput);
    nmax = atom->nmax;
    memory->create(numNeighsOutput, nmax, "smd/ulsph/num/neighs:numNeighsRefConfigOutput");
    vector_atom = numNeighsOutput;
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int itmp = 0;
  int *numNeighs = (int *) force->pair->extract("smd/ulsph/numNeighs_ptr", itmp);
  if (numNeighs == nullptr) {
    error->all(FLERR, "compute smd/ulsph/num/neighs failed to access numNeighs array");
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      numNeighsOutput[i] = numNeighs[i];
    } else {
      numNeighsOutput[i] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDULSPHNumNeighs::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
