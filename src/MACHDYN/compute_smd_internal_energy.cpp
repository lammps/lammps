// clang-format off
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

#include <cstring>
#include "compute_smd_internal_energy.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMDInternalEnergy::ComputeSMDInternalEnergy(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute smd/internal_energy command");
  if (atom->esph_flag != 1) error->all(FLERR,"compute smd/internal_energy command requires atom_style with internal_energy (e.g. smd)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  internal_energy_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMDInternalEnergy::~ComputeSMDInternalEnergy()
{
  memory->sfree(internal_energy_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDInternalEnergy::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"smd/internal_energy") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute smd/internal_energy");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDInternalEnergy::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rhoVector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(internal_energy_vector);
    nmax = atom->nmax;
    internal_energy_vector = (double *) memory->smalloc(nmax*sizeof(double),"atom:internal_energy_vector");
    vector_atom = internal_energy_vector;
  }

  double *esph = atom->esph;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              internal_energy_vector[i] = esph[i];
      }
      else {
              internal_energy_vector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSMDInternalEnergy::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
