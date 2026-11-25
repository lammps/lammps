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

#include "compute_smd_hourglass_error.h"

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

ComputeSMDHourglassError::ComputeSMDHourglassError(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, 2, "Illegal compute smd/hourglass_error command");
  if (atom->smd_flag != 1)
    error->all(FLERR, 2,
               "compute smd/hourglass/error command requires atom_style with hourglass_error");
  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  hourglass_error_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMDHourglassError::~ComputeSMDHourglassError()
{
  memory->sfree(hourglass_error_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDHourglassError::init()
{

  if ((comm->me == 0) && (modify->get_compute_by_style("^smd/hourglass/error").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDHourglassError::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow output Vector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(hourglass_error_vector);
    nmax = atom->nmax;
    hourglass_error_vector =
        (double *) memory->smalloc(nmax * sizeof(double), "atom:hourglass_error_vector");
    vector_atom = hourglass_error_vector;
  }

  int itmp = 0;
  auto *hourglass_error = (double *) force->pair->extract("smd/tlsph/hourglass_error_ptr", itmp);
  if (hourglass_error == nullptr) {
    error->all(FLERR, Error::NOLASTLINE, "compute smd/hourglass/error failed to access hourglass_error array");
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      hourglass_error_vector[i] = hourglass_error[i];
    } else {
      hourglass_error_vector[i] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDHourglassError::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
