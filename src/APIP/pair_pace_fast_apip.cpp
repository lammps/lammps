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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#include "pair_pace_fast_apip.h"

#include "atom.h"
#include "atom_vec_apip.h"

#include <cstring>

using namespace LAMMPS_NS;

PairPACEFastAPIP::PairPACEFastAPIP(LAMMPS *lmp) : PairPACEAPIP(lmp) {}

/**
  * Set lambda_required based on lambda and lambda_const
  * @return true if this calculation is not required
  */

int PairPACEFastAPIP::check_abort_condition(double *lambda, double *lambda_const,
                                            int *lambda_required, int i)
{
  if ((lambda[i] == 0) && ((!lambda_thermostat) || (lambda_thermostat && lambda_const[i] == 0))) {
    lambda_required[i] |= ApipLambdaRequired::NO_SIMPLE;
    return 1;
  }
  lambda_required[i] |= ApipLambdaRequired::SIMPLE;
  return 0;
}

/**
  * @return prefactor lambda which is used for a fast ACE potential
  */

double PairPACEFastAPIP::compute_factor_lambda(double lambda)
{
  return lambda;
}

/**
  * @return atom->apip_e_fast which is used for a fast ACE potential
  */

double *PairPACEFastAPIP::get_e_ref_ptr()
{
  return atom->apip_e_fast;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACEFastAPIP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  dim = 0;
  if (strcmp(str, "pace/fast/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }
  return nullptr;
}
