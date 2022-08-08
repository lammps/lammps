// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

   This class contains a series of tools for DEM contacts
   Multiple models can be defined and used to calculate forces
   and torques based on contact geometry
*/

#include "contact_sub_models.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace Contact;

SubModel::SubModel() :
  Pointers(lmp)
{
  allocated = 0;
  size_history = 0;
  history_index = 0;
  allow_limit_damping = 1;
  beyond_contact = 0;
  num_coeffs = 0;

  nondefault_history_transfer = 0;
  transfer_history_factor = nullptr;
}

/* ---------------------------------------------------------------------- */

SubModel::~SubModel()
{
  if (allocated) delete [] coeffs;
  delete [] transfer_history_factor;
}

/* ---------------------------------------------------------------------- */

void SubModel::allocate_coeffs()
{
  allocated = 1;
  coeffs = new double[num_coeffs];
}

/* ---------------------------------------------------------------------- */

int SubModel::parse_coeffs(char **arg, int iarg, int narg)
{
  if (iarg + num_coeffs >= narg)
    error->all(FLERR, "Insufficient arguments provided for {} model", name);

  for (int i = 0; i < num_coeffs; i++) {
    // A few parameters (eg.g kt for tangential mindlin) allow null
    if (strcmp(arg[iarg+i+1], "NULL") == 0) coeffs[i] = -1;
    else coeffs[i] = utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
  }
  coeffs_to_local();

  return iarg + num_coeffs;
}

/* ----------------------------------------------------------------------
   mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double SubModel::mix_stiffnessE(double E1, double E2,
                                    double pois1, double pois2)
{
  double factor1 = (1 - pois1 * pois1) / E1;
  double factor2 = (1 - pois2 * pois2) / E2;
  return 1 / (factor1 + factor2);
}

/* ----------------------------------------------------------------------
   mixing of shear modulus (G)
------------------------------------------------------------------------ */

double SubModel::mix_stiffnessG(double E1, double E2,
                                    double pois1, double pois2)
{
  double factor1 = 2 * (2 - pois1) * (1 + pois1) / E1;
  double factor2 = 2 * (2 - pois2) * (1 + pois2) / E2;
  return 1 / (factor1 + factor2);
}

/* ----------------------------------------------------------------------
   mixing of everything else
------------------------------------------------------------------------- */

double SubModel::mix_geom(double val1, double val2)
{
  return sqrt(val1 * val2);
}
