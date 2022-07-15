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
#include "contact_sub_model.h"
#include "pointers.h"
#include "utils.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

namespace Contact{

SubModel::SubModel()
{
  allocated = 0;
  size_history = 0;
  history_index = 0;
  allow_limit_damping = 1;
  beyond_contact = 0;
}

/* ---------------------------------------------------------------------- */

SubModel::~SubModel()
{
  if (allocated) delete [] coeffs;
}

/* ---------------------------------------------------------------------- */

void SubModel::allocate_coeffs()
{
  allocated = 1;
  coeffs = new double[num_coeffs];
}

/* ---------------------------------------------------------------------- */

void SubModel::parse_coeffs(char **arg, int iarg)
{
  for (int i = 0; i < num_coeffs; i++) {
    // A few parameters (eg.g kt for tangential mindlin) allow null
    if (strcmp(arg[iarg+i+1], "NULL") == 0) coeffs[i] = -1;
    else coeffs[i] = utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
  }
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void SubModel::write_restart(FILE *fp)
{
  fwrite(&model_name.length(),sizeof(int),1,fp);
  fwrite(model_name.data(),sizeof(char),model_name.length(),fp);
  fwrite(&num_coeffs,sizeof(int),1,fp);
  fwrite(coeffs,sizeof(int),num_coeffs,fp);
}

/* ---------------------------------------------------------------------- */

void SubModel::read_restart(FILE *fp, int num_char)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&num_coeffs,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(const_cast<char*>(model_name.data()), num_char, MPI_CHAR, world);
  allocate_coeffs();
}

/* ---------------------------------------------------------------------- */

void SubModel::read_restart(FILE *fp)
{
  int num_char;
  if (me == 0) {
    utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(&num_char, 1, MPI_INT, 0, world);
  read_restart(fp, num_char);
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
  double factor2 = 2 * (2 - pois2) * (1 + pois2) / E2)
  return 1 / (factor1 + factor2);
}

/* ----------------------------------------------------------------------
   mixing of everything else
------------------------------------------------------------------------- */

double SubModel::mix_geom(double val1, double val2)
{
  return sqrt(val1 * val2);
}


}