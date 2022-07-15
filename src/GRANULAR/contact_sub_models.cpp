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
  material_prop_flag = 0;
  size_history = 0;
  history_index = 0;
}

SubModel::~SubModel()
{
  if (allocated) delete [] coeffs;
}

void SubModel::allocate_coeffs()
{
  allocated = 1;
  coeffs = new double[num_coeffs];
}

void SubModel::parse_coeffs(char **arg, int iarg)
{
  for (int i = 0; i < num_coeffs; i++) {
    coeffs[i] = utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
  }
  coeffs_to_local();
}

void SubModel::write_restart(FILE *fp)
{
  fwrite(&model_name.length(),sizeof(int),1,fp);
  fwrite(model_name.data(),sizeof(char),model_name.length(),fp);
  fwrite(&num_coeffs,sizeof(int),1,fp);
  fwrite(coeffs,sizeof(int),num_coeffs,fp);
}

void SubModel::read_restart(FILE *fp, int num_char)
{
  if (comm->me == 0){
    utils::sfread(FLERR,&num_coeffs,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(const_cast<char*>(model_name.data()), num_char, MPI_CHAR, world);
  allocate_coeffs();
}

void SubModel::read_restart(FILE *fp)
{
  int num_char;
  if (me == 0){
    utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(&num_char, 1, MPI_INT, 0, world);
  read_restart(fp, num_char);
}




}

