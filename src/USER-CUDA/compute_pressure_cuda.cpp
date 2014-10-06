/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include <cstring>
#include <cstdlib>
#include "compute_pressure_cuda.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"
#include "user_cuda.h"

using namespace LAMMPS_NS;

enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

/* ---------------------------------------------------------------------- */

ComputePressureCuda::ComputePressureCuda(LAMMPS *lmp, int narg, char **arg) :
  ComputePressure(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");
  cudable = 1;

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  int n = strlen(arg[3]) + 1;
  char* id_temp = new char[n];
  strcpy(id_temp,arg[3]);

  int icompute = modify->find_compute(id_temp);
  delete [] id_temp;
  if (modify->compute[icompute]->cudable == 0)
  {
    error->warning(FLERR,"Compute pressure/cuda temperature ID is not cudable! Try a temp/cuda style.");
    cudable = 0;
  }

}

double ComputePressureCuda::compute_scalar()
{
  if(not temperature->cudable && cuda->finished_setup) cuda->downloadAll();
  return ComputePressure::compute_scalar();
}

void ComputePressureCuda::compute_vector()
{
  if(not temperature->cudable && cuda->finished_setup) cuda->downloadAll();
  ComputePressure::compute_vector();
}
