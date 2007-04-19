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
#include "string.h"
#include "stdlib.h"
#include "compute_variable.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{INDEX,LOOP,EQUAL,WORLD,UNIVERSE,ULOOP,ATOM};  // also in variable.cpp

/* ---------------------------------------------------------------------- */

ComputeVariable::ComputeVariable(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute variable command");

  // store variable name

  int n = strlen(arg[3]) + 1;
  varname = new char[n];
  strcpy(varname,arg[3]);

  scalar_flag = 1;
  extensive = 0;
}

/* ---------------------------------------------------------------------- */

ComputeVariable::~ComputeVariable()
{
  delete [] varname;
}

/* ---------------------------------------------------------------------- */

void ComputeVariable::init()
{
  // check if variable exists

  int ivariable = input->variable->find(varname);
  if (ivariable < 0)
    error->all("Could not find compute variable name");
}

/* ---------------------------------------------------------------------- */

double ComputeVariable::compute_scalar()
{
  scalar = atof(input->variable->retrieve(varname));
  return scalar;
}
