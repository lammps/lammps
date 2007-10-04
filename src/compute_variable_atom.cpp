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

#include "string.h"
#include "compute_variable_atom.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVariableAtom::ComputeVariableAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute variable/atom command");

  // store variable name

  int n = strlen(arg[3]) + 1;
  varname = new char[n];
  strcpy(varname,arg[3]);

  peratom_flag = 1;
  size_peratom = 0;

  nmax = 0;
  result = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeVariableAtom::~ComputeVariableAtom()
{
  delete [] varname;
  memory->sfree(result);
}

/* ---------------------------------------------------------------------- */

void ComputeVariableAtom::init()
{
  // set ivariable used by this compute and check if it exists

  ivariable = input->variable->find(varname);
  if (ivariable < 0)
    error->all("Could not find compute variable name");
}

/* ---------------------------------------------------------------------- */

void ComputeVariableAtom::compute_peratom()
{
  // grow result array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(result);
    nmax = atom->nmax;
    result = (double *) memory->smalloc(nmax*sizeof(double),
					"compute/variable/atom:result");
    scalar_atom = result;
  }

  // parse variable once to create parse tree
  // evaluate tree for all atoms, will be zero for atoms not in group
  // free parse tree memory stored by Variable

  input->variable->build_parse_tree(ivariable);
  input->variable->evaluate_parse_tree(igroup,result);
  input->variable->free_parse_tree();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeVariableAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
