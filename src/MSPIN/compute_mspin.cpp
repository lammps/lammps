/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#include "compute_mspin.h"
#include "error.h"
#include "fix_mspin_nh.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include <cstdio>
#include <cstring>

using namespace LAMMPS_NS;
using namespace std;

// @todo: change this into rigid/mspin/energy
ComputeMspin::ComputeMspin(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg), rfix(NULL)
{
  if (narg != 4) error->all(FLERR, "Illegal compute mspin command.");

  vector_flag = 1;
  size_vector = 2;

  // are the quantities extensive or intensive
  extvector = 0;
  extlist = new int[size_vector];
  extlist[0] = 0;
  extlist[1] = 0;

  int n = strlen(arg[3]) + 1;
  rfix = new char[n];
  strcpy(rfix, arg[3]);

  memory->create(vector, 2, "compute/mspin:vector");
}

ComputeMspin::~ComputeMspin()
{
  delete[] rfix;
  delete[] extlist;
  memory->destroy(vector);
}

void ComputeMspin::init()
{
  irfix = modify->find_fix(rfix);
  if (irfix < 0) error->all(FLERR, "Fix ID for compute mspin does not exist");

  if (strncmp(modify->fix[irfix]->style, "rigid/mspin", 11))
    error->all(FLERR, "Compute mspin with non-mspin fix-ID");
}

void ComputeMspin::compute_vector()
{
  double zeeman, dipolar;

  invoked_vector = update->ntimestep;

  if (strncmp(modify->fix[irfix]->style, "rigid/mspin", 11) == 0) {
    zeeman = ((FixMspinNH *) modify->fix[irfix])->extract_zeeman_pe();
    dipolar = ((FixMspinNH *) modify->fix[irfix])->extract_dipolar_pe();
  }

  vector[0] = dipolar;
  vector[1] = zeeman;
}