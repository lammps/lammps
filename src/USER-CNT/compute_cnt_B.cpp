/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "compute_cnt_B.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "bond.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include <typeinfo>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCNT_B::ComputeCNT_B(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  buckling(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute cnt/B command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peatomflag = 1;
  timeflag = 1;
  comm_reverse = 0;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCNT_B::~ComputeCNT_B()
{
  memory->destroy(buckling);
}

/* ---------------------------------------------------------------------- */

void ComputeCNT_B::compute_peratom()
{
  int i;

  // grow local buckling array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    memory->destroy(buckling);
    nmax = atom->nmax;
    memory->create(buckling,nmax,"cnt_B:buckling");
    vector_atom = buckling;
  }

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) buckling[i] = atom->buckling[i];
}

/* ---------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCNT_B::memory_usage()
{
  double bytes = nmax * sizeof(int);
  return bytes;
}
