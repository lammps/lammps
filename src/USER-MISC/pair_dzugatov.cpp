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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_dzugatov.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;



PairDzugatov::PairDzugatov(class LAMMPS *lmp) : Pair(lmp)
{
  write_data = 1;
}


PairDzugatov::~PairDzugatov()
{
  if( allocatd ) {
    
  }
}


void PairDzugatov::compute(int, int)
{

}

void PairDzugatov::settings(int, char **)
{}

void PairDzugatov::coeff(int, char **)
{}

double PairDzugatov::init_one(int, int)
{
  return 0.0;
}

void PairDzugatov::write_restart(FILE *)
{}

void PairDzugatov::read_restart(FILE *)
{}

void PairDzugatov::write_restart_settings(FILE *)
{}

void PairDzugatov::read_restart_settings(FILE *)
{}

void PairDzugatov::write_data(FILE *)
{}

void PairDzugatov::write_data_all(FILE *)
{}

double PairDzugatov::single(int, int, int, int, double, double, double, double &)
{
  return 0.0;
}

void PairDzugatov::*extract(const char *, int &)
{}
