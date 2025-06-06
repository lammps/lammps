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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#include "charge_solver.h"

using namespace LAMMPS_NS;

ChargeSolver::ChargeSolver()
{
  constraint = ChargeConstraint::NONE;
  mult_time = 0.;
}

/* ---------------------------------------------------------------------- */

ChargeSolver::~ChargeSolver() {}

/* ---------------------------------------------------------------------- */

void ChargeSolver::set_constraint(double qtotal)
{
  this->qtotal = qtotal;
  constraint = ChargeConstraint::SINGLE;
}
/* ---------------------------------------------------------------------- */

void ChargeSolver::set_constraint(std::vector<double> qtotal_group)
{
  this->qtotal_group = qtotal_group;
  constraint = ChargeConstraint::GROUP;
}

/* ---------------------------------------------------------------------- */

double ChargeSolver::get_mult_time()
{
  return mult_time;
}
