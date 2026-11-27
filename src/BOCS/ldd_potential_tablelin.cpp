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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_tablelin.h"

#include "error.h"
#include "memory.h"
#include "utils.h"

using namespace LAMMPS_NS;

LddPotentialTableLin::LddPotentialTableLin(class LAMMPS *lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
}

LddPotentialTableLin::~LddPotentialTableLin()
{
  if (allocated == 1) {
    memory->destroy(coeffs);
    memory->destroy(potl_table.r);
    memory->destroy(potl_table.u);
    memory->destroy(potl_table.f);
  }
  allocated = 0;
}

void LddPotentialTableLin::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_potential:coeffs");
  allocated = 1;
}

void LddPotentialTableLin::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg <= ipt + 1)
    error->all(FLERR, "Missing filename following the ldd table/lin keyword");

  read_table_file(arg[ipt + 2], false);
}

double LddPotentialTableLin::u(double rho)
{
  // Handle this case separately
  if (rho == potl_table.r[potl_table.n_pts - 1]) { return potl_table.u[potl_table.n_pts - 1]; }
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);
  double B = 1.0 - A;
  // If we didn't handle the first case separately, we'd try to access
  // potl_table.u[idx+1] here and it wouldn't work
  return (A * potl_table.u[idx] + B * potl_table.u[idx + 1]);
}

double LddPotentialTableLin::f(double rho)
{
  if (rho == potl_table.r[potl_table.n_pts - 1]) return potl_table.f[potl_table.n_pts - 1];
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);
  double B = 1.0 - A;
  return (A * potl_table.f[idx] + B * potl_table.f[idx + 1]);
}
