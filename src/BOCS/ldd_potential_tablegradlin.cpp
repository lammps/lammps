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
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_tablegradlin.h"

#include "error.h"
#include "memory.h"
#include "utils.h"

using namespace LAMMPS_NS;

LddPotentialTableGradLin::LddPotentialTableGradLin(class LAMMPS *lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
}

LddPotentialTableGradLin::~LddPotentialTableGradLin()
{
  if (allocated == 1) {
    memory->destroy(coeffs);
    memory->destroy(potl_table.r);
    memory->destroy(potl_table.u);
    memory->destroy(potl_table.f);
  }
  allocated = 0;
}

void LddPotentialTableGradLin::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_potential:coeffs");
  allocated = 1;
}

void LddPotentialTableGradLin::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg <= ipt + 1)
    error->all(FLERR, "Missing filename following the ldd table/gradlin keyword");

  read_table_file(arg[ipt + 2], false);
}

double LddPotentialTableGradLin::u(double rho)
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

double LddPotentialTableGradLin::f(double rho)
{
  if (rho == potl_table.r[potl_table.n_pts - 1]) {
    return 1 / potl_table.dr *
        potl_table.u[potl_table.n_pts - 1];    // potl_table.f[potl_table.n_pts-1];
  }
  int idx = get_table_index(rho);
  //  double A = calc_A_table(rho, idx);
  //  double B = 1.0 - A;
  if (rho == potl_table.r[0]) { return potl_table.f[0]; }
  double Bp =
      -1 / potl_table.dr;    // In BOCS we calculate the deriv, here I need the negative deriv
  double Ap = -Bp;
  return (Ap * potl_table.u[idx] +
          Bp * potl_table.u[idx + 1]);    // MCL Gotta use a delta interp for f if linear for u
  // I think this amounts to using a delta on the fwd dif at idx i
}
