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
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#include "bond_mesocnt.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

static constexpr double A_CC = 1.421;

/* ---------------------------------------------------------------------- */

BondMesoCNT::BondMesoCNT(LAMMPS *_lmp) : BondHarmonic(_lmp)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

BondMesoCNT::~BondMesoCNT()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondMesoCNT::coeff(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "bond_coeff", error);

  // units, eV to energy unit conversion

  double ang = force->angstrom;
  double eunit;
  if (strcmp(update->unit_style, "real") == 0)
    eunit = 23.06054966;
  else if (strcmp(update->unit_style, "metal") == 0)
    eunit = 1.0;
  else if (strcmp(update->unit_style, "si") == 0)
    eunit = 1.6021765e-19;
  else if (strcmp(update->unit_style, "cgs") == 0)
    eunit = 1.6021765e-12;
  else if (strcmp(update->unit_style, "electron") == 0)
    eunit = 3.674932248e-2;
  else if (strcmp(update->unit_style, "micro") == 0)
    eunit = 1.6021765e-4;
  else if (strcmp(update->unit_style, "nano") == 0)
    eunit = 1.6021765e2;
  else
    error->all(FLERR, "Bond style mesocnt does not support {} units", update->unit_style);

  // set parameters

  double k_one, r0_one;
  if (strcmp(arg[1], "custom") == 0) {
    if (narg != 4) error->all(FLERR, "Incorrect number of args for 'custom' bond coefficients");
    k_one = utils::numeric(FLERR, arg[2], false, lmp);
    r0_one = utils::numeric(FLERR, arg[3], false, lmp);
  } else if (strcmp(arg[1], "C") == 0) {
    if (narg != 5)
      error->all(FLERR, "Incorrect number of args for 'C' preset in bond coefficients");
    int n = utils::inumeric(FLERR, arg[2], false, lmp);
    int m = utils::inumeric(FLERR, arg[3], false, lmp);
    r0_one = utils::numeric(FLERR, arg[4], false, lmp);

    double r_ang = sqrt(3.0 * (n * n + n * m + m * m)) * A_CC / MY_2PI;
    k_one = 0.5 * (86.64 + 100.56 * r_ang) * eunit / (ang * r0_one);
  } else
    error->all(FLERR, "Unknown {} preset in bond coefficients", arg[1]);

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Invalid bond type {}", arg[0]);
}
