/* -*- c++ -*- ----------------------------------------------------------
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
#include "ldd_potential.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

LddPotential::LddPotential(class LAMMPS *lmp) : Pointers(lmp)
{
  allocated = n_coeffs = 0;
  coeffs = nullptr;
  potl_table.n_pts = 0;
  potl_table.dr = 0.0;
  potl_table.r = nullptr;
  potl_table.u = nullptr;
  potl_table.u2 = nullptr;
  potl_table.f = nullptr;
  potl_table.f2 = nullptr;
}

/* ---------------------------------------------------------------------- */

LddPotential::~LddPotential()
{
  allocated = 0;
}

/* ---------------------------------------------------------------------- */
void LddPotential::read_table_file(char *fnm, bool bspline)
{
  potl_table.n_pts = 0;
  FILE *fp = nullptr;
  char line[1000];
  if (comm->me == 0)
  {
  fp = utils::open_potential(fnm, lmp, nullptr);
  if (!fp)
    error->one(FLERR, "unable to open table file {} for reading: {}", fnm, utils::getsyserror());

  while (fgets(line, 1000, fp) != nullptr) { ++(potl_table.n_pts); }
  rewind(fp);
  }
  MPI_Bcast(&potl_table.n_pts, 1, MPI_INT, 0, world );

  memory->create(potl_table.r, potl_table.n_pts, "LDpotl:r");
  memory->create(potl_table.u, potl_table.n_pts, "LDpotl:u");
  memory->create(potl_table.f, potl_table.n_pts, "LDpotl:f");
  if (bspline) {
    memory->create(potl_table.u2, potl_table.n_pts, "LDpotl:u2");
    memory->create(potl_table.f2, potl_table.n_pts, "LDpotl:f2");
  }

  if (comm->me == 0)
  {
  for (int i = 0; i < potl_table.n_pts; ++i) {
    fgets(line, 1000, fp);
    float rt, ut, ft;
    int test_sscanf = sscanf(line, " %f %f %f ", &rt, &ut, &ft);
    if (test_sscanf != 3)
      error->one(FLERR, "unable to read rho u f from line {} in file {}", i, fnm);

    potl_table.r[i] = (double) rt;
    potl_table.u[i] = (double) ut;
    potl_table.f[i] = (double) ft;
  }
  fclose(fp);
  }
  MPI_Bcast(potl_table.r, potl_table.n_pts, MPI_DOUBLE, 0, world);
  MPI_Bcast(potl_table.u, potl_table.n_pts, MPI_DOUBLE, 0, world);
  MPI_Bcast(potl_table.f, potl_table.n_pts, MPI_DOUBLE, 0, world);

  potl_table.dr = potl_table.r[1] - potl_table.r[0];
  for (int i = 0; i < potl_table.n_pts - 1; ++i) {
    if (fabs(potl_table.r[i + 1] - potl_table.r[i] - potl_table.dr) > potl_table.dr / 100.0)
      error->all(FLERR, "in file {}\nr[1] - r[0] = {} - {} = {}\nr[{}] - r[{}] = {} - {} = {}", fnm,
                 potl_table.r[1], potl_table.r[0], potl_table.dr, i + 1, i, potl_table.r[i + 1],
                 potl_table.r[i], potl_table.r[i + 1] - potl_table.r[i]);
  }
}
/* ---------------------------------------------------------------------- */
int LddPotential::get_table_index(double r)
{
  if (r > potl_table.r[potl_table.n_pts - 1])
    error->all(FLERR, "local density = {} > {} = table hi", r, potl_table.r[potl_table.n_pts - 1]);

  if (r < potl_table.r[0])
    error->all(FLERR, "local density = {} < {} = table lo", r, potl_table.r[0]);

  return (int)floor((r - potl_table.r[0]) / potl_table.dr);
}

/* ---------------------------------------------------------------------- */
double LddPotential::calc_A_table(double r, int idx)
{
  return ((potl_table.r[idx + 1] - r) / potl_table.dr);
}
