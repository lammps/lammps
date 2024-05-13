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

#include "body.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Body::Body(LAMMPS *lmp, int /*narg*/, char **arg) : Pointers(lmp)
{
  style = utils::strdup(arg[0]);
  icp = nullptr;
  dcp = nullptr;
}

/* ---------------------------------------------------------------------- */

Body::~Body()
{
  delete[] style;
}
