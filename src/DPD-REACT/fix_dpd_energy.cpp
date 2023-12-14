// clang-format off
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

#include "fix_dpd_energy.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "pair_dpd_fdt_energy.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDPDenergy::FixDPDenergy(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3 ) error->all(FLERR,"Illegal fix dpd/energy command");

  pairDPDE = nullptr;
  pairDPDE = dynamic_cast<PairDPDfdtEnergy *>(force->pair_match("dpd/fdt/energy",1));
  if (pairDPDE == nullptr)
    pairDPDE = dynamic_cast<PairDPDfdtEnergy *>(force->pair_match("dpd/fdt/energy/kk",1));

  if (pairDPDE == nullptr)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy with fix dpd/energy");
  if (!(atom->dpd_flag))
    error->all(FLERR,"Must use atom_style dpd/fdt/energy with fix dpd/energy");
}

/* ---------------------------------------------------------------------- */

int FixDPDenergy::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixDPDenergy::initial_integrate(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *duCond = pairDPDE->duCond;
  double *duMech = pairDPDE->duMech;

  for (int i = 0; i < nlocal; i++) {
    uCond[i] += 0.5*update->dt*duCond[i];
    uMech[i] += 0.5*update->dt*duMech[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixDPDenergy::final_integrate()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *duCond = pairDPDE->duCond;
  double *duMech = pairDPDE->duMech;

  for (int i = 0; i < nlocal; i++) {
    uCond[i] += 0.5*update->dt*duCond[i];
    uMech[i] += 0.5*update->dt*duMech[i];
  }
}
