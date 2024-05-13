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

#include "compute_pe.h"

#include "angle.h"
#include "atom.h"
#include "atom_masks.h"
#include "bond.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePE::ComputePE(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal compute pe command");
  if (igroup) error->all(FLERR, "Compute pe must use group all");

  scalar_flag = 1;
  extscalar = 1;
  peflag = 1;
  timeflag = 1;

  if (narg == 3) {
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = 1;
    fixflag = 1;
  } else {
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = 0;
    fixflag = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "pair") == 0)
        pairflag = 1;
      else if (strcmp(arg[iarg], "bond") == 0)
        bondflag = 1;
      else if (strcmp(arg[iarg], "angle") == 0)
        angleflag = 1;
      else if (strcmp(arg[iarg], "dihedral") == 0)
        dihedralflag = 1;
      else if (strcmp(arg[iarg], "improper") == 0)
        improperflag = 1;
      else if (strcmp(arg[iarg], "kspace") == 0)
        kspaceflag = 1;
      else if (strcmp(arg[iarg], "fix") == 0)
        fixflag = 1;
      else
        error->all(FLERR, "Illegal compute pe command");
      iarg++;
    }
  }

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

double ComputePE::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR, "Energy was not tallied on needed timestep");

  double one = 0.0;
  if (pairflag && force->pair) one += force->pair->eng_vdwl + force->pair->eng_coul;

  if (atom->molecular != Atom::ATOMIC) {
    if (bondflag && force->bond) one += force->bond->energy;
    if (angleflag && force->angle) one += force->angle->energy;
    if (dihedralflag && force->dihedral) one += force->dihedral->energy;
    if (improperflag && force->improper) one += force->improper->energy;
  }

  MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  if (kspaceflag && force->kspace) scalar += force->kspace->energy;

  if (pairflag && force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    scalar += force->pair->etail / volume;
  }

  if (fixflag && modify->n_energy_global) scalar += modify->energy_global();

  return scalar;
}
