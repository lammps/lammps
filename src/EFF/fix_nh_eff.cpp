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

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "fix_nh_eff.h"

#include "atom.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixNHEff::FixNHEff(LAMMPS *lmp, int narg, char **arg) : FixNH(lmp, narg, arg)
{
  if (!atom->electron_flag)
    error->all(FLERR,"Fix {} requires atom style electron", style);
}

/* ----------------------------------------------------------------------
   perform half-step update of electron radial velocities
-----------------------------------------------------------------------*/

void FixNHEff::nve_v()
{
  // standard nve_v velocity update

  FixNH::nve_v();

  double *erforce = atom->erforce;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  double mefactor = domain->dimension/4.0;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfm;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (abs(spin[i])==1) {
        dtfm = dtf / mass[type[i]];
        ervel[i] += dtfm * erforce[i] / mefactor;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of electron radii
-----------------------------------------------------------------------*/

void FixNHEff::nve_x()
{
  // standard nve_x position update

  FixNH::nve_x();

  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *type = atom->type;
  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double mefactor = domain->dimension/4.0;
  double dtfm;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      if (abs(spin[i])==1) {
        ervel[i] += dtfm * erforce[i] / mefactor;
        eradius[i] += dtv * ervel[i];
      }
    }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of electron radial velocities
-----------------------------------------------------------------------*/

void FixNHEff::nh_v_temp()
{
  // standard nh_v_temp velocity scaling

  FixNH::nh_v_temp();

  double *ervel = atom->ervel;
  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (abs(spin[i])==1) ervel[i] *= factor_eta;
}
