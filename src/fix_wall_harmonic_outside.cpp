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

#include "fix_wall_harmonic_outside.h"
#include "atom.h"
#include "comm.h"
#include "citeme.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

static const char cite_fix_wall_harmonic_outside[] =
    "fix wall/harmonic/outside command: https://doi.org/10.1021/acs.langmuir.4c04293\n\n"
    "@Article{Barraud2025,\n"
    " author  = {Barraud, Eddy and Dalmazzone, Christine and Mouret, "
    "Aurelie and De Bruin, Theodorus and Creton, Benoit and Pasquier, David "
    "and Lachet, Veronique and Nieto-Draghi, Carlos},\n"
    " title        = {A Coarse-Grained Model Describing the Critical Micelle Concentration"
    " of Perfluoroalkyl Surfactants in Ionic Aqueous Phase},\n"
    " journal = {Langmuir},\n"
    " year         = {2025},\n"
    " number       = {11}, \n"
    " pages        = {7272--7282}, \n"
    " volume       = {41}, \n"
    "}\n\n";

FixWallHarmonicOutside::FixWallHarmonicOutside(LAMMPS *lmp, int narg, char **arg) : FixWall(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_wall_harmonic_outside);
  /* Raise error if using the old name */
  if ((strcmp(style, "wall/harmonic/returned") == 0) && (comm->me == 0))
    error->warning(FLERR,"Fix wall/harmonic/returned was renamed to fix wall/harmonic/outside. Please update your input");
  dynamic_group_allow = 1;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   recalling force applied if outside the control volume
   and within the interaction cutoff
   m = index of wall coeffs
   which = 0,1,..,5 (xlo,xhi,ylo,yhi,zlo,zhi)
   coord = wall coordinate on the dim
   dim = 0,1,2 (x,y,z)
   side = -1,1 (low, high)
   if side is the low boundary,
   no error if any particle is on or above the wall
------------------------------------------------------------------------- */


void FixWallHarmonicOutside::wall_particle(int m, int which, double coord)
{
  double dr, fwall;
  double vn;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  // iterate through the atoms owned by the proc
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      // calculate the distance (dr) of each atom from the wall
      if (side < 0)
        dr = coord - x[i][dim];
      else
        dr = x[i][dim] - coord;
      if (dr >= cutoff[m]) continue; // no force if above the interaction cutoff
      if (dr <= 0.0) {
        /* No force if the particle is inside the control volume */
        continue;
      }
      fwall = side * 2.0 * epsilon[m] * dr; // calculate the simple harmonic force
      f[i][dim] -= fwall; // apply the force over the atom in the same dimension as the wall
      ewall[0] += epsilon[m] * dr * dr; // sum the energies of the walls for record
      ewall[m + 1] += fwall; // sum the forces of the wall for record

      if (evflag) {
        if (side < 0)
          vn = -fwall * dr;
        else
          vn = fwall * dr;
        v_tally(dim, i, vn);
      }
    }

}
