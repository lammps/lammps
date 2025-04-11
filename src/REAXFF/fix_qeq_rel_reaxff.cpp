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
   Contributing authors:
      Navraj S Lalli, Imperial College London (navrajsinghlalli@gmail.com)

------------------------------------------------------------------------- */

#include "fix_qeq_rel_reaxff.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_efield.h"
#include "force.h"
#include "neighbor.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEqRelReaxFF::FixQEqRelReaxFF(LAMMPS *lmp, int narg, char **arg) : FixQtpieReaxFF(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixQEqRelReaxFF::calc_chi_eff()
{
  memset(&chi_eff[0], 0, atom->nmax * sizeof(double));

  const auto x = (const double *const *) atom->x;
  const int *type = atom->type;

  double dx, dy, dz, dist_sq, overlap, sum_n, sum_d, chia, phia, phib;
  int i, j;

  // check ghost atoms are stored up to the distance cutoff for overlap integrals
  const double comm_cutoff = MAX(neighbor->cutneighmax, comm->cutghostuser);
  if (comm_cutoff*comm_cutoff < dist_cutoff_sq) {
    error->all(FLERR, Error::NOLASTLINE,
               "Comm cutoff {} is smaller than distance cutoff {} for overlap integrals in fix {}. "
               "Increase accordingly using comm_modify cutoff",
               comm_cutoff, sqrt(dist_cutoff_sq), style);
  }

  // efield energy is in real units of kcal/mol, factor needed for conversion to eV
  const double qe2f = force->qe2f;
  const double factor = 1.0 / qe2f;

  if (efield) {
    if (efield->varflag != FixEfield::CONSTANT) efield->update_efield_variables();

    // compute chi_eff for each local atom
    for (i = 0; i < nn; i++) {
      chia = chi[type[i]];
      if (efield->varflag != FixEfield::ATOM) {
        phia = -factor * (x[i][0] * efield->ex + x[i][1] * efield->ey + x[i][2] * efield->ez);
      } else {    // atom-style potential from FixEfield
        phia = efield->efield[i][3];
      }

      sum_n = 0.0;
      sum_d = 0.0;

      for (j = 0; j < nt; j++) {
          dx = x[i][0] - x[j][0];
          dy = x[i][1] - x[j][1];
          dz = x[i][2] - x[j][2];
          dist_sq = (dx*dx + dy*dy + dz*dz);

        if (dist_sq < dist_cutoff_sq) {

          // overlap integral of two normalised 1s Gaussian type orbitals
          overlap = prefactor[type[i]][type[j]] * exp(-expfactor[type[i]][type[j]] * dist_sq);

          if (efield->varflag != FixEfield::ATOM) {
            phib = -factor * (x[j][0] * efield->ex + x[j][1] * efield->ey + x[j][2] * efield->ez);
          } else {    // atom-style potential from FixEfield
            phib = efield->efield[j][3];
          }
          sum_n += phib * overlap;
          sum_d += overlap;
        }
      }
      if (sum_d != 0.0)
        chi_eff[i] = chia + scale * (phia - sum_n / sum_d);
      else
        chi_eff[i] = chia;
    }
  } else {
    for (i = 0; i < nn; i++) { chi_eff[i] = chi[type[i]]; }
  }
}
