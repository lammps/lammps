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
   Package      FixPIMDBNVT
   Purpose      Path Integral Molecular Dynamics of Bosons with Nose-Hoover Thermostat
   Copyright    Hirshberg lab @ Tel Aviv University
   Authors      Ofir Blumer, Jacob Higer, Yotam Feldman (yotam.feldman at gmail.com), Barak Hirshberg (hirshb at tauex.tau.ac.il)

   Updated      Jan-06-2025
   Version      1.0
------------------------------------------------------------------------- */

#include "fix_pimd_nvt_bosonic.h"

#include "bosonic_exchange.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "universe.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::FixPIMDBNVT(LAMMPS *lmp, int narg, char **arg) : FixPIMDNVT(lmp, narg, arg)

{
  bosonic_exchange = new BosonicExchange(lmp, atom->nlocal, np, universe->me, true, true);
  virial = 0.0;
  prim = 0.0;
  spring_energy = 0.0;
  size_vector = 4;
  if (method != PIMD && method != NMPIMD) {
    error->universe_all(FLERR,
                        "Method not supported in fix pimdb/nvt; only methods PIMD and NMPIMD");
  }
  if (comm->nprocs != 1)
    error->universe_all(FLERR,
                        fmt::format("Fix {} only supports a single processor per bead", style));
}

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::~FixPIMDBNVT()
{
  delete bosonic_exchange;
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::pre_spring_force_estimators()
{
  FixPIMDNVT::pre_spring_force_estimators();
  spring_energy = bosonic_exchange->get_bead_spring_energy();
  prim = bosonic_exchange->prim_estimator();
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::prepare_coordinates()
{
  comm_exec(atom->x);
  double **x = atom->x;
  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];
  double ff = fbond * atom->mass[atom->type[0]];
  bosonic_exchange->prepare_with_coordinates(*x, xlast, xnext, beta, -ff);
}

void FixPIMDBNVT::spring_force()
{
  double **f = atom->f;

  bosonic_exchange->spring_force(f);
}

/* ---------------------------------------------------------------------- */

double FixPIMDBNVT::compute_vector(int n)
{
  if (0 <= n && n < 3) { return FixPIMDNVT::compute_vector(n); }

  if (n == 3) {
    return prim;
  } else {
    error->universe_all(FLERR, "Fix only has 4 outputs!");
  }

  return 0.0;
}
