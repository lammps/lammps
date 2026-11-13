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
----------------------------------
   Contributing author: Joel Clemmer (SNL)

   Ewald summation for Newtonian gravity in fully periodic systems.

   Point masses interact like charges of a one-component plasma: the
   per-atom mass is the source of the field, the prefactor is -G, and
   the uniform background (the mean density that is subtracted in
   comoving cosmological simulations) is the standard correction for
   non-neutral systems.
   Reference: L. Hernquist, F. R. Bouchet, Y. Suto, ApJS 75, 231 (1991).
------------------------------------------------------------------------- */

#include "ewald_grav.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

EwaldGrav::EwaldGrav(LAMMPS *lmp) : Ewald(lmp), gconst(0.0), mass_src(nullptr), nmax_src(0)
{
  require_q = 0;
  group_group_enable = 0;
}

/* ---------------------------------------------------------------------- */

EwaldGrav::~EwaldGrav()
{
  memory->destroy(mass_src);
}

/* ---------------------------------------------------------------------- */

void EwaldGrav::init()
{
  // the gravitational constant is set by the pair style, so that the
  // short-range and long-range parts of the interaction always agree

  if (force->pair == nullptr)
    error->all(FLERR, "KSpace style {} requires a pair style", force->kspace_style);
  int itmp;
  auto *p_gconst = (double *) force->pair->extract("G_grav", itmp);
  if (p_gconst == nullptr)
    error->all(FLERR, "KSpace style {} requires a gravity pair style like grav/long",
               force->kspace_style);
  gconst = *p_gconst;
  if (gconst <= 0.0)
    error->all(FLERR, "KSpace style {} requires a positive gravitational constant",
               force->kspace_style);

  if (slabflag)
    error->all(FLERR, "KSpace style {} does not support the slab correction",
               force->kspace_style);

  // base class computes the parameters using this->qsum_qsq() and
  // this->two_charge(), which use the masses and G

  Ewald::init();

  // gravity is attractive: masses act like charges with prefactor -G

  qqrd2e = -gconst;
}

/* ---------------------------------------------------------------------- */

void EwaldGrav::compute(int eflag, int vflag)
{
  update_masses();
  Ewald::compute(eflag, vflag);
  if (qsqsum == 0.0) return;

  // virial of the uniform background, E_bg = G pi M^2 / (2 g_ewald^2 V), which
  // Ewald::compute() adds to the energy but not to the virial.  Not in
  // Hernquist et al. (1991), which only uses forces.  Akin to virial_self for
  // TERM_COUL in EwaldDisp (global and per-atom), but with the opposite sign;
  // the sign used here matches a finite-difference pressure.  This may be an
  // error in EwaldDisp, although relatively unimportant there for nearly
  // charge neutral systems.

  const double bgfac = -qqrd2e * scale * MY_PI2 / (g_ewald * g_ewald * volume);

  if (vflag_global) {
    const double ebg = bgfac * qsum * qsum;
    virial[0] += ebg;
    virial[1] += ebg;
    virial[2] += ebg;
  }

  if (vflag_atom) {
    const int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      const double ebg = bgfac * mass_src[i] * qsum;
      vatom[i][0] += ebg;
      vatom[i][1] += ebg;
      vatom[i][2] += ebg;
    }
  }
}

/* ----------------------------------------------------------------------
   total mass and sum of squared masses
   q2 is used for the accuracy estimates and therefore must be positive
   a net mass is the normal case, so there is no neutrality check
------------------------------------------------------------------------- */

void EwaldGrav::qsum_qsq(int warning_flag)
{
  update_masses();

  const int nlocal = atom->nlocal;
  double msum_local = 0.0, msqsum_local = 0.0;
  for (int i = 0; i < nlocal; i++) {
    msum_local += mass_src[i];
    msqsum_local += mass_src[i] * mass_src[i];
  }

  MPI_Allreduce(&msum_local, &qsum, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&msqsum_local, &qsqsum, 1, MPI_DOUBLE, MPI_SUM, world);

  if ((qsqsum == 0.0) && (comm->me == 0) && warn_nocharge && warning_flag) {
    error->warning(FLERR, "Using kspace style {} on system with no mass", force->kspace_style);
    warn_nocharge = 0;
  }

  q2 = qsqsum * gconst;
}

/* ----------------------------------------------------------------------
   force between two unit masses at unit distance
   reference for the relative accuracy
------------------------------------------------------------------------- */

void EwaldGrav::two_charge()
{
  two_charge_force = gconst;
}

/* ----------------------------------------------------------------------
   copy per-atom masses of owned atoms into the source array
------------------------------------------------------------------------- */

void EwaldGrav::update_masses()
{
  if (atom->nmax > nmax_src) {
    memory->destroy(mass_src);
    nmax_src = atom->nmax;
    memory->create(mass_src, nmax_src, "ewald/grav:mass_src");
  }

  const int nlocal = atom->nlocal;
  if (atom->rmass) {
    const double *const rmass = atom->rmass;
    for (int i = 0; i < nlocal; i++) mass_src[i] = rmass[i];
  } else {
    const double *const mass = atom->mass;
    const int *const type = atom->type;
    for (int i = 0; i < nlocal; i++) mass_src[i] = mass[type[i]];
  }
}

/* ---------------------------------------------------------------------- */

double EwaldGrav::memory_usage()
{
  return Ewald::memory_usage() + (double) nmax_src * sizeof(double);
}
