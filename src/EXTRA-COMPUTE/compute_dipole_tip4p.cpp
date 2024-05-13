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

#include "compute_dipole_tip4p.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { MASSCENTER, GEOMCENTER };

/* ---------------------------------------------------------------------- */

ComputeDipoleTIP4P::ComputeDipoleTIP4P(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if ((narg < 3) || (narg > 4)) error->all(FLERR, "Illegal compute dipole/tip4p command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  extscalar = 0;
  extvector = 0;

  vector = new double[size_vector];
  vector[0] = vector[1] = vector[2] = 0.0;
  usecenter = MASSCENTER;

  if (narg == 4) {
    if (utils::strmatch(arg[3], "^geom"))
      usecenter = GEOMCENTER;
    else if (strcmp(arg[3], "mass") == 0)
      usecenter = MASSCENTER;
    else
      error->all(FLERR, "Illegal compute dipole/tip4p command");
  }
}

/* ---------------------------------------------------------------------- */

ComputeDipoleTIP4P::~ComputeDipoleTIP4P()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleTIP4P::init()
{
  if (!force->pair) error->all(FLERR, "Pair style must be defined for compute dipole/ti4p");

  int itmp;
  double *p_qdist = (double *) force->pair->extract("qdist", itmp);
  int *p_typeO = (int *) force->pair->extract("typeO", itmp);
  int *p_typeH = (int *) force->pair->extract("typeH", itmp);
  int *p_typeA = (int *) force->pair->extract("typeA", itmp);
  int *p_typeB = (int *) force->pair->extract("typeB", itmp);
  if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
    error->all(FLERR, "Pair style is incompatible with compute dipole/tip4p");
  typeO = *p_typeO;
  typeH = *p_typeH;
  int typeA = *p_typeA;
  int typeB = *p_typeB;

  if (!force->angle || !force->bond || !force->angle->setflag || !force->bond->setflag)
    error->all(FLERR, "Bond and angle potentials must be defined for compute dipole/tip4p");
  if ((typeA < 1) || (typeA > atom->nangletypes) || (force->angle->setflag[typeA] == 0))
    error->all(FLERR, "Bad TIP4P angle type for compute dipole/tip4p");
  if ((typeB < 1) || (typeB > atom->nbondtypes) || (force->bond->setflag[typeB] == 0))
    error->all(FLERR, "Bad TIP4P bond type for compute dipole/tip4p");
  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = *p_qdist / (cos(0.5 * theta) * blen);
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleTIP4P::compute_vector()
{
  invoked_vector = update->ntimestep;

  const auto x = atom->x;
  const auto mask = atom->mask;
  const auto type = atom->type;
  const auto image = atom->image;
  const auto mass = atom->mass;
  const auto rmass = atom->rmass;
  const auto q = atom->q;
  const auto mu = atom->mu;
  const auto nlocal = atom->nlocal;

  double dipole[3] = {0.0, 0.0, 0.0};
  double comproc[3] = {0.0, 0.0, 0.0};
  double com[3] = {0.0, 0.0, 0.0};
  double masstotal = 0.0;
  double chrgtotal = 0.0;
  double massproc = 0.0;
  double chrgproc = 0.0;
  double unwrap[3], xM[3];
  double *xi;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      double massone = 1.0;    // for usecenter == GEOMCENTER
      if (usecenter == MASSCENTER) {
        if (rmass)
          massone = rmass[i];
        else
          massone = mass[type[i]];
      }
      massproc += massone;
      if (atom->q_flag) chrgproc += q[i];
      domain->unmap(x[i], image[i], unwrap);
      comproc[0] += unwrap[0] * massone;
      comproc[1] += unwrap[1] * massone;
      comproc[2] += unwrap[2] * massone;
    }
  }
  MPI_Allreduce(&massproc, &masstotal, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&chrgproc, &chrgtotal, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(comproc, com, 3, MPI_DOUBLE, MPI_SUM, world);

  if (masstotal > 0.0) {
    com[0] /= masstotal;
    com[1] /= masstotal;
    com[2] /= masstotal;
  }

  // compute dipole moment

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == typeO) {
        find_M(i, xM);
        xi = xM;
      } else {
        xi = x[i];
      }
      domain->unmap(xi, image[i], unwrap);
      if (atom->q_flag) {
        dipole[0] += q[i] * unwrap[0];
        dipole[1] += q[i] * unwrap[1];
        dipole[2] += q[i] * unwrap[2];
      }
      if (atom->mu_flag) {
        dipole[0] += mu[i][0];
        dipole[1] += mu[i][1];
        dipole[2] += mu[i][2];
      }
    }
  }

  MPI_Allreduce(dipole, vector, 3, MPI_DOUBLE, MPI_SUM, world);

  // correct for position dependence with a net charged group
  vector[0] -= chrgtotal * com[0];
  vector[1] -= chrgtotal * com[1];
  vector[2] -= chrgtotal * com[2];
}

/* ---------------------------------------------------------------------- */

double ComputeDipoleTIP4P::compute_scalar()
{
  if (invoked_vector != update->ntimestep) compute_vector();

  invoked_scalar = update->ntimestep;
  scalar = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleTIP4P::find_M(int i, double *xM)
{
  double **x = atom->x;

  int iH1 = atom->map(atom->tag[i] + 1);
  int iH2 = atom->map(atom->tag[i] + 2);

  if ((iH1 == -1) || (iH2 == -1)) error->one(FLERR, "TIP4P hydrogen is missing");
  if ((atom->type[iH1] != typeH) || (atom->type[iH2] != typeH))
    error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

  // set iH1,iH2 to index of closest image to O

  iH1 = domain->closest_image(i, iH1);
  iH2 = domain->closest_image(i, iH2);

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];

  xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
}
