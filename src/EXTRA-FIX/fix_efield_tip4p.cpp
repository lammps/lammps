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

#include "fix_efield_tip4p.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEfieldTIP4P::FixEfieldTIP4P(LAMMPS *_lmp, int narg, char **arg) : FixEfield(_lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixEfieldTIP4P::init()
{
  FixEfield::init();

  if (atom->tag_enable == 0) error->all(FLERR, "Fix efield/tip4p requires atom IDs");
  if (!atom->q_flag) error->all(FLERR, "Fix efield/tip4p requires atom attribute q");
  if (!force->pair) error->all(FLERR, "A TIP4P pair style must be defined fix efield/tip4p");
  if (pstr) error->all(FLERR, "Fix efield/tip4p does not support the potential keyword");

  int itmp;
  double *p_qdist = (double *) force->pair->extract("qdist", itmp);
  int *p_typeO = (int *) force->pair->extract("typeO", itmp);
  int *p_typeH = (int *) force->pair->extract("typeH", itmp);
  int *p_typeA = (int *) force->pair->extract("typeA", itmp);
  int *p_typeB = (int *) force->pair->extract("typeB", itmp);
  if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
    error->all(FLERR, "Pair style is incompatible with compute {}", style);
  typeO = *p_typeO;
  typeH = *p_typeH;
  int typeA = *p_typeA;
  int typeB = *p_typeB;

  if (!force->angle || !force->bond || !force->angle->setflag || !force->bond->setflag)
    error->all(FLERR, "Bond and angle potentials must be defined for fix efield/tip4p");
  if ((typeA < 1) || (typeA > atom->nangletypes) || (force->angle->setflag[typeA] == 0))
    error->all(FLERR, "Bad TIP4P angle type for fix efield/tip4p");
  if ((typeB < 1) || (typeB > atom->nbondtypes) || (force->bond->setflag[typeB] == 0))
    error->all(FLERR, "Bad TIP4P bond type for fix efield/tip4p");

  // determine alpha parameter

  const double theta = force->angle->equilibrium_angle(typeA);
  const double blen = force->bond->equilibrium_distance(typeB);
  alpha = *p_qdist / (cos(0.5 * theta) * blen);
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixEfieldTIP4P::post_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // virial setup

  v_init(vflag);

  // reallocate efield array if necessary

  if ((varflag == ATOM) && (atom->nmax > maxatom)) {
    maxatom = atom->nmax;
    memory->destroy(efield);
    memory->create(efield, maxatom, 4, "efield:efield");
  }

  // update region if necessary

  if (region) region->prematch();

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double **x = atom->x;
  double fx, fy, fz, xM[3];
  double v[6], unwrap[3];
  int iO, iH1, iH2;

  // constant efield

  if (varflag == CONSTANT) {

    // charge interactions
    // force = qE, potential energy = F dot x in unwrapped coords

    if (qflag) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {

          // process *all* atoms belonging to a TIP4P molecule since we need
          // to consider molecules where parts of the molecule are ghost atoms

          if ((type[i] == typeO) || (type[i] == typeH)) {

            if (type[i] == typeO) {
              iO = i;
              iH1 = atom->map(tag[i] + 1);
              iH2 = atom->map(tag[i] + 2);
            } else {

              // set indices for first or second hydrogen

              iO = atom->map(tag[i] - 1);
              if ((iO != -1) && (type[iO] == typeO)) {
                iH1 = i;
                iH2 = atom->map(tag[i] + 1);
              } else {
                iO = atom->map(tag[i] - 2);
                iH1 = atom->map(tag[i] - 1);
                iH2 = i;
              }
            }
            if ((iH1 == -1) || (iH2 == -1)) error->one(FLERR, "TIP4P hydrogen is missing");
            if (iO == -1) error->one(FLERR, "TIP4P oxygen is missing");
            if ((type[iH1] != typeH) || (type[iH2] != typeH))
              error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
            if (type[iO] != typeO) error->one(FLERR, "TIP4P oxygen has incorrect atom type");
            iH1 = domain->closest_image(i, iH1);
            iH2 = domain->closest_image(i, iH2);

            find_M(x[iO], x[iH1], x[iH2], xM);

            // M site contributions

            if (!region || region->match(xM[0], xM[1], xM[2])) {

              fx = q[iO] * ex;
              fy = q[iO] * ey;
              fz = q[iO] * ez;

              // distribute and apply forces

              if (i == iO) {
                f[iO][0] += fx * (1.0 - alpha);
                f[iO][1] += fy * (1.0 - alpha);
                f[iO][2] += fz * (1.0 - alpha);
                if (iH1 < nlocal) {
                  f[iH1][0] += 0.5 * alpha * fx;
                  f[iH1][1] += 0.5 * alpha * fy;
                  f[iH1][2] += 0.5 * alpha * fz;
                }
                if (iH2 < nlocal) {
                  f[iH2][0] += 0.5 * alpha * fx;
                  f[iH2][1] += 0.5 * alpha * fy;
                  f[iH2][2] += 0.5 * alpha * fz;
                }
              } else {
                if (iO >= nlocal) {
                  if (i == iH1) {
                    f[iH1][0] += 0.5 * alpha * fx;
                    f[iH1][1] += 0.5 * alpha * fy;
                    f[iH1][2] += 0.5 * alpha * fz;
                  }
                  if (i == iH2) {
                    f[iH2][0] += 0.5 * alpha * fx;
                    f[iH2][1] += 0.5 * alpha * fy;
                    f[iH2][2] += 0.5 * alpha * fz;
                  }
                }
              }

              if (i == iO) {
                domain->unmap(xM, image[iO], unwrap);
                fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contribution with Oxygen

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iO, v);
                }
              }
            }

            // H1 site contributions

            if (!region || region->match(x[iH1][0], x[iH1][1], x[iH1][2])) {

              fx = q[iH1] * ex;
              fy = q[iH1] * ey;
              fz = q[iH1] * ez;

              if (i == iH1) {
                f[iH1][0] += fx;
                f[iH1][1] += fy;
                f[iH1][2] += fz;

                // tally global force

                domain->unmap(x[iH1], image[iH1], unwrap);
                fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contributions

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iH1, v);
                }
              }
            }

            // H2 site contributions

            if (!region || region->match(x[iH2][0], x[iH2][1], x[iH2][2])) {

              fx = q[iH2] * ex;
              fy = q[iH2] * ey;
              fz = q[iH2] * ez;

              if (i == iH2) {
                f[iH2][0] += fx;
                f[iH2][1] += fy;
                f[iH2][2] += fz;

                // tally global force

                domain->unmap(x[iH2], image[iH2], unwrap);
                fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contributions

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iH2, v);
                }
              }
            }

            // non-TIP4P atoms

          } else {

            if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
            fx = q[i] * ex;
            fy = q[i] * ey;
            fz = q[i] * ez;
            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;

            domain->unmap(x[i], image[i], unwrap);
            fsum[0] -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
            fsum[1] += fx;
            fsum[2] += fy;
            fsum[3] += fz;
            if (evflag) {
              v[0] = fx * unwrap[0];
              v[1] = fy * unwrap[1];
              v[2] = fz * unwrap[2];
              v[3] = fx * unwrap[1];
              v[4] = fx * unwrap[2];
              v[5] = fy * unwrap[2];
              v_tally(i, v);
            }
          }
        }
      }
    }

    // dipole interactions, no special TIP4P treatment needed
    // no force, torque = mu cross E, potential energy = -mu dot E

    if (muflag) {
      double **mu = atom->mu;
      double **t = atom->torque;
      double tx, ty, tz;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
          tx = ez * mu[i][1] - ey * mu[i][2];
          ty = ex * mu[i][2] - ez * mu[i][0];
          tz = ey * mu[i][0] - ex * mu[i][1];
          t[i][0] += tx;
          t[i][1] += ty;
          t[i][2] += tz;
          fsum[0] -= mu[i][0] * ex + mu[i][1] * ey + mu[i][2] * ez;
        }
    }

  } else {

    // variable efield, wrap with clear/add
    // potential energy = evar if defined, else 0.0

    modify->clearstep_compute();

    if (xstyle == EQUAL) {
      ex = qe2f * input->variable->compute_equal(xvar);
    } else if (xstyle == ATOM) {
      input->variable->compute_atom(xvar, igroup, &efield[0][0], 4, 0);
    }
    if (ystyle == EQUAL) {
      ey = qe2f * input->variable->compute_equal(yvar);
    } else if (ystyle == ATOM) {
      input->variable->compute_atom(yvar, igroup, &efield[0][1], 4, 0);
    }
    if (zstyle == EQUAL) {
      ez = qe2f * input->variable->compute_equal(zvar);
    } else if (zstyle == ATOM) {
      input->variable->compute_atom(zvar, igroup, &efield[0][2], 4, 0);
    }
    if (estyle == ATOM) input->variable->compute_atom(evar, igroup, &efield[0][3], 4, 0);

    modify->addstep_compute(update->ntimestep + 1);

    // charge interactions
    // force = qE

    if (qflag) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {

          // process *all* atoms belonging to a TIP4P molecule since we need
          // to consider molecules where parts of the molecule are ghost atoms

          if ((type[i] == typeO) || (type[i] == typeH)) {

            if (type[i] == typeO) {
              iO = i;
              iH1 = atom->map(tag[i] + 1);
              iH2 = atom->map(tag[i] + 2);
            } else {

              // set indices for first or second hydrogen

              iO = atom->map(tag[i] - 1);
              if ((iO != -1) && (type[iO] == typeO)) {
                iH1 = i;
                iH2 = atom->map(tag[i] + 1);
              } else {
                iO = atom->map(tag[i] - 2);
                iH1 = atom->map(tag[i] - 1);
                iH2 = i;
              }
            }
            if ((iH1 == -1) || (iH2 == -1)) error->one(FLERR, "TIP4P hydrogen is missing");
            if (iO == -1) error->one(FLERR, "TIP4P oxygen is missing");
            if ((type[iH1] != typeH) || (type[iH2] != typeH))
              error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
            if (type[iO] != typeO) error->one(FLERR, "TIP4P oxygen has incorrect atom type");
            iH1 = domain->closest_image(i, iH1);
            iH2 = domain->closest_image(i, iH2);

            find_M(x[iO], x[iH1], x[iH2], xM);

            // M site contributions

            if (!region || region->match(xM[0], xM[1], xM[2])) {

              if (xstyle == ATOM) {
                fx = qe2f * q[iO] * efield[iO][0];
              } else {
                fx = q[iO] * ex;
              }
              if (ystyle == ATOM) {
                fy = qe2f * q[iO] * efield[iO][1];
              } else {
                fy = q[iO] * ey;
              }
              if (zstyle == ATOM) {
                fz = qe2f * q[iO] * efield[iO][2];
              } else {
                fz = q[iO] * ez;
              }

              // distribute and apply forces, but only to local atoms

              if (i == iO) {
                f[iO][0] += fx * (1.0 - alpha);
                f[iO][1] += fy * (1.0 - alpha);
                f[iO][2] += fz * (1.0 - alpha);
                if (iH1 < nlocal) {
                  f[iH1][0] += 0.5 * alpha * fx;
                  f[iH1][1] += 0.5 * alpha * fy;
                  f[iH1][2] += 0.5 * alpha * fz;
                }
                if (iH2 < nlocal) {
                  f[iH2][0] += 0.5 * alpha * fx;
                  f[iH2][1] += 0.5 * alpha * fy;
                  f[iH2][2] += 0.5 * alpha * fz;
                }
              } else {
                if (iO >= nlocal) {
                  if (i == iH1) {
                    f[iH1][0] += 0.5 * alpha * fx;
                    f[iH1][1] += 0.5 * alpha * fy;
                    f[iH1][2] += 0.5 * alpha * fz;
                  }
                  if (i == iH2) {
                    f[iH2][0] += 0.5 * alpha * fx;
                    f[iH2][1] += 0.5 * alpha * fy;
                    f[iH2][2] += 0.5 * alpha * fz;
                  }
                }
              }

              if (i == iO) {
                if (estyle == ATOM) fsum[0] += efield[0][3];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contribution for point M with Oxygen

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iO, v);
                }
              }
            }

            // H1 site contributions

            if (!region || region->match(x[iH1][0], x[iH1][1], x[iH1][2])) {

              if (xstyle == ATOM) {
                fx = qe2f * q[iH1] * efield[iH1][0];
              } else {
                fx = q[iH1] * ex;
              }
              if (ystyle == ATOM) {
                fy = qe2f * q[iH1] * efield[iH1][1];
              } else {
                fy = q[iH1] * ey;
              }
              if (zstyle == ATOM) {
                fz = qe2f * q[iH1] * efield[iH1][2];
              } else {
                fz = q[iH1] * ez;
              }

              if (i == iH1) {
                f[iH1][0] += fx;
                f[iH1][1] += fy;
                f[iH1][2] += fz;

                // tally global force

                if (estyle == ATOM) fsum[0] += efield[0][3];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contribution

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iH1, v);
                }
              }
            }

            // H2 site contributions

            if (!region || region->match(x[iH2][0], x[iH2][1], x[iH2][2])) {

              if (xstyle == ATOM) {
                fx = qe2f * q[iH2] * efield[iH2][0];
              } else {
                fx = q[iH2] * ex;
              }
              if (ystyle == ATOM) {
                fy = qe2f * q[iH2] * efield[iH2][1];
              } else {
                fy = q[iH2] * ey;
              }
              if (zstyle == ATOM) {
                fz = qe2f * q[iH2] * efield[iH2][2];
              } else {
                fz = q[iH2] * ez;
              }

              if (i == iH2) {
                f[iH2][0] += fx;
                f[iH2][1] += fy;
                f[iH2][2] += fz;

                // tally global force

                if (estyle == ATOM) fsum[0] += efield[0][3];
                fsum[1] += fx;
                fsum[2] += fy;
                fsum[3] += fz;

                // tally virial contribution

                if (evflag) {
                  v[0] = fx * unwrap[0];
                  v[1] = fy * unwrap[1];
                  v[2] = fz * unwrap[2];
                  v[3] = fx * unwrap[1];
                  v[4] = fx * unwrap[2];
                  v[5] = fy * unwrap[2];
                  v_tally(iH2, v);
                }
              }
            }

          } else {

            // non-TIP4P charge interactions
            // force = qE

            if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

            if (xstyle == ATOM) {
              fx = qe2f * q[i] * efield[i][0];
            } else {
              fx = q[i] * ex;
            }
            f[i][0] += fx;
            fsum[1] += fx;
            if (ystyle == ATOM) {
              fy = qe2f * q[i] * efield[i][1];
            } else {
              fy = q[i] * ey;
            }
            f[i][1] += fy;
            fsum[2] += fy;
            if (zstyle == ATOM) {
              fz = qe2f * q[i] * efield[i][2];
            } else {
              fz = q[i] * ez;
            }
            f[i][2] += fz;
            fsum[3] += fz;

            if (estyle == ATOM) fsum[0] += efield[0][3];
          }
        }
      }
    }

    // dipole interactions
    // no force, torque = mu cross E

    if (muflag) {
      double **mu = atom->mu;
      double **t = atom->torque;
      double tx, ty, tz;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
          tx = ez * mu[i][1] - ey * mu[i][2];
          ty = ex * mu[i][2] - ez * mu[i][0];
          tz = ey * mu[i][0] - ex * mu[i][1];
          t[i][0] += tx;
          t[i][1] += ty;
          t[i][2] += tz;
        }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEfieldTIP4P::find_M(double *xO, double *xH1, double *xH2, double *xM)
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];

  xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
}
