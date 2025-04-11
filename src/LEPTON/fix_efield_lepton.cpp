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
   Contributing author: Gabriel Alkuino (Syracuse University) - gsalkuin@syr.edu
   Modified from fix_efield
------------------------------------------------------------------------- */

#include "fix_efield_lepton.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"

#include <array>

#include "Lepton.h"
#include "lepton_utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

FixEfieldLepton::FixEfieldLepton(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), region(nullptr)
{
  if (domain->xperiodic || domain->yperiodic || domain->zperiodic) {
    if (comm->me == 0)
      error->warning(FLERR, "Fix {} uses unwrapped coordinates", style);
  }
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  qe2f = force->qe2f;
  mue2e = qe2f;

  // optional args
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " region", error);
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix {} does not exist", arg[iarg + 1], style);
      delete[] idregion;
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "step") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + "step", error);
      h = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown keyword for fix {} command: {}", style, arg[iarg]);
    }
  }

  // check validity of Lepton expression
  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed without error
  expr = LeptonUtils::condense(arg[3]);
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
    auto phi = parsed.createCompiledExpression();
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixEfieldLepton::~FixEfieldLepton()
{
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixEfieldLepton::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEfieldLepton::init()
{
  if (!atom->q_flag && !atom->mu_flag)
    error->all(FLERR, "Fix {} requires atom attribute q or mu", style);
  if (atom->mu_flag && !atom->torque_flag)
    error->all(FLERR, "Dipoles must be finite-sized to rotate", style);

  // set index and check validity of region
  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix {} does not exist", idregion, style);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto respa = dynamic_cast<Respa *>(update->integrate);
    if (respa) ilevel_respa = respa->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }

  // unit conversion restrictions (see issue #1377)
  char *unit_style = update->unit_style;
  if (strcmp(unit_style, "electron") == 0 || strcmp(unit_style, "micro") == 0 ||
      strcmp(unit_style, "nano") == 0) {
    error->all(FLERR, "Fix {} does not support {} units", style, unit_style);
  }
}

/* ---------------------------------------------------------------------- */

void FixEfieldLepton::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto respa = dynamic_cast<Respa *>(update->integrate);
    if (respa) {
      respa->copy_flevel_f(ilevel_respa);
      post_force_respa(vflag, ilevel_respa, 0);
      respa->copy_f_flevel(ilevel_respa);
    }
  } else {
    post_force(vflag);
  }
}

/* ---------------------------------------------------------------------- */

void FixEfieldLepton::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
  Apply F = qE,
        F = (mu . D) E,
        T = mu x E
------------------------------------------------------------------------- */

void FixEfieldLepton::post_force(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp)).optimize();
  Lepton::CompiledExpression phi;
  auto dphi_x = parsed.differentiate("x").createCompiledExpression();
  auto dphi_y = parsed.differentiate("y").createCompiledExpression();
  auto dphi_z = parsed.differentiate("z").createCompiledExpression();
  std::array<Lepton::CompiledExpression *, 3> dphis = {&dphi_x, &dphi_y, &dphi_z};

  // array of vectors of ptrs to Lepton variable references
  std::array<std::vector<double *>, 3> var_ref_ptrs{};

  // fill ptr-vectors with Lepton refs as needed
  const char *DIM_NAMES[] = {"x", "y", "z"};
  if (atom->q_flag) {
    phi = parsed.createCompiledExpression();
    for (size_t d = 0; d < 3; d++) {
      try {
        double *ptr = &(phi.getVariableReference(DIM_NAMES[d]));
        var_ref_ptrs[d].push_back(ptr);
      } catch (Lepton::Exception &) {
        // do nothing
      }
    }
  }

  bool e_uniform = true;
  for (size_t j = 0; j < 3; j++)
    for (size_t d = 0; d < 3; d++) {
      try {
        double *ptr = &((*dphis[j]).getVariableReference(DIM_NAMES[d]));
        var_ref_ptrs[d].push_back(ptr);
        e_uniform = false;
      } catch (Lepton::Exception &) {
        // do nothing
      }
    }
  if (!e_uniform && atom->mu_flag && h < 0) {
    error->all(FLERR, "Fix {} requires keyword `step' for dipoles in a non-uniform electric field",
               style);
  }

  // virial setup
  v_init(vflag);

  // update region if necessary
  if (region) region->prematch();

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double ex, ey, ez;
  double fx, fy, fz;
  double v[6], unwrap[3], dstep[3];
  double exf, eyf, ezf, exb, eyb, ezb;
  double mu_norm, h_mu;

  double *q = atom->q;
  double **mu = atom->mu;
  double **t = atom->torque;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      fx = fy = fz = 0.0;
      domain->unmap(x[i], image[i], unwrap);

      // put unwrapped coords into Lepton variable refs
      for (size_t d = 0; d < 3; d++) {
        for (auto &var_ref_ptr : var_ref_ptrs[d]) { *var_ref_ptr = unwrap[d]; }
      }

      // evaluate e-field, used by q and mu
      ex = -dphi_x.evaluate();
      ey = -dphi_y.evaluate();
      ez = -dphi_z.evaluate();

      // charges
      // force = q E
      if (atom->q_flag) {
        fx = qe2f * q[i] * ex;
        fy = qe2f * q[i] * ey;
        fz = qe2f * q[i] * ez;
        // potential energy = q phi
        fsum[0] += qe2f * q[i] * phi.evaluate();
      }

      if (atom->mu_flag) {
        // dipoles
        mu_norm = sqrt(mu[i][0] * mu[i][0] + mu[i][1] * mu[i][1] + mu[i][2] * mu[i][2]);
        if (mu_norm > EPSILON) {
          // torque = mu cross E
          t[i][0] += mue2e * (ez * mu[i][1] - ey * mu[i][2]);
          t[i][1] += mue2e * (ex * mu[i][2] - ez * mu[i][0]);
          t[i][2] += mue2e * (ey * mu[i][0] - ex * mu[i][1]);
          // potential energy = - mu dot E
          fsum[0] -= mue2e * (mu[i][0] * ex + mu[i][1] * ey + mu[i][2] * ez);

          // force = (mu dot D) E for non-uniform E
          // using central difference method
          if (!e_uniform) {
            h_mu = h / mu_norm;
            dstep[0] = h_mu * mu[i][0];
            dstep[1] = h_mu * mu[i][1];
            dstep[2] = h_mu * mu[i][2];

            // one step forwards, two steps back ;)
            for (size_t d = 0; d < 3; d++) {
              for (auto &var_ref_ptr : var_ref_ptrs[d]) { *var_ref_ptr += dstep[d]; }
            }

            exf = -dphi_x.evaluate();
            eyf = -dphi_y.evaluate();
            ezf = -dphi_z.evaluate();

            for (size_t d = 0; d < 3; d++) {
              for (auto &var_ref_ptr : var_ref_ptrs[d]) { *var_ref_ptr -= 2 * dstep[d]; }
            }

            exb = -dphi_x.evaluate();
            eyb = -dphi_y.evaluate();
            ezb = -dphi_z.evaluate();

            fx += qe2f * (exf - exb) / 2.0 / h_mu;
            fy += qe2f * (eyf - eyb) / 2.0 / h_mu;
            fz += qe2f * (ezf - ezb) / 2.0 / h_mu;
          }
        }
      }

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

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

/* ---------------------------------------------------------------------- */

void FixEfieldLepton::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEfieldLepton::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixEfieldLepton::compute_scalar()
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum, fsum_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsum_all[0];
}

/* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

double FixEfieldLepton::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum, fsum_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsum_all[n + 1];
}
