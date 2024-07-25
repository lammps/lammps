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

#include "fix_epot_lepton.h"
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

#include "Lepton.h"
#include "lepton_utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEpotLepton::FixEpotLepton(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), 
    idregion(nullptr), region(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
  
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // optional region keyword
  if (narg > 4) {
    if (strcmp(arg[4], "region") == 0) {
      if (narg != 6)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " region", error);
      region = domain->get_region_by_id(arg[5]);
      if (!region) error->all(FLERR, "Region {} for fix epot/lepton does not exist.", arg[5]);
      idregion = utils::strdup(arg[5]);
    }
  }

  // check validity of lepton expression
  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed without error
  expr = LeptonUtils::condense(arg[3]);
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
    auto phi = parsed.createCompiledExpression();
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }
}

/* ---------------------------------------------------------------------- */

FixEpotLepton::~FixEpotLepton()
{
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixEpotLepton::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEpotLepton::init()
{
  if (!atom->q_flag && !atom->mu_flag) error->all(FLERR, "Fix {} requires atom attribute q or mu", style);
  if (atom->mu_flag && !atom->torque_flag) error->all(FLERR, "Dipoles must be finite-sized to rotate", style);

  // set index and check validity of region
  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix {} does not exist", idregion, style);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixEpotLepton::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    respa->copy_f_flevel(ilevel_respa);
  } else {
    post_force(vflag);
  }
}

/* ---------------------------------------------------------------------- */

void FixEpotLepton::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
  Apply F = qE,
        F = (mu . D) E, 
        T = mu x E 
------------------------------------------------------------------------- */

void FixEpotLepton::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  const double qe2f = force->qe2f;
  const double qqr2e = force->qqr2e;

  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp)).optimize();
  auto phi = parsed.createCompiledExpression();
  auto dphi_x = parsed.differentiate("x").createCompiledExpression();
  auto dphi_y = parsed.differentiate("y").createCompiledExpression();
  auto dphi_z = parsed.differentiate("z").createCompiledExpression();
  auto dphi_xx = parsed.differentiate("x").differentiate("x").createCompiledExpression();
  auto dphi_xy = parsed.differentiate("x").differentiate("y").createCompiledExpression();
  auto dphi_xz = parsed.differentiate("x").differentiate("z").createCompiledExpression();
  auto dphi_yy = parsed.differentiate("y").differentiate("y").createCompiledExpression();
  auto dphi_yz = parsed.differentiate("y").differentiate("z").createCompiledExpression();
  auto dphi_zz = parsed.differentiate("z").differentiate("z").createCompiledExpression();

  std::vector<Lepton::CompiledExpression*> cexprs = {&phi, &dphi_x, &dphi_y, &dphi_z, &dphi_xx, &dphi_xy, &dphi_xz, &dphi_yy, &dphi_yz, &dphi_zz}; 
  
  // virial setup
  v_init(vflag);

  // update region if necessary
  if (region) region->prematch();

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double *q = nullptr;
  if (atom->q_flag) q = atom->q;
  
  double **mu = nullptr, **t = nullptr;
  if (atom->mu_flag) {mu = atom->mu; t = atom->torque;}

  double **x = atom->x;
  double fx, fy, fz;
  double v[6], unwrap[3];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      domain->unmap(x[i], image[i], unwrap);
      
      // Substitute x, y, z if they exist           
      for (auto &cexpr : cexprs) {
        try { 
          (*cexpr).getVariableReference("x") = unwrap[0];
        }
        catch (Lepton::Exception &) {
          // ignore
        }
        try { 
          (*cexpr).getVariableReference("y") = unwrap[1];
        }
        catch (Lepton::Exception &) {
          // ignore
        }
        try { 
          (*cexpr).getVariableReference("z") = unwrap[2];
        }
        catch (Lepton::Exception &) {
          // ignore
        }
      }

      // charges
      // force = q E, potential energy = q phi
      if (atom->q_flag) {
        fx = - qe2f * q[i] * dphi_x.evaluate();
        fy = - qe2f * q[i] * dphi_y.evaluate();
        fz = - qe2f * q[i] * dphi_z.evaluate();
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        fsum[0] += qqr2e * q[i] * phi.evaluate();
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

      // dipoles
      // force = (mu dot D) E, torque = mu cross E, potential energy = - mu dot E
      if (atom->mu_flag) {
        double tx, ty, tz;
        auto ex = - dphi_x.evaluate();
        auto ey = - dphi_y.evaluate();
        auto ez = - dphi_z.evaluate();
        tx = ez * mu[i][1] - ey * mu[i][2];
        ty = ex * mu[i][2] - ez * mu[i][0];
        tz = ey * mu[i][0] - ex * mu[i][1];
        t[i][0] += qqr2e * tx;
        t[i][1] += qqr2e * ty;
        t[i][2] += qqr2e * tz;
        fsum[0] -= qqr2e * (mu[i][0] * ex + mu[i][1] * ey + mu[i][2] * ez);

        fx = - qe2f * ( mu[i][0] * dphi_xx.evaluate() + mu[i][1] * dphi_xy.evaluate() + mu[i][2] * dphi_xz.evaluate() );
        fy = - qe2f * ( mu[i][0] * dphi_xy.evaluate() + mu[i][1] * dphi_yy.evaluate() + mu[i][2] * dphi_yz.evaluate() );
        fz = - qe2f * ( mu[i][0] * dphi_xz.evaluate() + mu[i][1] * dphi_yz.evaluate() + mu[i][2] * dphi_zz.evaluate() );
        
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        fsum[1] += fx;
        fsum[2] += fy;
        fsum[3] += fz;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEpotLepton::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEpotLepton::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixEpotLepton::compute_scalar()
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

double FixEpotLepton::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum, fsum_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fsum_all[n + 1];
}
