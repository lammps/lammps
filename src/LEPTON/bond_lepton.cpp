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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "bond_lepton.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>
#include <exception>

#include "Lepton.h"
#include "lepton_utils.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondLepton::BondLepton(LAMMPS *_lmp) :
    Bond(_lmp), r0(nullptr), type2expression(nullptr), offset(nullptr)
{
  writedata = 1;
  reinitflag = 0;
  auto_offset = 1;
}

/* ---------------------------------------------------------------------- */

BondLepton::~BondLepton()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(type2expression);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void BondLepton::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  ev_init(eflag, vflag);
  if (evflag) {
    if (eflag) {
      if (force->newton_bond)
        eval<1, 1, 1>();
      else
        eval<1, 1, 0>();
    } else {
      if (force->newton_bond)
        eval<1, 0, 1>();
      else
        eval<1, 0, 0>();
    }
  } else {
    if (force->newton_bond)
      eval<0, 0, 1>();
    else
      eval<0, 0, 0>();
  }
}

/* ---------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG, int NEWTON_BOND> void BondLepton::eval()
{
  std::vector<Lepton::CompiledExpression> bondforce;
  std::vector<Lepton::CompiledExpression> bondpot;
  std::vector<bool> has_ref;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
      bondforce.emplace_back(parsed.differentiate("r").createCompiledExpression());
      has_ref.push_back(true);
      try {
        bondforce.back().getVariableReference("r");
      } catch (Lepton::Exception &) {
        has_ref.back() = false;
      }
      if (EFLAG) bondpot.emplace_back(parsed.createCompiledExpression());
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const int *const *const bondlist = neighbor->bondlist;
  const int nbondlist = neighbor->nbondlist;
  const int nlocal = atom->nlocal;

  for (int n = 0; n < nbondlist; n++) {
    const int i1 = bondlist[n][0];
    const int i2 = bondlist[n][1];
    const int type = bondlist[n][2];

    const double delx = x[i1][0] - x[i2][0];
    const double dely = x[i1][1] - x[i2][1];
    const double delz = x[i1][2] - x[i2][2];

    const double rsq = delx * delx + dely * dely + delz * delz;
    const double r = sqrt(rsq);
    const double dr = r - r0[type];
    const int idx = type2expression[type];

    // force and energy

    double fbond = 0.0;
    if (r > 0.0) {
      if (has_ref[idx]) bondforce[idx].getVariableReference("r") = dr;
      fbond = -bondforce[idx].evaluate() / r;
    }

    // apply force to each of 2 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    double ebond = 0.0;
    if (EFLAG) {
      try {
        bondpot[idx].getVariableReference("r") = dr;
      } catch (Lepton::Exception &) {
        ;    // ignore -> constant potential
      }
      ebond = bondpot[idx].evaluate() - offset[type];
    }
    if (EVFLAG) ev_tally(i1, i2, nlocal, NEWTON_BOND, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondLepton::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(r0, np1, "bond:r0");
  memory->create(type2expression, np1, "bond:type2expression");
  memory->create(offset, np1, "bond:offset");
  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void BondLepton::settings(int narg, char **arg)
{
  auto_offset = 1;
  if (narg > 0) {
    if (strcmp(arg[0],"auto_offset") == 0) {
      auto_offset = 1;
    } else if (strcmp(arg[0],"no_offset") == 0) {
      auto_offset = 0;
    } else {
      error->all(FLERR, "Unknown bond style lepton setting {}", arg[0]);
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondLepton::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect number of args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double r0_one = utils::numeric(FLERR, arg[1], false, lmp);

  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed and evaluated without error
  std::string exp_one = LeptonUtils::condense(arg[2]);
  double offset_one = 0.0;
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
    auto bondpot = parsed.createCompiledExpression();
    auto bondforce = parsed.differentiate("r").createCompiledExpression();
    try {
      bondpot.getVariableReference("r") = 0.0;
    } catch (Lepton::Exception &e) {
      if (comm->me == 0)
        error->warning(FLERR, "Lepton potential expression {} does not depend on 'r'", exp_one);
    }
    try {
      bondforce.getVariableReference("r") = 0.0;
    } catch (Lepton::Exception &e) {
      if (comm->me == 0)
        error->warning(FLERR, "Force from Lepton expression {} does not depend on 'r'", exp_one);
    }
    if (auto_offset) offset_one = bondpot.evaluate();
    bondforce.evaluate();
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  std::size_t idx = 0;
  for (const auto &exp : expressions) {
    if (exp == exp_one) break;
    ++idx;
  }

  // if not found, add to list
  if ((expressions.size() == 0) || (idx == expressions.size())) expressions.push_back(exp_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    r0[i] = r0_one;
    type2expression[i] = idx;
    offset[i] = offset_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondLepton::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondLepton::write_restart(FILE *fp)
{
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&type2expression[1], sizeof(int), atom->nbondtypes, fp);
  fwrite(&offset[1], sizeof(double), atom->nbondtypes, fp);

  int num = expressions.size();
  int maxlen = 0;
  for (const auto &exp : expressions) maxlen = MAX(maxlen, (int) exp.size());
  ++maxlen;

  fwrite(&num, sizeof(int), 1, fp);
  fwrite(&maxlen, sizeof(int), 1, fp);
  for (const auto &exp : expressions) {
    int n = exp.size() + 1;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(exp.c_str(), sizeof(char), n, fp);
  }
  fwrite(&auto_offset, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondLepton::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &type2expression[1], sizeof(int), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &offset[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&type2expression[1], atom->nbondtypes, MPI_INT, 0, world);
  MPI_Bcast(&offset[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;

  int num, maxlen, len;
  if (comm->me == 0) {
    utils::sfread(FLERR, &num, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &maxlen, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&num, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxlen, 1, MPI_INT, 0, world);

  char *buf = new char[maxlen];

  for (int i = 0; i < num; ++i) {
    if (comm->me == 0) {
      utils::sfread(FLERR, &len, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, buf, sizeof(char), len, fp, nullptr, error);
    }
    MPI_Bcast(buf, maxlen, MPI_CHAR, 0, world);
    expressions.emplace_back(buf);
  }

  if (comm->me == 0) utils::sfread(FLERR, &auto_offset, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&auto_offset, 1, MPI_INT, 0, world);

  delete[] buf;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondLepton::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g %s\n", i, r0[i], expressions[type2expression[i]].c_str());
}

/* ---------------------------------------------------------------------- */

double BondLepton::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  const double r = sqrt(rsq);
  const double dr = r - r0[type];

  const auto &expr = expressions[type2expression[type]];
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
  auto bondpot = parsed.createCompiledExpression();
  auto bondforce = parsed.differentiate("r").createCompiledExpression();
  try {
    bondpot.getVariableReference("r") = dr;
    bondforce.getVariableReference("r") = dr;
  } catch (Lepton::Exception &) {
    ;    // ignore -> constant potential or force
  }

  // force and energy

  fforce = 0.0;
  double ebond = 0.0;
  if (r > 0.0) {
    fforce = -bondforce.evaluate() / r;
    ebond = bondpot.evaluate() - offset[type];
  }
  return ebond;
}

/* ---------------------------------------------------------------------- */

void *BondLepton::extract(const char *str, int &dim)
{
  dim = 1;
  if (str) {
    std::string keyword(str);
    if (keyword == "r0") return (void *) r0;
  }
  return nullptr;
}
