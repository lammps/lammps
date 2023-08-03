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

#include "angle_lepton.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>

#include "Lepton.h"
#include "lepton_utils.h"

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleLepton::AngleLepton(LAMMPS *_lmp) :
    Angle(_lmp), theta0(nullptr), type2expression(nullptr), offset(nullptr)
{
  writedata = 1;
  reinitflag = 0;
}

/* ---------------------------------------------------------------------- */

AngleLepton::~AngleLepton()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(theta0);
    memory->destroy(type2expression);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void AngleLepton::compute(int eflag, int vflag)
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

template <int EVFLAG, int EFLAG, int NEWTON_BOND> void AngleLepton::eval()
{
  std::vector<Lepton::CompiledExpression> angleforce;
  std::vector<Lepton::CompiledExpression> anglepot;
  for (const auto &expr : expressions) {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
    angleforce.emplace_back(parsed.differentiate("theta").createCompiledExpression());
    if (EFLAG) anglepot.emplace_back(parsed.createCompiledExpression());
  }

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const int *const *const anglelist = neighbor->anglelist;
  const int nanglelist = neighbor->nanglelist;
  const int nlocal = atom->nlocal;

  for (int n = 0; n < nanglelist; n++) {
    const int i1 = anglelist[n][0];
    const int i2 = anglelist[n][1];
    const int i3 = anglelist[n][2];
    const int type = anglelist[n][3];

    // 1st bond

    const double delx1 = x[i1][0] - x[i2][0];
    const double dely1 = x[i1][1] - x[i2][1];
    const double delz1 = x[i1][2] - x[i2][2];

    const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    const double r1 = sqrt(rsq1);

    // 2nd bond

    const double delx2 = x[i3][0] - x[i2][0];
    const double dely2 = x[i3][1] - x[i2][1];
    const double delz2 = x[i3][2] - x[i2][2];

    const double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    const double r2 = sqrt(rsq2);

    // angle (cos and sin)

    double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    // force and energy

    const double dtheta = acos(c) - theta0[type];
    const int idx = type2expression[type];
    angleforce[idx].getVariableReference("theta") = dtheta;

    const double a = -angleforce[idx].evaluate() * s;
    const double a11 = a * c / rsq1;
    const double a12 = -a / (r1 * r2);
    const double a22 = a * c / rsq2;

    double f1[3], f3[3];
    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    // apply force to each of 3 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    double eangle = 0.0;
    if (EFLAG) {
      anglepot[idx].getVariableReference("theta") = dtheta;
      eangle = anglepot[idx].evaluate() - offset[type];
    }
    if (EVFLAG)
      ev_tally(i1, i2, i3, nlocal, NEWTON_BOND, eangle, f1, f3, delx1, dely1, delz1, delx2, dely2,
               delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleLepton::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(theta0, np1, "angle:theta0");
  memory->create(type2expression, np1, "angle:type2expression");
  memory->create(offset, np1, "angle:offset");
  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleLepton::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect number of args for angle coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double theta0_one = utils::numeric(FLERR, arg[1], false, lmp);

  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed and evaluated without error
  std::string exp_one = LeptonUtils::condense(arg[2]);
  double offset_one = 0.0;
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp));
    auto anglepot = parsed.createCompiledExpression();
    auto angleforce = parsed.differentiate("theta").createCompiledExpression();
    anglepot.getVariableReference("theta") = 0.0;
    angleforce.getVariableReference("theta") = 0.0;
    offset_one = anglepot.evaluate();
    angleforce.evaluate();
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

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    theta0[i] = DEG2RAD * theta0_one;
    type2expression[i] = idx;
    offset[i] = offset_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleLepton::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleLepton::write_restart(FILE *fp)
{
  fwrite(&theta0[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&type2expression[1], sizeof(int), atom->nangletypes, fp);
  fwrite(&offset[1], sizeof(double), atom->nangletypes, fp);

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
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleLepton::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &theta0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &type2expression[1], sizeof(int), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &offset[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&theta0[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&type2expression[1], atom->nangletypes, MPI_INT, 0, world);
  MPI_Bcast(&offset[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;

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

  delete[] buf;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleLepton::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %g %s\n", i, RAD2DEG * theta0[i], expressions[type2expression[i]].c_str());
}

/* ---------------------------------------------------------------------- */

double AngleLepton::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1, dely1, delz1);
  double r1 = sqrt(delx1 * delx1 + dely1 * dely1 + delz1 * delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2, dely2, delz2);
  double r2 = sqrt(delx2 * delx2 + dely2 * dely2 + delz2 * delz2);

  double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
  c /= r1 * r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  auto expr = expressions[type2expression[type]];
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp));
  auto anglepot = parsed.createCompiledExpression();
  anglepot.getVariableReference("theta") = dtheta;
  return anglepot.evaluate() - offset[type];
}

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *AngleLepton::extract(const char *str, int &dim)
{
  dim = 1;
  if (str) {
    std::string keyword(str);
    if (keyword == "theta0") return (void *) theta0;
  }
  return nullptr;
}
