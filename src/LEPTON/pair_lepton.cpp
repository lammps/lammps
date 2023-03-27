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

#include "pair_lepton.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include "Lepton.h"
#include "lepton_utils.h"
#include <cmath>
#include <map>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLepton::PairLepton(LAMMPS *lmp) :
    Pair(lmp), cut(nullptr), type2expression(nullptr), offset(nullptr)
{
  respa_enable = 0;
  single_enable = 1;
  writedata = 1;
  restartinfo = 1;
  reinitflag = 0;
  cut_global = 0.0;
  centroidstressflag = CENTROID_SAME;

  functions["zbl"] = new Lepton::ZBLFunction(force->qqr2e, force->angstrom, force->qelectron);
}

/* ---------------------------------------------------------------------- */

PairLepton::~PairLepton()
{
  for (auto &f : functions) delete f.second;
  if (allocated) {
    memory->destroy(cut);
    memory->destroy(cutsq);
    memory->destroy(setflag);
    memory->destroy(type2expression);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLepton::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  if (evflag) {
    if (eflag) {
      if (force->newton_pair)
        eval<1, 1, 1>();
      else
        eval<1, 1, 0>();
    } else {
      if (force->newton_pair)
        eval<1, 0, 1>();
      else
        eval<1, 0, 0>();
    }
  } else {
    if (force->newton_pair)
      eval<0, 0, 1>();
    else
      eval<0, 0, 0>();
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void PairLepton::eval()
{
  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *const special_lj = force->special_lj;

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  double fxtmp, fytmp, fztmp;

  std::vector<Lepton::CompiledExpression> pairforce;
  std::vector<Lepton::CompiledExpression> pairpot;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp), functions);
      pairforce.emplace_back(parsed.differentiate("r").createCompiledExpression());
      pairforce.back().getVariableReference("r");
      if (EFLAG) pairpot.emplace_back(parsed.createCompiledExpression());
    }
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  // loop over neighbors of my atoms

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const int jtype = type[j];

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        const int idx = type2expression[itype][jtype];
        double &r_for = pairforce[idx].getVariableReference("r");
        r_for = r;
        const double fpair = -pairforce[idx].evaluate() / r * factor_lj;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || (j < nlocal)) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        double evdwl = 0.0;
        if (EFLAG) {
          pairpot[idx].getVariableReference("r") = r;
          evdwl = pairpot[idx].evaluate() - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (EVFLAG) ev_tally(i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLepton::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cut, np1, np1, "pair:cut");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(type2expression, np1, np1, "pair:type2expression");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLepton::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Incorrect number of arguments for pair_style lepton command");
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
------------------------------------------------------------------------- */

void PairLepton::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all(FLERR, "Incorrect number of args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = cut_global;
  if (narg == 4) {
    if (pppmflag || ewaldflag || msmflag || dispersionflag || tip4pflag) {
      error->all(FLERR, "Only a global cutoff is allowed with Kspace compatibility enabled");
    } else {
      cut_one = utils::numeric(FLERR, arg[3], false, lmp);
    }
  }

  // remove whitespace and quotes from expression string and then
  // check if the expression can be parsed and evaluated without error
  auto exp_one = LeptonUtils::condense(arg[2]);
  try {
    auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(exp_one, lmp), functions);
    auto pairforce = parsed.differentiate("r").createCompiledExpression();
    auto pairpot = parsed.createCompiledExpression();
    pairpot.getVariableReference("r") = 1.0;
    pairforce.getVariableReference("r") = 1.0;
    pairpot.evaluate();
    pairforce.evaluate();
  } catch (std::exception &e) {
    error->all(FLERR, e.what());
  }

  std::size_t idx = 0;
  for (const auto &exp : expressions) {
    if (exp == exp_one) break;
    ++idx;
  }

  // not found, add to list
  if ((expressions.size() == 0) || (idx == expressions.size())) expressions.push_back(exp_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      type2expression[i][j] = idx;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

double PairLepton::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  offset[i][j] = 0.0;
  if (offset_flag) {
    try {
      auto expr = LeptonUtils::substitute(expressions[type2expression[i][j]], lmp);
      auto pairpot = Lepton::Parser::parse(expr, functions).createCompiledExpression();
      pairpot.getVariableReference("r") = cut[i][j];
      offset[i][j] = pairpot.evaluate();
    } catch (std::exception &) {
    }
  }

  cut[j][i] = cut[i][j];
  type2expression[j][i] = type2expression[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLepton::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j], sizeof(double), 1, fp);
        fwrite(&type2expression[i][j], sizeof(int), 1, fp);
      }
    }

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
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLepton::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();
  expressions.clear();

  const int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &type2expression[i][j], sizeof(int), 1, fp, nullptr, error);
        }
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&type2expression[i][j], 1, MPI_INT, 0, world);
      }
    }

  int num, maxlen, len;
  if (me == 0) {
    utils::sfread(FLERR, &num, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &maxlen, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&num, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxlen, 1, MPI_INT, 0, world);

  char *buf = new char[maxlen];

  for (int i = 0; i < num; ++i) {
    if (me == 0) {
      utils::sfread(FLERR, &len, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, buf, sizeof(char), len, fp, nullptr, error);
    }
    MPI_Bcast(buf, maxlen, MPI_CHAR, 0, world);
    expressions.emplace_back(buf);
  }

  delete[] buf;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLepton::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLepton::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLepton::write_data(FILE *fp)
{
  if (pppmflag || ewaldflag || msmflag || dispersionflag || tip4pflag) {
    for (int i = 1; i <= atom->ntypes; i++)
      fprintf(fp, "%d %s\n", i, expressions[type2expression[i][i]].c_str());
  } else {
    for (int i = 1; i <= atom->ntypes; i++)
      fprintf(fp, "%d %s %g\n", i, expressions[type2expression[i][i]].c_str(), cut[i][i]);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLepton::write_data_all(FILE *fp)
{
  if (pppmflag || ewaldflag || msmflag || dispersionflag || tip4pflag) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        fprintf(fp, "%d %d %s\n", i, j, expressions[type2expression[i][j]].c_str());
  } else {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        fprintf(fp, "%d %d %s %g\n", i, j, expressions[type2expression[i][j]].c_str(), cut[i][j]);
  }
}

/* ---------------------------------------------------------------------- */

double PairLepton::single(int /* i */, int /* j */, int itype, int jtype, double rsq,
                          double /* factor_coul */, double factor_lj, double &fforce)
{
  auto expr = expressions[type2expression[itype][jtype]];
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp), functions);
  auto pairpot = parsed.createCompiledExpression();
  auto pairforce = parsed.differentiate("r").createCompiledExpression();

  const double r = sqrt(rsq);
  pairpot.getVariableReference("r") = r;
  pairforce.getVariableReference("r") = r;

  fforce = -pairforce.evaluate() / r * factor_lj;
  return (pairpot.evaluate() - offset[itype][jtype]) * factor_lj;
}
