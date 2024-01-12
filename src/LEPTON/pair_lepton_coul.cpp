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

#include "pair_lepton_coul.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include "Lepton.h"
#include "lepton_utils.h"
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void PairLeptonCoul::compute(int eflag, int vflag)
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

template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void PairLeptonCoul::eval()
{
  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const double *const q = atom->q;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *const special_coul = force->special_coul;

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  double fxtmp, fytmp, fztmp;

  const double q2e = sqrt(force->qqrd2e);

  std::vector<Lepton::CompiledExpression> pairforce;
  std::vector<Lepton::CompiledExpression> pairpot;
  std::vector<std::pair<bool, bool>> have_q;
  try {
    for (const auto &expr : expressions) {
      auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp), functions);
      pairforce.emplace_back(parsed.differentiate("r").createCompiledExpression());
      if (EFLAG) pairpot.emplace_back(parsed.createCompiledExpression());
      pairforce.back().getVariableReference("r");
      have_q.emplace_back(true, true);

      // check if there are references to charges
      try {
        pairforce.back().getVariableReference("qi");
      } catch (std::exception &) {
        have_q.back().first = false;
      }
      try {
        pairforce.back().getVariableReference("qj");
      } catch (std::exception &) {
        have_q.back().second = false;
      }
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
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      const int jtype = type[j];

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        const int idx = type2expression[itype][jtype];
        pairforce[idx].getVariableReference("r") = r;
        if (have_q[idx].first) pairforce[idx].getVariableReference("qi") = q2e * q[i];
        if (have_q[idx].second) pairforce[idx].getVariableReference("qj") = q2e * q[j];
        const double fpair = -pairforce[idx].evaluate() / r * factor_coul;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || (j < nlocal)) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        double ecoul = 0.0;
        if (EFLAG) {
          pairpot[idx].getVariableReference("r") = r;
          if (have_q[idx].first) pairpot[idx].getVariableReference("qi") = q2e * q[i];
          if (have_q[idx].second) pairpot[idx].getVariableReference("qj") = q2e * q[j];
          ecoul = pairpot[idx].evaluate();
          ecoul *= factor_coul;
        }

        if (EVFLAG) ev_tally(i, j, nlocal, NEWTON_PAIR, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLeptonCoul::settings(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "pair_style lepton/coul", error);
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // optional keywords
  // assert the pair style is compatible with a specific long-range solver

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ewald") == 0)
      ewaldflag = 1;
    else if (strcmp(arg[iarg], "pppm") == 0)
      pppmflag = 1;
    else if (strcmp(arg[iarg], "msm") == 0)
      msmflag = 1;
    else if (strcmp(arg[iarg], "dispersion") == 0)
      dispersionflag = 1;
    else if (strcmp(arg[iarg], "tip4p") == 0)
      tip4pflag = 1;
    else
      error->all(FLERR, "Unknown pair_style lepton/coul keyword: {}", arg[iarg]);
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

void PairLeptonCoul::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style lepton/coul requires atom attribute q");
  if (offset_flag) error->all(FLERR, "Pair style lepton/coul does not suport pair_modify shift");
  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLeptonCoul::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&ewaldflag, sizeof(int), 1, fp);
  fwrite(&pppmflag, sizeof(int), 1, fp);
  fwrite(&msmflag, sizeof(int), 1, fp);
  fwrite(&dispersionflag, sizeof(int), 1, fp);
  fwrite(&tip4pflag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLeptonCoul::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &ewaldflag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &pppmflag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &msmflag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &dispersionflag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tip4pflag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ewaldflag, 1, MPI_INT, 0, world);
  MPI_Bcast(&pppmflag, 1, MPI_INT, 0, world);
  MPI_Bcast(&msmflag, 1, MPI_INT, 0, world);
  MPI_Bcast(&dispersionflag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tip4pflag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairLeptonCoul::single(int i, int j, int itype, int jtype, double rsq, double factor_coul,
                              double /* factor_lj */, double &fforce)
{
  const auto &expr = expressions[type2expression[itype][jtype]];
  auto parsed = Lepton::Parser::parse(LeptonUtils::substitute(expr, lmp), functions);
  auto pairpot = parsed.createCompiledExpression();
  auto pairforce = parsed.differentiate("r").createCompiledExpression();

  const double r = sqrt(rsq);
  const double q2e = sqrt(force->qqrd2e);
  pairpot.getVariableReference("r") = r;
  pairforce.getVariableReference("r") = r;
  try {
    pairpot.getVariableReference("qi") = q2e * atom->q[i];
    pairforce.getVariableReference("qi") = q2e * atom->q[i];
  } catch (std::exception &) {
    /* ignore */
  }
  try {
    pairpot.getVariableReference("qj") = q2e * atom->q[j];
    pairforce.getVariableReference("qj") = q2e * atom->q[j];
  } catch (std::exception &) {
    /* ignore */
  }

  fforce = -pairforce.evaluate() / r * factor_coul;
  return pairpot.evaluate() * factor_coul;
}

/* ---------------------------------------------------------------------- */

void *PairLeptonCoul::extract(const char *str, int &dim)
{
  if (pppmflag || ewaldflag || msmflag || dispersionflag || tip4pflag) {
    if (strcmp(str, "cut_coul") == 0) {
      dim = 0;
      return (void *) &cut_global;
    }
  } else {
    if (strcmp(str, "cut_coul") == 0) {
      dim = 2;
      return (void *) &cut;
    }
  }
  return nullptr;
}
