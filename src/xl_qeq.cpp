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

#include "xl_qeq.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "lammps.h"
#include "update.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;

static const char cite_xl_qeq[] =
    "extended-Lagrangian charge equilibration: doi:10.1016/j.cpc.2015.02.023\n\n"
    "@Article{Nomura2015,\n"
    " author =  {K. Nomura and P. E. Small and R. K. Kalia and A. Nakano and P. Vashishta},\n"
    " title =   {An Extended-Lagrangian Scheme for Charge Equilibration in Reactive Molecular\n"
    "    Dynamics Simulations},\n"
    " journal = {Computer Physics Communications},\n"
    " year =    2015,\n"
    " volume =  192,\n"
    " pages =   {91--96}\n"
    "}\n\n";

static const char cite_xl_qeq_dissipation[] =
    "extended-Lagrangian charge equilibration with dissipation: doi:10.1063/1.3148075\n\n"
    "@Article{Niklasson2009,\n"
    " author =  {A. M. N. Niklasson and P. Steneteg and A. Odell and N. Bock and\n"
    "    M. Challacombe and C. J. Tymczak and E. Holmstr{\\\"o}m and G. Zheng and V. Weber},\n"
    " title =   {Extended {Lagrangian} {Born--Oppenheimer} Molecular Dynamics with Dissipation},\n"
    " journal = {Journal of Chemical Physics},\n"
    " year =    2009,\n"
    " volume =  130,\n"
    " pages =   {214109}\n"
    "}\n\n";

// dissipative modified-Verlet coefficients kappa, alpha, c_0 ... c_K for the
// auxiliary variable propagation, from Niklasson et al., JCP 130, 214109 (2009), Table I.
// row index = K - KDIS_MIN

static constexpr int KDIS_MIN = 3;
static constexpr int KDIS_MAX = 9;
static const struct {
  double kappa, alpha;
  double c[XLQEq::MAXHIST];
} dis_table[] = {
    {1.69, 0.150, {-2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {1.75, 0.057, {-3.0, 6.0, -2.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {1.82, 0.018, {-6.0, 14.0, -8.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0, 0.0}},
    {1.84, 0.0055, {-14.0, 36.0, -27.0, -2.0, 12.0, -6.0, 1.0, 0.0, 0.0, 0.0}},
    {1.86, 0.0016, {-36.0, 99.0, -88.0, 11.0, 32.0, -25.0, 8.0, -1.0, 0.0, 0.0}},
    {1.88, 0.00044, {-99.0, 286.0, -286.0, 78.0, 78.0, -90.0, 42.0, -10.0, 1.0, 0.0}},
    {1.89, 0.00012, {-286.0, 858.0, -936.0, 364.0, 168.0, -300.0, 184.0, -63.0, 12.0, -1.0}},
};

/* ---------------------------------------------------------------------- */

int XLQEq::parse_keyword(LAMMPS *lmp, const char *style, int narg, char **arg, int iarg)
{
  if (strcmp(arg[iarg], "xlcg") == 0) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix ") + style + " xlcg", lmp->error);
    ncg = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if (ncg < 1) lmp->error->all(FLERR, iarg + 1, "Fix {} xlcg value {} must be >= 1", style, ncg);
    return 2;
  }
  if (strcmp(arg[iarg], "xlkappa") == 0) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix ") + style + " xlkappa", lmp->error);
    kappa = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    kappa_set = 1;
    if ((kappa <= 0.0) || (kappa > 2.0))
      lmp->error->all(FLERR, iarg + 1, "Fix {} xlkappa value {} must be > 0 and <= 2", style,
                      kappa);
    return 2;
  }
  if (strcmp(arg[iarg], "xlguard") == 0) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix ") + style + " xlguard", lmp->error);
    guard = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    return 2;
  }
  if (strcmp(arg[iarg], "xldamp") == 0) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix ") + style + " xldamp", lmp->error);
    kdis = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if ((kdis != 0) && ((kdis < KDIS_MIN) || (kdis > KDIS_MAX)))
      lmp->error->all(FLERR, iarg + 1, "Fix {} xldamp value {} must be 0 or between {} and {}",
                      style, kdis, KDIS_MIN, KDIS_MAX);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   must be called after all keywords were parsed and only when enabled
------------------------------------------------------------------------- */

void XLQEq::finalize(LAMMPS *lmp, const char *style, int nevery)
{
  if (nevery != 1)
    lmp->error->all(FLERR, "Fix {} requires Nevery = 1 for extended-Lagrangian charge propagation",
                    style);

  // default is dissipation of order 5: the plain time-reversible propagation
  // accumulates noise in the auxiliary variables and eventually turns unstable.
  // a custom xlkappa implies the undamped scheme (kappa comes from the table
  // otherwise); combining it with a nonzero dissipation order makes no sense

  if (kdis < 0)
    kdis = kappa_set ? 0 : 5;
  else if (kdis && kappa_set)
    lmp->error->all(FLERR, "Fix {} keywords xlkappa and xldamp may not be combined", style);

  if (kdis) {
    const auto &row = dis_table[kdis - KDIS_MIN];
    kappa = row.kappa;
    alpha = row.alpha;
    for (int k = 0; k <= kdis; ++k) c[k] = row.c[k];
  }
  nhist = (kdis + 1 > 2) ? kdis + 1 : 2;

  if (lmp->citeme) {
    lmp->citeme->add(cite_xl_qeq);
    if (kdis) lmp->citeme->add(cite_xl_qeq_dissipation);
  }
}

/* ----------------------------------------------------------------------
   the auxiliary variables can seed the truncated solve only during dynamics,
   never in the setup phase, and only with a completely seeded history from
   consecutive timesteps for the same set of atoms
------------------------------------------------------------------------- */

bool XLQEq::ready(const LAMMPS *lmp) const
{
  return flag && !bypass && (lmp->update->whichflag == 1) && (nseed >= nhist) &&
      (lmp->update->ntimestep == laststep + 1) && (lmp->atom->natoms == natoms);
}

/* ----------------------------------------------------------------------
   a valid history consists of converged solutions from consecutive steps
   for the same set of atoms.  otherwise restart the warm-up by filling all
   slots with the current solution.  atoms created during the run have no
   usable history, and the decision must be the same on all MPI ranks.
------------------------------------------------------------------------- */

bool XLQEq::restart(const LAMMPS *lmp) const
{
  return (nseed == 0) || (lmp->update->ntimestep != laststep + 1) || (lmp->atom->natoms != natoms);
}

/* ---------------------------------------------------------------------- */

void XLQEq::advance(const LAMMPS *lmp, bool active, bool restart)
{
  if (!active) {
    if (restart)
      nseed = 1;
    else
      nseed = (nseed + 1 < nhist) ? nseed + 1 : nhist;
  }
  laststep = lmp->update->ntimestep;
  natoms = lmp->atom->natoms;
}

/* ----------------------------------------------------------------------
   a large residual after the truncated solve indicates that the auxiliary
   variables no longer follow the converged solution.  the caller has
   completed the solve to convergence; the history is seeded again from it
------------------------------------------------------------------------- */

void XLQEq::fallback(LAMMPS *lmp, const char *style, double resid, int warn)
{
  nseed = 0;
  ++nfallback;
  if (warn && (lmp->comm->me == 0))
    lmp->error->warning(FLERR,
                        "Fix {} extended-Lagrangian charge propagation became unstable on step "
                        "{} (residual {:.3g} > xlguard {:.3g}). Restarting it from fully "
                        "converged charges",
                        style, lmp->update->ntimestep, resid, guard);
}
