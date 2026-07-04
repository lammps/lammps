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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#include "pair_ldd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include "ldd_indicator.h"
#include "ldd_potential.h"
#include "ldd_indicator_styles.h"
#include "ldd_potential_styles.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define DOT_PROD_GRAD(a, ta, b, tb) \
  (a[3 * ta] * b[3 * tb] + a[3 * ta + 1] * b[3 * tb + 1] + a[3 * ta + 2] * b[3 * tb + 2])
#define GRADTYPE(a) (3 * a)

namespace {
constexpr int MAXLINE = 1024;

constexpr char KEY_LDD_IND[] = "indicator";
constexpr char KEY_LDD_POTL[] = "potential";
constexpr char KEY_LDD_GRAD[] = "gradient";
constexpr char KEY_LDD_SELF[] = "self";
constexpr char KEY_LDD_IGNORE[] = "ignore";

constexpr char cite_pair_ldd1_c[] = "pair ldd command: doi:10.1063/1.5128665\n\n"
                                    "@Article{DeLyser1,\n"
                                    " author = {Michael R. DeLyser and W. G. Noid},\n"
                                    " title = {Analysis of local density potentials},\n"
                                    " journal = {The journal of chemical physics},\n"
                                    " year =    2019,\n"
                                    " volume =  151,\n"
                                    " pages =   {22:224106}\n"
                                    "}\n\n";

constexpr char cite_pair_ldd2_c[] =
    "pair ldd command gradient keyword: doi:10.1063/5.0075291\n\n"
    "@Article{DeLyser2,\n"
    " author = {Michael R. DeLyser and W. G. Noid},\n"
    " title = {Coarse-grained models for local density gradients},\n"
    " journal = {The Journal of Chemical Physics},\n"
    " year =    2021,\n"
    " volume =  156,\n"
    " pages =   {034106}\n"
    "}\n\n";
}    // namespace

/* ---------------------------------------------------------------------- */

PairLdd::PairLdd(LAMMPS *lmp) :
    Pair(lmp), indicator_map(nullptr), potential_map(nullptr), self_interaction(nullptr),
    ignore_pair(nullptr), ignore_me(nullptr), bGradient(nullptr), Inds(nullptr), Potls(nullptr),
    GradPotls(nullptr), local_density(nullptr), grad_density(nullptr), ld_energy(nullptr),
    ld_grad_energy(nullptr), total_energy(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_ldd1_c);
  if (lmp->citeme) lmp->citeme->add(cite_pair_ldd2_c);

  restartinfo = 0;
  // the force depends on per-atom local densities that are only available
  // after compute(), so there is no meaningful single() pair contribution
  single_enable = 0;
  one_coeff = 1;
  // we communicate the local density and 3 gradient components for each species.
  // the species count is only known after reading the file in coeff(), but
  // nspecies <= ntypes, so size the comm buffers here with the ntypes upper bound
  // (must be set in the constructor so pair_style hybrid picks it up correctly).
  comm_forward = 4 * atom->ntypes;
  comm_reverse = 4 * atom->ntypes;

  nmax = 0;

  map = new int[atom->ntypes + 1];

  LDD_factory();
}

// analogous to void _noopt Force::create_factories()

void PairLdd::LDD_factory()
{

  indicator_map = new IndicatorCreatorMap();

#define LDD_INDICATOR_CLASS
#define LddIndicatorStyle(key, Class) (*indicator_map)[#key] = &indicator_creator<Class>;

#include "ldd_indicator_styles.h"

#undef LddIndicatorStyle
#undef LDD_INDICATOR_CLASS

  potential_map = new PotentialCreatorMap();

#define LDD_POTENTIAL_CLASS
#define LddPotentialStyle(key, Class) (*potential_map)[#key] = &potential_creator<Class>;

#include "ldd_potential_styles.h"

#undef LddPotentialStyle
#undef LDD_POTENTIAL_CLASS
}

/* ---------------------------------------------------------------------- */

PairLdd::~PairLdd()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    // delete the per-species-pair indicator/potential objects and their grids
    if (Inds) {
      for (int i = 0; i < nelements; i++) {
        for (int j = 0; j < nelements; j++) {
          delete Inds[i][j];
          delete Potls[i][j];
          delete GradPotls[i][j];
        }
        delete[] Inds[i];
        delete[] Potls[i];
        delete[] GradPotls[i];
      }
      delete[] Inds;
      delete[] Potls;
      delete[] GradPotls;
    }
    memory->destroy(ignore_pair);
    memory->destroy(ignore_me);
    memory->destroy(bGradient);
    memory->destroy(self_interaction);
    allocated = 0;
  }

  memory->destroy(local_density);
  memory->destroy(grad_density);
  memory->destroy(ld_energy);
  memory->destroy(ld_grad_energy);
  memory->destroy(total_energy);

  delete indicator_map;
  delete potential_map;
}

/* ----------------------------------------------------------------------
   grow per-atom local-density arrays to current atom->nmax
------------------------------------------------------------------------- */

void PairLdd::grow_peratom()
{
  if (atom->nmax <= nmax) return;
  nmax = atom->nmax;
  memory->destroy(local_density);
  memory->destroy(grad_density);
  memory->destroy(ld_energy);
  memory->destroy(ld_grad_energy);
  memory->destroy(total_energy);
  memory->create(local_density, nmax, nelements, "ldd:local_density");
  memory->create(grad_density, nmax, 3 * nelements, "ldd:grad_density");
  memory->create(ld_energy, nmax, nelements, "ldd:ld_energy");
  memory->create(ld_grad_energy, nmax, nelements, "ldd:ld_grad_energy");
  memory->create(total_energy, nmax, "ldd:total_energy");
}

/* ---------------------------------------------------------------------- */

void PairLdd::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  }
  memory->create(cutsq, n + 1, n + 1, "LDD:cutsq");
}

/* ----------------------------------------------------------------------
   allocate the per-species-pair interaction data, sized by nelements
   (called from coeff() once the species count is known from the type map)
------------------------------------------------------------------------- */

void PairLdd::allocate_species()
{
  const int n = nelements;
  memory->create(ignore_me, n, "LDD:ignore_me");
  memory->create(ignore_pair, n, n, "LDD:ignore_pair");
  memory->create(bGradient, n, n, "LDD:bGradient");
  memory->create(self_interaction, n, n, "LDD:self_interaction");

  // 2D grids of owned polymorphic interaction objects (managed with new/delete,
  // not memory->create, since each cell holds a heap-allocated subclass object)
  Inds = new LddIndicator **[n];
  Potls = new LddPotential **[n];
  GradPotls = new LddPotential **[n];

  // default: every species pair is ignored until an entry defines it
  for (int i = 0; i < n; i++) {
    Inds[i] = new LddIndicator *[n];
    Potls[i] = new LddPotential *[n];
    GradPotls[i] = new LddPotential *[n];
    ignore_me[i] = true;
    for (int j = 0; j < n; j++) {
      ignore_pair[i][j] = true;
      bGradient[i][j] = false;
      self_interaction[i][j] = false;
      Inds[i][j] = nullptr;
      Potls[i][j] = nullptr;
      GradPotls[i][j] = nullptr;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLdd::compute(int eflag, int vflag)
{
  // Standard stuff from other force file
  int i, j, ii, jj, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq;
  double fpair_LD;

  evdwl = 0.0;
  ev_init(eflag,
          vflag);    // MCL 09.24.25, this lets per atom energies be set up so we can talk to them

  // grow per-atom local-density arrays to include ghost atoms
  grow_peratom();

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  double fxtmp, fytmp, fztmp;
  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  // ldd stuff
  double r_pair;
  double **const local_dens = local_density;
  double *const LD_ttl_nrg = total_energy;

  // gradient stuff
  // gf1 = grad force 1. used for the part of the force in the ij dir
  // gf2 = grad force 2. used for the part of the force in the ld gradient dir
  double gf1, gf2[3], eij[3];
  double **const grad_dens = grad_density;

  LDD_calculate_LDs();
  LDD_calculate_energies();

  // loop over this processor's atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int iatype = type[i];
    const int itype = map[iatype];    // central-atom species
    // if this species has no ldd potentials, skip this atom
    if (!ignore_me[itype]) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      fxtmp = fytmp = fztmp = 0.0;

      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      // loop over this atom's neighbors
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        const int jatype = type[j];
        jtype = map[jatype];    // neighbor species
        // make sure there's a ld potential for this pair type
        if ((!ignore_pair[itype][jtype]) || (!ignore_pair[jtype][itype])) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cutsq[iatype][jatype]) {
            r_pair = sqrt(rsq);
            // The pair force decomposition for ldd potentials, divided by r_pair
            fpair_LD = 0.0;
            if (!ignore_pair[itype][jtype]) {
              fpair_LD +=
                  Inds[itype][jtype]->wp(r_pair) * Potls[itype][jtype]->f(local_dens[i][jtype]);
            }
            if (!ignore_pair[jtype][itype]) {
              fpair_LD +=
                  Inds[jtype][itype]->wp(r_pair) * Potls[jtype][itype]->f(local_dens[j][itype]);
            }
            fpair_LD /= r_pair;

            f[i][0] += delx * fpair_LD;
            f[i][1] += dely * fpair_LD;
            f[i][2] += delz * fpair_LD;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx * fpair_LD;
              f[j][1] -= dely * fpair_LD;
              f[j][2] -= delz * fpair_LD;
            }
            // set fpair as the ld pair force.
            // needed for virial/pressure calculation later.
            // will be overwritten if there's a gradient force.
            fpair = fpair_LD;
            // check if there's a gradient force for this pair of atom types
            if (bGradient[itype][jtype] || bGradient[jtype][itype]) {
              // calculate unit vector from atom j to atom i
              eij[0] = delx / r_pair;
              eij[1] = dely / r_pair;
              eij[2] = delz / r_pair;
              gf1 = gf2[0] = gf2[1] = gf2[2] = 0.0;

              if (bGradient[itype][jtype]) {
                gf1 += GradPotls[itype][jtype]->f(local_dens[i][jtype]) *
                        DOT_PROD_GRAD(grad_dens[i], jtype, grad_dens[i], jtype) *
                        Inds[itype][jtype]->wp(r_pair) -
                    2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                        (Inds[itype][jtype]->wp2(r_pair) -
                         Inds[itype][jtype]->wp(r_pair) / r_pair) *
                        DOT_PROD_GRAD(eij, 0, grad_dens[i], jtype);

                gf2[0] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype)];
                gf2[1] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 1];

                gf2[2] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 2];
              }
              if (bGradient[jtype][itype]) {
                gf1 += GradPotls[jtype][itype]->f(local_dens[j][itype]) *
                        DOT_PROD_GRAD(grad_dens[j], itype, grad_dens[j], itype) *
                        Inds[jtype][itype]->wp(r_pair) +
                    2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                        (Inds[jtype][itype]->wp2(r_pair) -
                         Inds[jtype][itype]->wp(r_pair) / r_pair) *
                        DOT_PROD_GRAD(eij, 0, grad_dens[j], itype);

                gf2[0] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype)];
                gf2[1] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 1];
                gf2[2] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 2];
              }
              // total gradient force broken down into components
              fxtmp = gf1 * eij[0] + gf2[0];
              fytmp = gf1 * eij[1] + gf2[1];
              fztmp = gf1 * eij[2] + gf2[2];

              f[i][0] += fxtmp;
              f[i][1] += fytmp;
              f[i][2] += fztmp;
              if (newton_pair || j < nlocal) {
                f[j][0] -= fxtmp;
                f[j][1] -= fytmp;
                f[j][2] -= fztmp;
              }
              // add original ld force components
              fxtmp += fpair_LD * delx;
              fytmp += fpair_LD * dely;
              fztmp += fpair_LD * delz;
              // The gradient force contains components not directed along r_ij.
              // Accordingly, we have to call ev_tally_xyz to properly update the virial tensor.
              if (evflag)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fxtmp, fytmp, fztmp, delx, dely,
                             delz);
            }
            // If there's no gradient term, call ev_tally like normal.
            else {
              if (evflag) {
                ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
              }
            }
          }
        }
      }
      // Increment van der Waals' energy with ld potl energy for this atom
      evdwl = LD_ttl_nrg[i];
      if (eflag_atom)    //(eatom != NULL) // We need to talk to the per atom energy like we would've in ev_tally so that compute peratom pe works
      {
        eatom[i] += LD_ttl_nrg[i];
      }
      if (evflag) eng_vdwl += evdwl;
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLdd::settings(int narg, char ** /*arg*/)
{
  // ldd derives its type map and per-atom communication buffer sizes from the
  // number of atom types, which is only known once the simulation box exists.
  // Because ldd is a one_coeff (manybody) style it never reads per-type Pair
  // Coeffs from a data file, so there is no need to define it before the box.
  // Require the box up front and fail clearly instead of corrupting memory.
  // This check belongs in settings() (run only for the pair_style command), not
  // the constructor, so it is not triggered when the style is rebuilt during a
  // read_restart (which restores settings from the file instead of settings()).
  if (domain->box_exist == 0)
    error->all(FLERR, "Pair style ldd must be defined after the simulation box is created");

  // pair_style ldd takes no arguments; the interaction cutoff is the indicator
  // rc of each species pair, read from the potential file in coeff()
  if (narg != 0) error->all(FLERR, "Illegal pair_style ldd command: expected no arguments");
}

/* ---------------------------------------------------------------------- */

// Again, these are analogous to what is done in force.h
LddIndicator *PairLdd::new_indicator(const std::string &wtype)
{
  if (indicator_map->find(wtype) != indicator_map->end()) {
    IndicatorCreator indicator_creator = (*indicator_map)[wtype];
    return indicator_creator(lmp);
  }
  error->all(FLERR, utils::check_packages_for_style("LddIndicator", wtype, lmp));
  return nullptr;
}

template <typename T> LddIndicator *PairLdd::indicator_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

LddPotential *PairLdd::new_potential(const std::string &ptype)
{
  if (potential_map->find(ptype) != potential_map->end()) {
    PotentialCreator potential_creator = (*potential_map)[ptype];
    return potential_creator(lmp);
  }
  error->all(FLERR, utils::check_packages_for_style("LddPotential", ptype, lmp));
  return nullptr;
}

template <typename T> LddPotential *PairLdd::potential_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

void PairLdd::ErrorDoubleKeyword(const char *keyword)
{
  error->all(FLERR, "Found ldd pair_coeff keyword {} twice", keyword);
}

void PairLdd::ErrorNumKeywordArgs(const char *keyword, const char *arglist)
{
  error->all(FLERR, "ldd pair_coeff keyword {} must be followed by: {}", keyword, arglist);
}

/* ---------------------------------------------------------------------- */

void PairLdd::coeff_ldd(int si, int sj, int narg, char **arg)
{
  /* arg holds the keyword portion of a potential file entry (the two leading
   * species names have already been consumed):
   *      indicator wtype r0 rc
   *      self yes/no
   *      potential potl_type *potl_coeffs*
   *      gradient grad_type *grad_coeffs*    (optional)
   *      ignore                              (optional)
   * indicator, self, and potential are required unless the pair is ignored.
   */
  int iarg = 0;
  bool bSelf = false, bIgnore = false;

  // set true once the corresponding keyword is seen (used to catch duplicates)
  bool bkInd = false, bkPotl = false, bkSelf = false, bkGrad = false;

  while (iarg < narg) {
    if (strcmp(arg[iarg], KEY_LDD_IND) == 0) {
      if (bkInd) ErrorDoubleKeyword(KEY_LDD_IND);
      bkInd = true;
      if (iarg + 3 >= narg) ErrorNumKeywordArgs(KEY_LDD_IND, "wtype r0 rc");
      Inds[si][sj] = new_indicator(arg[iarg + 1]);
      Inds[si][sj]->init_coeffs(utils::numeric(FLERR, arg[iarg + 2], false, lmp),
                                utils::numeric(FLERR, arg[iarg + 3], false, lmp),
                                domain->dimension);
      iarg += 4;
    } else if (strcmp(arg[iarg], KEY_LDD_SELF) == 0) {
      if (bkSelf) ErrorDoubleKeyword(KEY_LDD_SELF);
      bkSelf = true;
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_SELF, "yes/no");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        bSelf = true;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        bSelf = false;
      else
        error->all(FLERR, "Expected yes/no after ldd keyword {}, found {}", KEY_LDD_SELF,
                   arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], KEY_LDD_POTL) == 0) {
      if (bkPotl) ErrorDoubleKeyword(KEY_LDD_POTL);
      bkPotl = true;
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_POTL, "type *args*");
      Potls[si][sj] = new_potential(arg[iarg + 1]);
      // pass the whole keyword line so the potential can extract its arguments
      Potls[si][sj]->setup_potl(iarg, narg, arg);
      iarg += (Potls[si][sj]->n_coeffs + 2);
    } else if (strcmp(arg[iarg], KEY_LDD_GRAD) == 0) {
      if (bkGrad) ErrorDoubleKeyword(KEY_LDD_GRAD);
      bkGrad = true;
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_GRAD, "type *args*");
      GradPotls[si][sj] = new_potential(arg[iarg + 1]);
      GradPotls[si][sj]->setup_potl(iarg, narg, arg);
      iarg += (GradPotls[si][sj]->n_coeffs + 2);
    } else if (strcmp(arg[iarg], KEY_LDD_IGNORE) == 0) {
      bIgnore = true;
      iarg += 1;
    } else {
      error->all(FLERR, "Unknown ldd pair_coeff keyword {} (valid: {} {} {} {} {})", arg[iarg],
                 KEY_LDD_IND, KEY_LDD_SELF, KEY_LDD_POTL, KEY_LDD_GRAD, KEY_LDD_IGNORE);
    }
  }

  bGradient[si][sj] = bkGrad;
  ignore_pair[si][sj] = bIgnore;

  if (!bIgnore) {
    if (!bkInd)
      error->all(FLERR, "ldd species pair {} {}: missing required keyword {}", elements[si],
                 elements[sj], KEY_LDD_IND);
    if (!bkPotl)
      error->all(FLERR, "ldd species pair {} {}: missing required keyword {}", elements[si],
                 elements[sj], KEY_LDD_POTL);
    if (!bkSelf)
      error->all(FLERR, "ldd species pair {} {}: missing required keyword {}", elements[si],
                 elements[sj], KEY_LDD_SELF);
    if (Inds[si][sj]->r0 >= Inds[si][sj]->rc)
      error->all(FLERR, "ldd species pair {} {}: r0 ({}) must be less than rc ({})", elements[si],
                 elements[sj], Inds[si][sj]->r0, Inds[si][sj]->rc);

    if (bSelf) {
      if (si == sj)
        self_interaction[si][sj] = true;
      else
        error->warning(FLERR,
                       "ldd self interaction requested for distinct species {} {}; ignoring",
                       elements[si], elements[sj]);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLdd::init_one(int i, int j)
{
  // the cutoff for an atom-type pair is the largest indicator rc over the two
  // ordered species interactions it maps to (both directions use the same,
  // symmetric, neighbor-list cutoff)
  const int si = map[i], sj = map[j];
  double cut = 0.0;
  if (!ignore_pair[si][sj] && Inds[si][sj] && (Inds[si][sj]->rc > cut)) cut = Inds[si][sj]->rc;
  if (!ignore_pair[sj][si] && Inds[sj][si] && (Inds[sj][si]->rc > cut)) cut = Inds[sj][si]->rc;
  return cut;
}

/* ---------------------------------------------------------------------- */

void PairLdd::init_style()    //initialization specific to this pair style
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style ldd requires newton pair on");
  neighbor->request(this, instance_me);
}

/* ----------------------------------------------------------------------
   expose the per-atom local-density data to the rest of LAMMPS (e.g. fix pair)
------------------------------------------------------------------------- */

void *PairLdd::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str, "local_density") == 0) {
    ncol = nelements;
    return (void *) local_density;
  } else if (strcmp(str, "grad_density") == 0) {
    ncol = 3 * nelements;
    return (void *) grad_density;
  } else if (strcmp(str, "energy") == 0) {
    ncol = nelements;
    return (void *) ld_energy;
  } else if (strcmp(str, "grad_energy") == 0) {
    ncol = nelements;
    return (void *) ld_grad_energy;
  } else if (strcmp(str, "total_energy") == 0) {
    ncol = 0;
    return (void *) total_energy;
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLdd::LDD_calculate_LDs()    //
{
  int i, j, m, ii, jj, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq;

  const double *const *const x = atom->x;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  double r_pair, wprime;
  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  double **const local_dens = local_density;
  double **const grad_dens = grad_density;

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) {
      for (int tidx = 0; tidx < nelements; tidx++) {
        local_dens[i][tidx] = 0.0;
        grad_dens[i][GRADTYPE(tidx)] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 1] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 2] = 0.0;
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      for (int tidx = 0; tidx < nelements; tidx++) {
        local_dens[i][tidx] = 0.0;
        grad_dens[i][GRADTYPE(tidx)] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 1] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 2] = 0.0;
      }
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int iatype = type[i];
    const int itype = map[iatype];    // central-atom species
    if (!ignore_me[itype]) {

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      // self interaction
      if (self_interaction[itype][itype]) local_dens[i][itype] += Inds[itype][itype]->invnorm;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        const int jatype = type[j];
        jtype = map[jatype];    // neighbor species
        // MCL We need the displacements if relevent to either interaction
        if ((!ignore_pair[itype][jtype]) || (!ignore_pair[jtype][itype])) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cutsq[iatype][jatype] && (!ignore_pair[itype][jtype])) {
            // Compute local density
            r_pair = sqrt(rsq);
            local_dens[i][jtype] += Inds[itype][jtype]->w(r_pair);
            wprime = Inds[itype][jtype]->wp(r_pair);
            grad_dens[i][GRADTYPE(jtype)] += wprime * delx / r_pair;
            grad_dens[i][GRADTYPE(jtype) + 1] += wprime * dely / r_pair;
            grad_dens[i][GRADTYPE(jtype) + 2] += wprime * delz / r_pair;
          }
          if (newton_pair || j < nlocal) {
            if (!ignore_pair[jtype][itype]) {
              if (rsq < cutsq[jatype][iatype]) {
                r_pair = sqrt(rsq);
                local_dens[j][itype] += Inds[jtype][itype]->w(r_pair);
                wprime = Inds[jtype][itype]->wp(r_pair);
                grad_dens[j][GRADTYPE(itype)] -= wprime * delx / r_pair;
                grad_dens[j][GRADTYPE(itype) + 1] -= wprime * dely / r_pair;
                grad_dens[j][GRADTYPE(itype) + 2] -= wprime * delz / r_pair;
              }
            }
          }
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm(this);    //MCL 4.10.10 on LAMMPS DEV
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void PairLdd::LDD_calculate_energies()
{
  int i, ii;

  const int *const type = atom->type;
  const int inum = list->inum;
  const int *const ilist = list->ilist;

  double **const local_dens = local_density;
  double **const grad_dens = grad_density;
  double **const ld_nrg = ld_energy;
  double **const ld_grad_nrg = ld_grad_energy;
  double *const ld_ttl_nrg = total_energy;
  double Ai2;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int itype = map[type[i]];    // central-atom species
    ld_ttl_nrg[i] = 0;
    if (!ignore_me[itype]) {
      for (int tidx = 0; tidx < nelements; tidx++) {
        ld_nrg[i][tidx] = 0;
        ld_grad_nrg[i][tidx] = 0;
        if (!ignore_pair[itype][tidx]) {
          ld_nrg[i][tidx] = Potls[itype][tidx]->u(local_dens[i][tidx]);
          if (bGradient[itype][tidx]) {
            Ai2 = DOT_PROD_GRAD(grad_dens[i], tidx, grad_dens[i], tidx);
            ld_grad_nrg[i][tidx] = Ai2 * GradPotls[itype][tidx]->u(local_dens[i][tidx]);
          } else {
            ld_grad_nrg[i][tidx] = 0;
          }

          ld_ttl_nrg[i] += ld_nrg[i][tidx] + ld_grad_nrg[i][tidx];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairLdd::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, k, m;
  double **ld = local_density;
  double **ldg = grad_density;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nelements; k++) buf[m++] = ld[j][k];
    for (k = 0; k < nelements; k++) {
      buf[m++] = ldg[j][GRADTYPE(k)];
      buf[m++] = ldg[j][GRADTYPE(k) + 1];
      buf[m++] = ldg[j][GRADTYPE(k) + 2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairLdd::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double **ld = local_density;
  double **ldg = grad_density;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nelements; k++) ld[i][k] = buf[m++];
    for (k = 0; k < nelements; ++k) {
      ldg[i][GRADTYPE(k)] = buf[m++];
      ldg[i][GRADTYPE(k) + 1] = buf[m++];
      ldg[i][GRADTYPE(k) + 2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairLdd::pack_reverse_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double **ld = local_density;
  double **ldg = grad_density;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nelements; k++) buf[m++] = ld[i][k];
    for (k = 0; k < nelements; k++) {
      buf[m++] = ldg[i][GRADTYPE(k)];
      buf[m++] = ldg[i][GRADTYPE(k) + 1];
      buf[m++] = ldg[i][GRADTYPE(k) + 2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairLdd::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, m;
  double **ld = local_density;
  double **ldg = grad_density;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nelements; k++) ld[j][k] += buf[m++];
    for (k = 0; k < nelements; k++) {
      ldg[j][GRADTYPE(k)] += buf[m++];
      ldg[j][GRADTYPE(k) + 1] += buf[m++];
      ldg[j][GRADTYPE(k) + 2] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLdd::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // pair_coeff * * <file> <sp1> <sp2> ... : one species name per atom type, in the
  // manybody style.  The set (and number) of species is taken from these names.
  if (narg != 3 + atom->ntypes)
    error->all(FLERR,
               "Incorrect args for pair_coeff ldd: expected '* * <file> <species> ...' "
               "with one species name per atom type");

  map_element2type(narg - 3, arg + 3);
  for (int i = 1; i <= atom->ntypes; i++)
    if (map[i] < 0)
      error->all(FLERR, "Pair style ldd does not allow a NULL species in the type map");

  allocate_species();
  read_file(arg[2]);
}

/* ----------------------------------------------------------------------
   read the ldd potential file: each entry is

     <species_i> <species_j> indicator <w> r0 rc self yes/no
                 potential <type> <coeffs> [gradient <type> <coeffs>] [ignore]

   the file is read on rank 0 and each line is broadcast so that every rank
   builds the identical (polymorphic) per-species-pair interaction objects.
------------------------------------------------------------------------- */

void PairLdd::read_file(char *filename)
{
  FILE *fp = nullptr;
  if (comm->me == 0) {
    fp = utils::open_potential(filename, lmp, nullptr);
    if (!fp)
      error->one(FLERR, "Cannot open ldd potential file {}: {}", filename, utils::getsyserror());
  }

  // track which ordered species pairs have been defined (duplicate/completeness checks)
  std::vector<int> seen((std::size_t)nelements * nelements, 0);

  char line[MAXLINE];
  int done = 0;
  while (true) {
    if (comm->me == 0)
      if (fgets(line, MAXLINE, fp) == nullptr) done = 1;
    MPI_Bcast(&done, 1, MPI_INT, 0, world);
    if (done) break;
    MPI_Bcast(line, MAXLINE, MPI_CHAR, 0, world);

    // strip comments, then split into whitespace-separated tokens
    if (char *hash = strchr(line, '#')) *hash = '\0';
    std::vector<std::string> words = utils::split_words(line);
    if (words.size() < 2) continue;

    // the first two tokens are species names; skip the entry unless both are mapped
    int si = -1, sj = -1;
    for (int e = 0; e < nelements; e++) {
      if (words[0] == elements[e]) si = e;
      if (words[1] == elements[e]) sj = e;
    }
    if (si < 0 || sj < 0) continue;

    if (seen[si * nelements + sj])
      error->all(FLERR, "Duplicate entry for species pair {} {} in ldd potential file {}", words[0],
                 words[1], filename);
    seen[si * nelements + sj] = 1;

    // hand the keyword portion (everything after the two species names) to the parser
    int nkw = (int) words.size() - 2;
    auto **kw = new char *[nkw];
    for (int k = 0; k < nkw; k++) kw[k] = const_cast<char *>(words[k + 2].c_str());
    coeff_ldd(si, sj, nkw, kw);
    delete[] kw;
  }
  if (comm->me == 0) fclose(fp);

  // every ordered species pair must be specified
  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      if (!seen[i * nelements + j])
        error->all(FLERR, "Missing entry for species pair {} {} in ldd potential file {}",
                   elements[i], elements[j], filename);

  // a species is inactive as a central atom if all of its outgoing pairs are ignored
  for (int i = 0; i < nelements; i++) {
    bool all_ignored = true;
    for (int j = 0; j < nelements; j++)
      if (!ignore_pair[i][j]) all_ignored = false;
    ignore_me[i] = all_ignored;
  }
}
