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
   Contributing author: Jacopo Bilotto (EPFL), Jibril B. Coulibaly
------------------------------------------------------------------------- */

#include "pair_granular_superellipsoid.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_extra.h"
#include "math_extra_superellipsoids.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathExtra;

enum { HOOKE, HERTZ };
enum { MASS_VELOCITY, VISCOELASTIC };
enum { CLASSIC, LINEAR_HISTORY };

static constexpr int NUMSTEP_INITIAL_GUESS = 5;
static constexpr double EPSILON = 1e-10;
static constexpr double MIN_RADIUS_RATIO = 1e-4;
static constexpr double MIN_CURVATURE = 1e-12;

/* ---------------------------------------------------------------------- */

PairGranularSuperellipsoid::PairGranularSuperellipsoid(LAMMPS *lmp) :
    Pair(lmp), onerad_dynamic(nullptr), onerad_frozen(nullptr), maxrad_dynamic(nullptr),
    maxrad_frozen(nullptr), fix_dummy(nullptr), fix_history(nullptr), fix_rigid(nullptr),
    mass_rigid(nullptr), normal_model(nullptr), damping_model(nullptr), tangential_model(nullptr),
    limit_damping(nullptr), kn(nullptr), gamman(nullptr), kt(nullptr), xt(nullptr), xmu(nullptr),
    xi(nullptr), xj(nullptr), vi(nullptr), vj(nullptr), quati(nullptr), quatj(nullptr),
    angmomi(nullptr), angmomj(nullptr), inertiai(nullptr), inertiaj(nullptr), history_data(nullptr),
    xref(nullptr), cutoff_type(nullptr)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  finitecutflag = 1;

  single_extra = 17;
  svector = new double[single_extra];

  // Currently only option, generalize if more added

  neighprev = 0;
  nmax = 0;
  mass_rigid = nullptr;

  onerad_dynamic = nullptr;
  onerad_frozen = nullptr;
  maxrad_dynamic = nullptr;
  maxrad_frozen = nullptr;

  cutoff_type = nullptr;

  limit_damping = nullptr;
  normal_model = nullptr;
  damping_model = nullptr;
  tangential_model = nullptr;

  kn = nullptr;
  gamman = nullptr;
  kt = nullptr;
  xt = nullptr;
  xmu = nullptr;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  default_hist_size = 5;
  size_history = default_hist_size;    // default of 5 values, x0[4] and separating axis

  beyond_contact = 0;
  nondefault_history_transfer = 1;
  heat_flag = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy =
      dynamic_cast<FixDummy *>(modify->add_fix("NEIGH_HISTORY_GRANULAR_SE_DUMMY all DUMMY"));

  contact_formulation = MathExtraSuperellipsoids::FORMULATION_ALGEBRAIC;
}

/* ---------------------------------------------------------------------- */

PairGranularSuperellipsoid::~PairGranularSuperellipsoid()
{
  delete[] svector;

  if (!fix_history)
    modify->delete_fix("NEIGH_HISTORY_GRANULAR_SE_DUMMY");
  else
    modify->delete_fix("NEIGH_HISTORY_GRANULAR_SE");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutoff_type);
    memory->destroy(limit_damping);
    memory->destroy(normal_model);
    memory->destroy(damping_model);
    memory->destroy(tangential_model);
    memory->destroy(kn);
    memory->destroy(gamman);
    memory->destroy(kt);
    memory->destroy(xt);
    memory->destroy(xmu);

    // model variables

    delete[] onerad_dynamic;
    delete[] onerad_frozen;
    delete[] maxrad_dynamic;
    delete[] maxrad_frozen;
  }

  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

void PairGranularSuperellipsoid::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, inum, jnum;
  double factor_lj, mi, mj;

  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *history, *allhistory, **firsthistory;

  bool touchflag = false;
  history_update = update->setupflag == 0;

  ev_init(eflag, vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body", tmp);
    auto *mass_body = (double *) fix_rigid->extract("masstotal", tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid, nmax, "pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0)
        mass_rigid[i] = mass_body[body[i]];
      else
        mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;

  auto *avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::BonusSuper *bonus = avec_ellipsoid->bonus_super;
  int *ellipsoid = atom->ellipsoid;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      jtype = type[j];

      // Reset model and copy initial geometric data

      xi = x[i];
      xj = x[j];
      radi = radius[i];
      radj = radius[j];
      history_data = &allhistory[size_history * jj];
      xref = (tag[i] < tag[j]) ? xi : xj;
      tagi = tag[i];
      tagj = tag[j];
      flagi = bonus[ellipsoid[i]].type;
      flagj = bonus[ellipsoid[j]].type;

      radsum = radi + radj;
      sub3(xi, xj, dx);
      rsq = dot3(dx, dx);

      MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei0);
      MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej0);
      MathExtra::copy2(bonus[ellipsoid[i]].block, blocki0);
      MathExtra::copy2(bonus[ellipsoid[j]].block, blockj0);
      MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei);
      MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej);
      MathExtra::copy2(bonus[ellipsoid[i]].block, blocki);
      MathExtra::copy2(bonus[ellipsoid[j]].block, blockj);
      MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat, Ri);
      MathExtra::quat_to_mat(bonus[ellipsoid[j]].quat, Rj);

      touchjj = touch[jj];

      touchflag = check_contact();

      if (!touchflag) {
        // unset non-touching neighbors
        touch[jj] = 0;
        history = &allhistory[size_history * jj];
        for (k = 0; k < size_history; k++) {
          if (bounding_box && k == 4) continue;    // Do not delete cached axis information
          history[k] = 0.0;
        }
        continue;
      }

      touch[jj] = 1;

      // meff = effective mass of pair of particles
      // if I or J part of rigid body, use body mass
      // if I or J is frozen, meff is other particle
      mi = rmass[i];
      mj = rmass[j];
      if (fix_rigid) {
        if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
        if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
      }
      meff = mi * mj / (mi + mj);
      if (mask[i] & freeze_group_bit) meff = mj;
      if (mask[j] & freeze_group_bit) meff = mi;

      // Copy additional information and prepare force calculations

      vi = v[i];
      vj = v[j];
      angmomi = angmom[i];
      angmomj = angmom[j];
      quati = bonus[ellipsoid[i]].quat;
      quatj = bonus[ellipsoid[j]].quat;
      inertiai = bonus[ellipsoid[i]].inertia;
      inertiaj = bonus[ellipsoid[j]].inertia;

      calculate_forces();

      // apply forces & torques
      scale3(factor_lj, forces);
      add3(f[i], forces, f[i]);

      scale3(factor_lj, torquesi);
      add3(torque[i], torquesi, torque[i]);

      if (force->newton_pair || j < nlocal) {
        sub3(f[j], forces, f[j]);
        scale3(factor_lj, torquesj);
        add3(torque[j], torquesj, torque[j]);
      }

      if (evflag)
        ev_tally_xyz(i, j, nlocal, force->newton_pair, 0.0, 0.0, forces[0], forces[1], forces[2],
                     dx[0], dx[1], dx[2]);    // Correct even for non-spherical particles
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cutoff_type, n + 1, n + 1, "pair:cutoff_type");

  memory->create(limit_damping, n + 1, n + 1, "pair:limit_damping");
  memory->create(normal_model, n + 1, n + 1, "pair:normal_model");
  memory->create(damping_model, n + 1, n + 1, "pair:damping_model");
  memory->create(tangential_model, n + 1, n + 1, "pair:tangential_model");

  memory->create(kn, n + 1, n + 1, "pair:kn");
  memory->create(gamman, n + 1, n + 1, "pair:gamman");
  memory->create(kt, n + 1, n + 1, "pair:kt");
  memory->create(xt, n + 1, n + 1, "pair:xt");
  memory->create(xmu, n + 1, n + 1, "pair:xmu");

  onerad_dynamic = new double[n + 1];
  onerad_frozen = new double[n + 1];
  maxrad_dynamic = new double[n + 1];
  maxrad_frozen = new double[n + 1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::settings(int narg, char **arg)
{
  cutoff_global = -1;    // default: will be set based on particle sizes, model choice
  curvature_model = MathExtraSuperellipsoids::CURV_MEAN;
  bounding_box = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "no_bounding_box") == 0) {
      bounding_box = 0;
      iarg++;
    } else if (strcmp(arg[iarg], "geometric") == 0) {
      contact_formulation = MathExtraSuperellipsoids::FORMULATION_GEOMETRIC;
      iarg++;
    } else if (strcmp(arg[iarg], "curvature_gaussian") == 0) {
      curvature_model = MathExtraSuperellipsoids::CURV_GAUSSIAN;
      iarg++;
    } else if (iarg == 0) {
      // if it is the first argument and not a keyword, assume it is a cutoff
      cutoff_global = utils::numeric(FLERR, arg[iarg], false, lmp);
      iarg++;
    } else
      error->all(FLERR, "Illegal pair_style command");
  }

  if (bounding_box == 0) {
    default_hist_size--;
    size_history--;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::coeff(int narg, char **arg)
{
  double cutoff_one = -1;

  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));

  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int normal_one, damping_one, tangential_one, limit_one;
  double kn_one, gamman_one, kt_one, xt_one, xmu_one;

  int iarg = 2;
  if (strcmp(arg[iarg], "hooke") == 0) {
    normal_one = HOOKE;
    if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "pair granular/superellipsoid", error);
    kn_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    gamman_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
    if (kn_one < 0.0 || gamman_one < 0.0) error->all(FLERR, "Illegal linear normal model");
    iarg += 3;
  } else if (strcmp(arg[iarg], "hertz") == 0) {
    normal_one = HERTZ;
    if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "pair granular/superellipsoid", error);
    kn_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    gamman_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
    if (kn_one < 0.0 || gamman_one < 0.0) error->all(FLERR, "Illegal linear normal model");
    iarg += 3;
  } else {
    error->all(FLERR, "Unknown normal model {}", arg[iarg]);
  }

  damping_one = -1;
  limit_one = 0;

  //Parse optional arguments
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tangential") == 0) {
      iarg++;
      if (strcmp(arg[iarg], "linear_history") == 0) {
        tangential_one = LINEAR_HISTORY;
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "pair granular/superellipsoid", error);
        kt_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        xt_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        xmu_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (kt_one < 0.0 || xt_one < 0.0 || xmu_one < 0.0)
          error->all(FLERR, "Illegal linear tangential model");
        iarg += 4;
      } else if (strcmp(arg[iarg], "classic") == 0) {
        tangential_one = CLASSIC;
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "pair granular/superellipsoid", error);
        kt_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        xt_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        xmu_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (kt_one < 0.0 || xt_one < 0.0 || xmu_one < 0.0)
          error->all(FLERR, "Illegal linear tangential model");
        iarg += 4;
      } else {
        error->all(FLERR, "Unknown tangential model {}", arg[iarg]);
      }
    } else if (strcmp(arg[iarg], "damping") == 0) {
      iarg++;
      if (strcmp(arg[iarg], "mass_velocity") == 0) {
        damping_one = MASS_VELOCITY;
        iarg += 1;
      } else if (strcmp(arg[iarg], "viscoelastic") == 0) {
        damping_one = VISCOELASTIC;
        iarg += 1;
      } else {
        error->all(FLERR, "Unknown normal model {}", arg[iarg]);
      }
    } else if (strcmp(arg[iarg], "rolling") == 0) {
      error->all(FLERR, "Rolling models not yet implemented for superellipsoids");
    } else if (strcmp(arg[iarg], "twisting") == 0) {
      error->all(FLERR, "Twisting models not yet implemented for superellipsoids");
    } else if (strcmp(arg[iarg], "heat") == 0) {
      error->all(FLERR, "Heat models not yet implemented for superellipsoids");
      heat_flag = 1;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters for cutoff keyword");
      cutoff_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      limit_one = 1;
      iarg += 1;
    } else
      error->all(FLERR, "Illegal pair_coeff command {}", arg[iarg]);
  }

  // Define default damping sub model if unspecified, has no coeffs
  if (damping_one == -1) damping_one = VISCOELASTIC;

  // granular model init
  contact_radius_flag = 0;
  if (normal_one == HERTZ || damping_one == VISCOELASTIC) contact_radius_flag = 1;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cutoff_type[i][j] = cutoff_type[j][i] = cutoff_one;
      limit_damping[i][j] = limit_damping[j][i] = limit_one;

      normal_model[i][j] = normal_model[j][i] = normal_one;
      damping_model[i][j] = damping_model[j][i] = damping_one;
      tangential_model[i][j] = tangential_model[j][i] = tangential_one;

      kn[i][j] = kn[j][i] = kn_one;
      gamman[i][j] = gamman[j][i] = gamman_one;

      kt[i][j] = kt[j][i] = kt_one;
      xt[i][j] = xt[j][i] = xt_one;
      xmu[i][j] = xmu[j][i] = xmu_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag || !atom->angmom_flag || !atom->superellipsoid_flag)
    error->all(FLERR,
               "Pair granular/superellipsoid requires atom attributes radius, rmass, "
               "angmom and superellipsoid flag");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Pair granular/superellipsoid requires ghost atoms store velocity");

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR,
                 "Heat conduction in pair granular/superellipsoid requires atom style with "
                 "temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,
                 "Heat conduction in pair granular/superellipsoid requires atom style with "
                 "heatflow property");
  }

  for (i = 0; i < atom->nlocal; i++)
    if (atom->ellipsoid[i] < 0)
      error->one(FLERR, "Pair granular/superellipsoid requires all atoms are ellipsoids");

  // need a granular neighbor list

  neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_HISTORY);

  dt = update->dt;

  // grow history for contact models, right now this is superfluous and is just a placeholder

  int size_history_tangential = 0;
  for (int itype = 1; itype <= atom->ntypes; itype++)
    for (int jtype = 1; jtype <= atom->ntypes; jtype++)
      if (tangential_model[itype][jtype] == CLASSIC ||
          tangential_model[itype][jtype] == LINEAR_HISTORY)
        size_history_tangential = 3;
  size_history += size_history_tangential;

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (fix_history == nullptr) {
    fix_history =
        dynamic_cast<FixNeighHistory *>(modify->replace_fix("NEIGH_HISTORY_GRANULAR_SE_DUMMY",
                                                            "NEIGH_HISTORY_GRANULAR_SE"
                                                            " all NEIGH_HISTORY " +
                                                                std::to_string(size_history),
                                                            1));
    fix_history->pair = this;
  } else {
    fix_history =
        dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR_SE"));
    if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else
        fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (int i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranularSuperellipsoid::init_one(int i, int j)
{
  double cutoff = 0.0;

  if (setflag[i][j] == 0) {

    limit_damping[i][j] = MAX(limit_damping[i][i], limit_damping[j][j]);

    if (normal_model[i][i] != normal_model[j][j] ||
        tangential_model[i][i] != tangential_model[j][j] ||
        damping_model[i][i] != damping_model[j][j])
      error->all(FLERR,
                 "Granular pair style functional forms are different, "
                 "cannot mix coefficients for types {} and {}.\n"
                 "This combination must be set explicitly via a "
                 "pair_coeff command",
                 i, j);

    kn[i][j] = mix_geom(kn[i][i], kn[j][j]);
    gamman[i][j] = mix_geom(gamman[i][i], gamman[j][j]);
    kt[i][j] = mix_geom(kt[i][i], kt[j][j]);
    xt[i][j] = mix_geom(xt[i][i], xt[j][j]);
    xmu[i][j] = mix_geom(xmu[i][i], xmu[j][j]);

    cutoff_type[i][j] = cutoff_type[j][i] = MAX(cutoff_type[i][i], cutoff_type[j][j]);
  }

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  if (cutoff_type[i][j] < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) && (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0) && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
      cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);
      cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
    } else {
      // radius info about either i or j does not exist
      // (i.e. not present and not about to get poured;
      // set to largest value to not interfere with neighbor list)

      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax, 2.0 * maxrad_dynamic[k]);
        cutmax = MAX(cutmax, 2.0 * maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  } else if (cutoff_type[i][j] > 0) {
    cutoff = cutoff_type[i][j];
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

  dt = update->dt;
  return cutoff;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::write_restart(FILE *fp)
{
  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&cutoff_type[i][j], sizeof(double), 1, fp);
        fwrite(&limit_damping[i][j], sizeof(int), 1, fp);
        fwrite(&normal_model[i][j], sizeof(int), 1, fp);
        fwrite(&tangential_model[i][j], sizeof(int), 1, fp);
        fwrite(&damping_model[i][j], sizeof(int), 1, fp);

        fwrite(&kn[i][j], sizeof(double), 1, fp);
        fwrite(&gamman[i][j], sizeof(double), 1, fp);
        fwrite(&kt[i][j], sizeof(double), 1, fp);
        fwrite(&xt[i][j], sizeof(double), 1, fp);
        fwrite(&xmu[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::read_restart(FILE *fp)
{
  allocate();
  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &cutoff_type[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &limit_damping[i][j], sizeof(int), 1, fp, nullptr, error);
          utils::sfread(FLERR, &normal_model[i][j], sizeof(int), 1, fp, nullptr, error);
          utils::sfread(FLERR, &tangential_model[i][j], sizeof(int), 1, fp, nullptr, error);
          utils::sfread(FLERR, &damping_model[i][j], sizeof(int), 1, fp, nullptr, error);

          utils::sfread(FLERR, &kn[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &gamman[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &kt[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &xt[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &xmu[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&cutoff_type[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&limit_damping[i][j], 1, MPI_INT, 0, world);
        MPI_Bcast(&normal_model[i][j], 1, MPI_INT, 0, world);
        MPI_Bcast(&tangential_model[i][j], 1, MPI_INT, 0, world);
        MPI_Bcast(&damping_model[i][j], 1, MPI_INT, 0, world);

        MPI_Bcast(&kn[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&gamman[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&kt[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&xt[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&xmu[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGranularSuperellipsoid::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranularSuperellipsoid::single(int i, int j, int /*itype*/, int /*jtype*/,
                                          double /*rsq*/, double /*factor_coul*/, double factor_lj,
                                          double &fforce)
{
  if (factor_lj == 0) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  int nall = atom->nlocal + atom->nghost;
  if ((i >= nall) || (j >= nall))
    error->all(FLERR, "Not enough atoms for pair granular single function");

  // Reset model and copy initial geometric data

  double *allhistory;
  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];

  if ((fix_history == nullptr) || (fix_history->firstvalue == nullptr))
    error->one(FLERR, "Pair granular single computation needs history");
  allhistory = fix_history->firstvalue[i];
  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }
  touchjj = fix_history->firstflag[i][neighprev];

  xi = atom->x[i];
  xj = atom->x[j];
  radi = atom->radius[i];
  radj = atom->radius[j];
  history_data = &allhistory[size_history * neighprev];
  int indx_ref = (atom->tag[i] < atom->tag[j]) ? i : j;
  xref = atom->x[indx_ref];
  tagi = atom->tag[i];
  tagj = atom->tag[j];
  history_update = 0;    // Don't update history

  auto *avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::BonusSuper *bonus = avec_ellipsoid->bonus_super;
  int *ellipsoid = atom->ellipsoid;

  flagi = bonus[ellipsoid[i]].type;
  flagj = bonus[ellipsoid[j]].type;

  MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei0);
  MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej0);
  MathExtra::copy2(bonus[ellipsoid[i]].block, blocki0);
  MathExtra::copy2(bonus[ellipsoid[j]].block, blockj0);
  MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei);
  MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej);
  MathExtra::copy2(bonus[ellipsoid[i]].block, blocki);
  MathExtra::copy2(bonus[ellipsoid[j]].block, blockj);
  MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat, Ri);
  MathExtra::quat_to_mat(bonus[ellipsoid[j]].quat, Rj);

  int touchflag = check_contact();

  if (!touchflag) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle
  double *rmass = atom->rmass;
  int *mask = atom->mask;

  double mi = rmass[i];
  double mj = rmass[j];
  if (fix_rigid) {
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }
  meff = mi * mj / (mi + mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // Copy additional information and calculate forces

  vi = atom->v[i];
  vj = atom->v[j];
  angmomi = atom->angmom[i];
  angmomj = atom->angmom[j];
  quati = bonus[ellipsoid[i]].quat;
  quatj = bonus[ellipsoid[j]].quat;
  inertiai = bonus[ellipsoid[i]].inertia;
  inertiaj = bonus[ellipsoid[j]].inertia;

  calculate_forces();

  // set single_extra quantities
  svector[0] = fs[0];
  svector[1] = fs[1];
  svector[2] = fs[2];
  svector[3] = MathExtra::len3(fs);
  svector[4] = 0.0;
  svector[5] = 0.0;
  svector[6] = 0.0;
  svector[7] = 0.0;
  svector[8] = 0.0;
  svector[9] = dx[0];
  svector[10] = dx[1];
  svector[11] = dx[2];

  // Superellipsoid specific values - were these included?

  svector[12] = 0.0;    //contact_point_and_Lagrange_multiplier[0]
  svector[13] = 0.0;    //contact_point_and_Lagrange_multiplier[1]
  svector[14] = 0.0;    //contact_point_and_Lagrange_multiplier[2]
  svector[15] = 0.0;    //contact_point_and_Lagrange_multiplier[3]
  svector[16] = 0.0;    //bounding_box_separating_axis_index

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranularSuperellipsoid::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                                  int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranularSuperellipsoid::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   Transfer history
------------------------------------------------------------------------- */

void PairGranularSuperellipsoid::transfer_history(double *source, double *target, int itype,
                                                  int jtype)
{
  // copy of all history variables (shear, contact point, axis)

  for (int i = 0; i < size_history; i++) {
    if (i >= default_hist_size && tangential_model[itype][jtype] == CLASSIC) {
      target[i] = -source[i];    //shear
    } else {
      target[i] = source[i];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGranularSuperellipsoid::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

double PairGranularSuperellipsoid::mix_geom(double val1, double val2)
{
  return sqrt(val1 * val2);
}

/* ---------------------------------------------------------------------- */

double PairGranularSuperellipsoid::mix_mean(double val1, double val2)
{
  return 0.5 * (val1 + val2);
}

/* ---------------------------------------------------------------------- */

bool PairGranularSuperellipsoid::check_contact()
{
  bool touching = false;
  if (rsq >= radsum * radsum) {
    touching = false;
  } else {
    bool skip_contact_detection(false);
    if (bounding_box) {
      int separating_axis = (int) (history_data[4]);
      int new_axis = MathExtraSuperellipsoids::check_oriented_bounding_boxes(
          xi, Ri, shapei, xj, Rj, shapej, separating_axis);
      if (new_axis != -1) {
        skip_contact_detection = true;
        if (history_update) history_data[4] = (double) new_axis;
      }
    }
    if (skip_contact_detection) {
      touching = false;
      return touching;
    }

    double *X0_prev = history_data;

    // superellipsoid contact detection between atoms i and j

    if (touchjj == 1) {
      // Continued contact: use grain true shape and last contact point with respect to grain i
      X0[0] = X0_prev[0] + xref[0];
      X0[1] = X0_prev[1] + xref[1];
      X0[2] = X0_prev[2] + xref[2];
      X0[3] = X0_prev[3];
      int status = MathExtraSuperellipsoids::determine_contact_point(xi, Ri, shapei, blocki, flagi,
                                                                     xj, Rj, shapej, blockj, flagj,
                                                                     X0, nij, contact_formulation);
      if (status == 0) {
        touching = true;
      } else if (status == 1) {
        touching = false;
      } else {
        error->warning(FLERR,
                       "Ellipsoid contact detection (old contact) failed "
                       "between particle {} and particle {} ",
                       tagi, tagj);
      }
    } else {
      // New contact: Build initial guess incrementally by morphing the particles from spheres to actual shape

      // There might be better heuristic for the "volume equivalent spheres" suggested in the paper
      // but this is good enough. We might even be able to use radi and radj which is cheaper
      // MathExtra::scaleadd3(radj / radsum, x[i], radi /radsum, x[j], X0);

      double reqi = std::cbrt(shapei[0] * shapei[1] * shapei[2]);
      double reqj = std::cbrt(shapej[0] * shapej[1] * shapej[2]);
      double rsuminv = 1.0 / (reqi + reqj);
      MathExtra::scaleadd3(reqj * rsuminv, xi, reqi * rsuminv, xj, X0);
      X0[3] = reqj / reqi;    // Lagrange multiplier mu^2
      for (int iter_ig = 1; iter_ig <= NUMSTEP_INITIAL_GUESS; iter_ig++) {
        double frac = iter_ig / double(NUMSTEP_INITIAL_GUESS);
        shapei[0] = shapei[1] = shapei[2] = reqi;
        shapej[0] = shapej[1] = shapej[2] = reqj;
        MathExtra::scaleadd3(1.0 - frac, shapei, frac, shapei0, shapei);
        MathExtra::scaleadd3(1.0 - frac, shapej, frac, shapej0, shapej);
        blocki[0] = 2.0 + frac * (blocki0[0] - 2.0);
        blocki[1] = 2.0 + frac * (blocki0[1] - 2.0);
        blockj[0] = 2.0 + frac * (blockj0[0] - 2.0);
        blockj[1] = 2.0 + frac * (blockj0[1] - 2.0);

        // force ellipsoid flag for first initial guess iteration.
        // Avoid incorrect values of n1/n2 - 2 in second derivatives.
        int status = MathExtraSuperellipsoids::determine_contact_point(
            xi, Ri, shapei, blocki, iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagi,
            xj, Rj, shapej, blockj, iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagj,
            X0, nij, contact_formulation);

        if (status == 0) {
          touching = true;
        } else if (status == 1) {
          touching = false;
        } else if (iter_ig == NUMSTEP_INITIAL_GUESS) {
          // keep trying until last iteration to avoid erroring out too early
          error->warning(FLERR,
                         "Ellipsoid contact detection (new contact) failed "
                         "between particle {} and particle {}",
                         tagi, tagj);
        }
      }
    }
  }

  return touching;
}

/* ---------------------------------------------------------------------- */

void PairGranularSuperellipsoid::calculate_forces()
{
  // Store contact point with respect to grain i for next time step
  // This is crucial for periodic BCs when grains can move by large amount in one time step
  // Keeping the previous contact point relative to global frame would lead to bad initial guess

  if (history_update) {
    double *X0_prev = history_data;
    X0_prev[0] = X0[0] - xref[0];
    X0_prev[1] = X0[1] - xref[1];
    X0_prev[2] = X0[2] - xref[2];
    X0_prev[3] = X0[3];
  }

  double nji[3] = {-nij[0], -nij[1], -nij[2]};
  // compute overlap depth along normal direction for each grain
  // overlap is positive for both grains
  double overlap_i =
      MathExtraSuperellipsoids::compute_overlap_distance(shapei, blocki, Ri, flagi, X0, nij, xi);
  double overlap_j =
      MathExtraSuperellipsoids::compute_overlap_distance(shapej, blockj, Rj, flagj, X0, nji, xj);

  // branch vectors
  double cr_i[3], cr_j[3];
  MathExtra::sub3(X0, xi, cr_i);
  MathExtra::sub3(X0, xj, cr_j);

  // we need to take the cross product of omega

  double ex_space[3], ey_space[3], ez_space[3], omegai[3], omegaj[3];
  MathExtra::q_to_exyz(quati, ex_space, ey_space, ez_space);
  MathExtra::angmom_to_omega(angmomi, ex_space, ey_space, ez_space, inertiai, omegai);
  MathExtra::q_to_exyz(quatj, ex_space, ey_space, ez_space);
  MathExtra::angmom_to_omega(angmomj, ex_space, ey_space, ez_space, inertiaj, omegaj);

  double omega_cross_ri[3], omega_cross_rj[3];
  MathExtra::cross3(omegai, cr_i, omega_cross_ri);
  MathExtra::cross3(omegaj, cr_j, omega_cross_rj);

  // relative translational velocity
  // compute directly the sum of relative translational velocity at contact point
  // since rotational velocity contribution is different for superellipsoids
  double cv_i[3], cv_j[3];
  add3(vi, omega_cross_ri, cv_i);
  add3(vj, omega_cross_rj, cv_j);

  // total relative velocity at contact point
  sub3(cv_i, cv_j, vr);

  // normal component

  vnnr = dot3(vr, nij);
  scale3(vnnr, nij, vn);

  // tangential component

  sub3(vr, vn, vtr);

  vrel = len3(vtr);

  // Approximate contact radius

  // hertzian contact radius approximation
  if (contact_radius_flag) {
    double surf_point_i[3], surf_point_j[3], curvature_i, curvature_j;
    MathExtra::scaleadd3(overlap_i, nij, X0, surf_point_i);
    MathExtra::scaleadd3(overlap_j, nji, X0, surf_point_j);

    if (curvature_model == MathExtraSuperellipsoids::CURV_MEAN) {
      curvature_i = MathExtraSuperellipsoids::mean_curvature_superellipsoid(shapei, blocki, flagi,
                                                                            Ri, surf_point_i, xi);
      curvature_j = MathExtraSuperellipsoids::mean_curvature_superellipsoid(shapej, blockj, flagj,
                                                                            Rj, surf_point_j, xj);
    } else {
      curvature_i = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
          shapei, blocki, flagi, Ri, surf_point_i, xi);
      curvature_j = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
          shapej, blockj, flagj, Rj, surf_point_j, xj);
    }
    double sum_curvature = curvature_i + curvature_j;

    // Physical upper bound smallest particle's bounding sphere radius
    double max_physical_radius = MIN(radi, radj);
    double min_physical_radius = MIN_RADIUS_RATIO * max_physical_radius;

    if (sum_curvature > MIN_CURVATURE) {
      contact_radius = sqrt((overlap_i + overlap_j) / sum_curvature);
      // Cap the maximum radius (flat faces)
      contact_radius = MIN(contact_radius, max_physical_radius);
      // Cap the minimum radius (sharp corners) to prevent force collapse
      contact_radius = MAX(contact_radius, min_physical_radius);
    } else {
      contact_radius = max_physical_radius;
    }

    // hertzian contact radius approximation
    contact_radius = sqrt((overlap_i + overlap_j) / (curvature_i + curvature_j));
  }

  if (normal_model[itype][jtype] == HOOKE) {
    // assuming we get the overlap depth
    Fnormal = kn[itype][jtype] * (overlap_i + overlap_j);
  } else if (normal_model[itype][jtype] == HERTZ) {
    Fnormal = kn[itype][jtype] * (overlap_i + overlap_j) * contact_radius;
  }

  double damp = gamman[itype][jtype];
  double damp_prefactor, Fdamp;
  if (damping_model[itype][jtype] == MASS_VELOCITY) {
    damp_prefactor = damp * meff;
    Fdamp = damp_prefactor * vnnr;
  } else {
    damp_prefactor = damp * meff * contact_radius;
    Fdamp = damp_prefactor * vnnr;
  }

  // normal forces = elastic contact + normal velocity damping

  Fntot = Fnormal + Fdamp;
  if (limit_damping[itype][jtype] && (Fntot < 0.0)) Fntot = 0.0;
  double Fncrit = fabs(Fntot);

  // Tangential model

  double hist_increment[3], fdamp[3];
  double *history = &history_data[default_hist_size];
  double Fscrit = Fncrit * xmu[itype][jtype];
  double dampt = xt[itype][jtype] * damp_prefactor;
  if (tangential_model[itype][jtype] == LINEAR_HISTORY) {
    // rotate and update displacements / force.
    // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235

    int frame_update = 0;
    if (history_update) {
      double rsht = dot3(history, nij);
      frame_update = (fabs(rsht) * kt[itype][jtype]) > (EPSILON * Fscrit);

      if (frame_update) rotate_rescale_vec(history, nij);

      // update history, tangential force using velocities at half step
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(dt, vtr, hist_increment);
      add3(history, hist_increment, history);
    }

    // tangential forces = history + tangential velocity damping
    scale3(-kt[itype][jtype], history, fs);

    scale3(-dampt, vtr, fdamp);
    add3(fs, fdamp, fs);

    // rescale frictional displacements and forces if needed
    double magfs = len3(fs);
    if (magfs > Fscrit) {
      double shrmag = len3(history);
      if (shrmag != 0.0) {
        double magfs_inv = 1.0 / magfs;
        scale3(Fscrit * magfs_inv, fs, history);
        sub3(history, fdamp, history);
        scale3(-1.0 / kt[itype][jtype], history);
        scale3(Fscrit * magfs_inv, fs);
      } else {
        zero3(fs);
      }
    }

  } else if (tangential_model[itype][jtype] == CLASSIC) {

    // shear history effects

    if (history_update) {
      scale3(dt, vtr, hist_increment);
      add3(history, hist_increment, history);
    }
    double shrmag = len3(history);

    if (history_update) {
      // rotate shear displacements
      double rsht = dot3(history, nij);
      scale3(rsht, nij, hist_increment);
      sub3(history, hist_increment, history);
    }

    // tangential forces = history + tangential velocity damping
    if (contact_radius_flag)
      scale3(-kt[itype][jtype] * contact_radius, history, fs);
    else
      scale3(-kt[itype][jtype], history, fs);

    scale3(-dampt, vtr, fdamp);
    add3(fs, fdamp, fs);

    // rescale frictional displacements and forces if needed

    double magfs = len3(fs);

    if (magfs > Fscrit) {
      if (shrmag != 0.0) {
        // Rescale shear force
        scale3(Fscrit / magfs, fs);

        // Set shear to elastic component of rescaled force
        //  has extra factor of kt (+ contact radius)
        sub3(fs, fdamp, history);

        // Remove extra prefactors from shear history
        if (contact_radius_flag)
          scale3(-1.0 / (kt[itype][jtype] * contact_radius), history);
        else
          scale3(-1.0 / kt[itype][jtype], history);
      } else
        zero3(fs);
    }
  }

  // forces & torques

  scale3(Fntot, nji, forces);
  add3(forces, fs, forces);

  cross3(cr_i, forces, torquesi);
  cross3(forces, cr_j, torquesj);
}

/* ----------------------------------------------------------------------
  rotate-rescale vector v so it is perpendicular to unit vector n
  and has the same magnitude as before
    Copied from GranSubMod
  ---------------------------------------------------------------------- */
void PairGranularSuperellipsoid::rotate_rescale_vec(double *v, double *n)
{
  double rsht, shrmag, prjmag, temp_dbl, temp_array[3];

  rsht = dot3(v, n);
  shrmag = len3(v);

  scale3(rsht, n, temp_array);
  sub3(v, temp_array, v);

  // also rescale to preserve magnitude
  prjmag = len3(v);
  if (prjmag > 0)
    temp_dbl = shrmag / prjmag;
  else
    temp_dbl = 0;
  scale3(temp_dbl, v);
}
