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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "fix_rheo_thermal.h"

#include "atom.h"
#include "atom_vec.h"
#include "citeme.h"
#include "comm.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "fix_rheo.h"
#include "fix_update_special_bonds.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
enum { NONE, CONSTANT };

static const char cite_rheo_oxide[] =
    "@article{ApplMathModel.130.310,\n"
    " title = {A hybrid smoothed-particle hydrodynamics model of oxide skins on molten aluminum},\n"
    " journal = {Applied Mathematical Modelling},\n"
    " volume = {130},\n"
    " pages = {310-326},\n"
    " year = {2024},\n"
    " issn = {0307-904X},\n"
    " doi = {https://doi.org/10.1016/j.apm.2024.02.027},\n"
    " author = {Joel T. Clemmer and Flint Pierce and Thomas C. O'Connor and Thomas D. Nevins and "
    "Elizabeth M.C. Jones and Jeremy B. Lechman and John Tencer},\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixRHEOThermal::FixRHEOThermal(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), cv(nullptr), Tc(nullptr), kappa(nullptr), L(nullptr), cv_style(nullptr),
    Tc_style(nullptr), kappa_style(nullptr), L_style(nullptr), list(nullptr), fix_rheo(nullptr),
    compute_grad(nullptr), compute_vshift(nullptr), fix_update_special_bonds(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  cut_bond = 0;
  comm_forward = 0;

  // Currently can only have one instance of fix rheo/thermal
  if (igroup != 0) error->all(FLERR, "fix rheo/thermal command requires group all");

  int i, nlo, nhi;
  int n = atom->ntypes;

  memory->create(Tc_style, n + 1, "rheo:Tc_style");
  memory->create(kappa_style, n + 1, "rheo:kappa_style");
  memory->create(cv_style, n + 1, "rheo:cv_style");
  memory->create(L_style, n + 1, "rheo:L_style");

  memory->create(Tc, n + 1, "rheo:Tc");
  memory->create(kappa, n + 1, "rheo:kappa");
  memory->create(cv, n + 1, "rheo:cv");
  memory->create(L, n + 1, "rheo:L");

  for (i = 1; i <= n; i++) {
    Tc_style[i] = NONE;
    kappa_style[i] = NONE;
    cv_style[i] = NONE;
    L_style[i] = NONE;
  }

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "conductivity") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/thermal conductivity", error);
      utils::bounds(FLERR, arg[iarg + 1], 1, n, nlo, nhi, error);

      // Conductivity arguments
      if (strcmp(arg[iarg + 2], "constant") == 0) {
        if (iarg + 3 >= narg)
          utils::missing_cmd_args(FLERR, "fix rheo/thermal conductivity constant", error);

        double kappa_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (kappa_one < 0.0) error->all(FLERR, "The conductivity must be positive");
        iarg += 2;

        for (i = nlo; i <= nhi; i++) {
          kappa_style[i] = CONSTANT;
          kappa[i] = kappa_one;
        }
      } else {
        error->all(FLERR, "Illegal fix command, {}", arg[iarg + 2]);
      }

      iarg += 2;
    } else if (strcmp(arg[iarg], "specific/heat") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/thermal specific/heat", error);
      utils::bounds(FLERR, arg[iarg + 1], 1, n, nlo, nhi, error);

      // Cv arguments
      if (strcmp(arg[iarg + 2], "constant") == 0) {
        if (iarg + 3 >= narg)
          utils::missing_cmd_args(FLERR, "fix rheo/thermal specific/heat constant", error);

        double cv_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (cv_one < 0.0) error->all(FLERR, "The specific heat must be positive");
        iarg += 2;

        for (i = nlo; i <= nhi; i++) {
          cv_style[i] = CONSTANT;
          cv[i] = cv_one;
        }

      } else {
        error->all(FLERR, "Illegal fix command, {}", arg[iarg + 2]);
      }

      iarg += 2;
    } else if (strcmp(arg[iarg], "Tfreeze") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/thermal Tfreeze", error);
      utils::bounds(FLERR, arg[iarg + 1], 1, n, nlo, nhi, error);

      // T freeze arguments
      if (strcmp(arg[iarg + 2], "constant") == 0) {
        if (iarg + 3 >= narg)
          utils::missing_cmd_args(FLERR, "fix rheo/thermal Tfreeze constant", error);

        double Tc_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        iarg += 2;

        for (i = nlo; i <= nhi; i++) {
          Tc_style[i] = CONSTANT;
          Tc[i] = Tc_one;
        }

      } else {
        error->all(FLERR, "Illegal fix command, {}", arg[iarg + 2]);
      }

      iarg += 2;
    } else if (strcmp(arg[iarg], "latent/heat") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/thermal latent/heat", error);
      utils::bounds(FLERR, arg[iarg + 1], 1, n, nlo, nhi, error);

      // Cv arguments
      if (strcmp(arg[iarg + 2], "constant") == 0) {
        if (iarg + 3 >= narg)
          utils::missing_cmd_args(FLERR, "fix rheo/thermal latent/heat constant", error);

        double L_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (L_one < 0.0) error->all(FLERR, "The latent heat must be positive");
        iarg += 2;

        for (i = nlo; i <= nhi; i++) {
          L_style[i] = CONSTANT;
          L[i] = L_one;
        }

      } else {
        error->all(FLERR, "Illegal fix command, {}", arg[iarg + 2]);
      }

      iarg += 2;
    } else if (strcmp(arg[iarg], "react") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/thermal react", error);
      cut_bond = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      btype = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      comm_forward = 4;
      if (cut_bond <= 0.0) error->all(FLERR, "Illegal max bond length must be greater than zero");
      if ((btype < 1) || (btype > atom->nbondtypes))
        error->all(FLERR, "Illegal value {} for bond type", btype);

      cutsq_bond = cut_bond * cut_bond;
      iarg += 3;
    } else {
      error->all(FLERR, "Unknown fix rheo/thermal keyword: {}", arg[iarg]);
    }
  }

  for (i = 1; i <= n; i++) {
    if (cv_style[i] == NONE)
      error->all(FLERR, "Must specify specific/heat for atom type {} in fix/rheo/thermal", i);
    if (kappa_style[i] == NONE)
      error->all(FLERR, "Must specify conductivity for atom type {} in fix/rheo/thermal", i);
    if (Tc_style[i] == NONE && L_style[i] != NONE)
      error->all(FLERR,
                 "Must specify critical temperature for atom type {} to use latent heat in fix "
                 "rheo/thermal",
                 i);
  }

  if (lmp->citeme) lmp->citeme->add(cite_rheo_oxide);
}

/* ---------------------------------------------------------------------- */

FixRHEOThermal::~FixRHEOThermal()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_conductivity", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  memory->destroy(cv_style);
  memory->destroy(Tc_style);
  memory->destroy(kappa_style);
  memory->destroy(L_style);
  memory->destroy(cv);
  memory->destroy(Tc);
  memory->destroy(kappa);
  memory->destroy(L);
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_NEIGHBOR;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);
  cut_kernel = fix_rheo->cut;

  if (cut_bond > cut_kernel) error->all(FLERR, "Bonding length exceeds kernel cutoff");

  if (!fix_rheo->thermal_flag) error->all(FLERR, "Need to define thermal setting in fix rheo");
  compute_grad = fix_rheo->compute_grad;
  compute_vshift = fix_rheo->compute_vshift;

  dt = update->dt;
  dth = 0.5 * update->dt;

  if (atom->esph_flag != 1)
    error->all(FLERR, "fix rheo/thermal command requires atom property esph");
  if (atom->temperature_flag != 1)
    error->all(FLERR, "fix rheo/thermal command requires atom property temperature");
  if (atom->heatflow_flag != 1)
    error->all(FLERR, "fix rheo/thermal command requires atom property heatflow");
  if (atom->conductivity_flag != 1)
    error->all(FLERR, "fix rheo/thermal command requires atom property conductivity");

  if (cut_bond > 0.0) {
    if (!force->bond)
      error->all(FLERR,
                 "Must define a bond style to use reactive bond generation with fix rheo/thermal");
    if (!atom->avec->bonds_allow)
      error->all(FLERR, "Reactive bond generation in fix rheo/thermal requires atom bonds");

    // all special weights must be 1.0 (no special neighbors) or there must be an instance of fix update/special/bonds
    if (force->special_lj[0] != 1.0 || force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 ||
        force->special_lj[3] != 1.0) {
      auto fixes = modify->get_fix_by_style("UPDATE_SPECIAL_BONDS");
      if (fixes.size() == 0)
        error->all(FLERR,
                   "Without fix update/special/bonds, reactive bond generation in fix rheo/thermal "
                   "requires special weights of 1.0");
      fix_update_special_bonds = dynamic_cast<FixUpdateSpecialBonds *>(fixes[0]);
    }

    // must have newton off so both processors will search nlist to build bonds
    if (force->newton_pair) error->all(FLERR, "Need Newton off for reactive bond generation");

    // need a half neighbor list, built only when particles freeze
    auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    req->set_cutoff(cut_kernel);

    // find instances of bond history to delete/shift data
    histories = modify->get_fix_by_style("BOND_HISTORY");
    n_histories = histories.size();
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::setup_pre_force(int /*vflag*/)
{
  fix_rheo->thermal_fix_defined = 1;

  if (modify->get_fix_by_style("rheo/thermal").size() > 1)
    error->all(FLERR, "More than one fix rheo/thermal defined");

  post_neighbor();
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::initial_integrate(int /*vflag*/)
{
  // update temperature from shifting
  if (!fix_rheo->shift_flag) return;
  int i, a;

  int *status = atom->rheo_status;
  double *energy = atom->esph;
  double **grade = compute_grad->grade;
  double **vshift = compute_vshift->vshift;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (status[i] & STATUS_NO_SHIFT) continue;
    for (a = 0; a < dim; a++) energy[i] += dt * vshift[i][a] * grade[i][a];
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::post_integrate()
{
  int i, itype;
  double cvi, Tci, Ti, Li;

  int *status = atom->rheo_status;
  double *energy = atom->esph;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  int *type = atom->type;

  int n_melt = 0;
  int n_freeze = 0;

  //Integrate energy and check status
  for (i = 0; i < atom->nlocal; i++) {
    if (status[i] & STATUS_NO_INTEGRATION) continue;

    itype = type[i];
    cvi = calc_cv(itype);
    energy[i] += dth * heatflow[i];
    temperature[i] = energy[i] / cvi;

    if (Tc_style[itype] != NONE) {
      Ti = temperature[i];
      Tci = calc_Tc(itype);

      if (L_style[itype] != NONE) {
        Li = calc_L(itype);
        if (Ti > Tci) Ti = MAX(Tci, (energy[i] - Li) / cvi);
        temperature[i] = Ti;
      }

      // Check phase change if Ti != Tci

      if (Ti > Tci) {
        // If solid, melt
        if (status[i] & STATUS_SOLID) {
          status[i] &= PHASEMASK;
          status[i] |= STATUS_MELTING;
          n_melt += 1;
        }
      }

      if (Ti < Tci) {
        // If fluid, freeze
        if (!(status[i] & STATUS_SOLID)) {
          status[i] &= PHASEMASK;
          status[i] |= STATUS_SOLID;
          status[i] |= STATUS_FREEZING;
          n_freeze += 1;
        }
      }
    }
  }

  int n_melt_all, n_freeze_all;
  MPI_Allreduce(&n_melt, &n_melt_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&n_freeze, &n_freeze_all, 1, MPI_INT, MPI_SUM, world);

  if (cut_bond > 0 && (n_melt_all || n_freeze_all)) {

    // If a particle freezes, check if it already has bonds of that type
    // If so, assume it was inserted as a solid particle
    // Note: inserted solid particle may still shift one timestep
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;
    for (i = 0; i < atom->nlocal; i++) {
      if (status[i] & STATUS_FREEZING) {
        for (int n = 0; n < num_bond[i]; n++) {
          if (bond_type[i][n] == btype) {
            status[i] &= ~STATUS_FREEZING;
            break;
          }
        }
      }
    }

    // Forward status + positions (after inititial integrate, before comm)
    comm->forward_comm(this);

    if (n_freeze_all) create_bonds();
    if (n_melt_all) break_bonds();

    next_reneighbor = update->ntimestep;
  }
}

/* ----------------------------------------------------------------------
  Only need to update non-evolving conductivity styles after atoms exchange
------------------------------------------------------------------------- */

void FixRHEOThermal::post_neighbor()
{
  int i, itype;
  int *type = atom->type;
  double *conductivity = atom->conductivity;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    itype = type[i];
    if (kappa_style[itype] == CONSTANT) conductivity[i] = kappa[itype];
  }
}

/* ----------------------------------------------------------------------
  Calculate temperature
  In the future, update & forward evolving conductivity styles every timestep
------------------------------------------------------------------------- */

void FixRHEOThermal::pre_force(int /*vflag*/)
{
  double cvi, Tci, Ti, Li;

  double *energy = atom->esph;
  double *temperature = atom->temperature;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;

  // Calculate temperature
  for (int i = 0; i < nall; i++) {
    int itype = type[i];
    cvi = calc_cv(itype);
    temperature[i] = energy[i] / cvi;

    if (Tc_style[itype] != NONE) {
      Ti = temperature[i];
      Tci = calc_Tc(itype);

      if (L_style[itype] != NONE) {
        Li = calc_L(itype);
        if (Ti > Tci) Ti = MAX(Tci, (energy[i] - Li) / cvi);
        temperature[i] = Ti;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::final_integrate()
{
  int *status = atom->rheo_status;
  double *energy = atom->esph;
  double *heatflow = atom->heatflow;

  //Integrate energy
  for (int i = 0; i < atom->nlocal; i++) {
    if (status[i] & STATUS_NO_INTEGRATION) continue;
    energy[i] += dth * heatflow[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::reset_dt()
{
  dt = update->dt;
  dth = 0.5 * update->dt;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::break_bonds()
{
  int m, n, nmax, i, j, melti, meltj;

  tagint *tag = atom->tag;
  int *status = atom->rheo_status;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int nlocal = atom->nlocal;

  // Delete all bonds for local atoms that melt of a given type
  for (int i = 0; i < nlocal; i++) {
    melti = status[i] & STATUS_MELTING;
    if (!melti) continue;
    for (m = (num_bond[i] - 1); m >= 0; m--) {
      if (bond_type[i][m] != btype) continue;

      j = atom->map(bond_atom[i][m]);
      meltj = status[j] & STATUS_MELTING;

      nmax = num_bond[i] - 1;
      if (m == nmax) {
        if (n_histories > 0)
          for (auto &ihistory : histories)
            dynamic_cast<FixBondHistory *>(ihistory)->delete_history(i, m);
      } else {
        bond_type[i][m] = bond_type[i][nmax];
        bond_atom[i][m] = bond_atom[i][nmax];
        if (n_histories > 0) {
          for (auto &ihistory : histories) {
            auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
            fix_bond_history->shift_history(i, m, nmax);
            fix_bond_history->delete_history(i, nmax);
          }
        }
      }
      bond_type[i][nmax] = 0;
      num_bond[i]--;

      // Update special unless two owned atoms melt simultaneously then
      //  only update for atom with lower tag
      if (fix_update_special_bonds) {
        if ((i < nlocal) && (j < nlocal) && melti && meltj) {
          if (tag[i] < tag[j]) { fix_update_special_bonds->add_broken_bond(i, j); }
        } else {
          fix_update_special_bonds->add_broken_bond(i, j);
        }
      }
    }
  }

  // Update bond list and break solid-melted bonds
  for (n = 0; n < nbondlist; n++) {

    // skip bond if not correct type
    if (bondlist[n][2] != btype) continue;
    i = bondlist[n][0];
    j = bondlist[n][1];

    melti = status[i] & STATUS_MELTING;
    meltj = status[j] & STATUS_MELTING;

    if (!melti && !meltj) continue;

    bondlist[n][2] = 0;

    // Delete bonds for non-melted local atoms (shifting)
    if (i < nlocal && !melti) {
      for (m = 0; m < num_bond[i]; m++) {
        if ((bond_atom[i][m] == tag[j]) && (bond_type[i][m] == btype)) {
          nmax = num_bond[i] - 1;
          bond_type[i][m] = bond_type[i][nmax];
          bond_atom[i][m] = bond_atom[i][nmax];
          if (n_histories > 0)
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->shift_history(i, m, nmax);
              fix_bond_history->delete_history(i, nmax);
            }
          bond_type[i][nmax] = 0;
          num_bond[i]--;
          break;
        }
      }
    }

    if (j < nlocal && !meltj) {
      for (m = 0; m < num_bond[j]; m++) {
        if ((bond_atom[j][m] == tag[i]) && (bond_type[j][m] == btype)) {
          nmax = num_bond[j] - 1;
          bond_type[j][m] = bond_type[j][nmax];
          bond_atom[j][m] = bond_atom[j][nmax];
          if (n_histories > 0)
            for (auto &ihistory : histories) {
              auto fix_bond_history = dynamic_cast<FixBondHistory *>(ihistory);
              fix_bond_history->shift_history(j, m, nmax);
              fix_bond_history->delete_history(j, nmax);
            }
          bond_type[j][nmax] = 0;
          num_bond[j]--;
          break;
        }
      }
    }

    // Unless both atoms melt simultaneously, need to remove special bond if the melted atom is a ghost
    if (melti && meltj) continue;
    if (fix_update_special_bonds)
      if (((i >= nlocal) && melti) || ((j >= nlocal) && meltj))
        fix_update_special_bonds->add_broken_bond(i, j);
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::create_bonds()
{
  int i, j, ii, jj, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double delx, dely, delz, rsq;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int *status = atom->rheo_status;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  double **x = atom->x;

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // might be faster to do a full list and just act on the atom that freezes
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(status[i] & STATUS_SOLID)) continue;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(status[j] & STATUS_SOLID)) continue;
      if (!(status[i] & STATUS_FREEZING) && !(status[j] & STATUS_FREEZING)) continue;

      // Ensure pair is always ordered to ensure numerical operations
      // are identical to minimize the possibility that a bond straddling
      // an mpi grid (newton off) isn't created on one proc but not the other
      if (tag[i] < tag[j]) {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
      } else {
        delx = x[j][0] - x[i][0];
        dely = x[j][1] - x[i][1];
        delz = x[j][2] - x[i][2];
      }
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq_bond) continue;

      // Add bonds to owned atoms
      // If newton bond off, add to both, otherwise add to whichever has a smaller tag
      if ((i < nlocal) && (!newton_bond || (tag[i] < tag[j]))) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom in fix rheo/thermal");
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tag[j];
        num_bond[i]++;
      }

      if ((j < nlocal) && (!newton_bond || (tag[j] < tag[i]))) {
        if (num_bond[j] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom in fix rheo/thermal");
        bond_type[j][num_bond[j]] = btype;
        bond_atom[j][num_bond[j]] = tag[i];
        num_bond[j]++;
      }

      if (fix_update_special_bonds) fix_update_special_bonds->add_created_bond(i, j);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_cv(int itype)
{
  if (cv_style[itype] == CONSTANT) { return cv[itype]; }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_Tc(int itype)
{
  if (Tc_style[itype] == CONSTANT) { return Tc[itype]; }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_L(int itype)
{
  if (L_style[itype] == CONSTANT) { return L[itype]; }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                      int * /*pbc*/)
{
  int *status = atom->rheo_status;
  double **x = atom->x;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = ubuf(status[j]).d;
    buf[m++] = x[j][0];
    buf[m++] = x[j][1];
    buf[m++] = x[j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::unpack_forward_comm(int n, int first, double *buf)
{
  int *status = atom->rheo_status;
  double **x = atom->x;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    status[i] = (int) ubuf(buf[m++]).i;
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}
