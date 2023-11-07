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
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "fix_rheo_thermal.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "fix_bond_history.h"
#include "fix_update_special_bonds.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
enum {NONE, CONSTANT, TYPE};

/* ---------------------------------------------------------------------- */

FixRHEOThermal::FixRHEOThermal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), fix_rheo(nullptr), compute_grad(nullptr), compute_vshift(nullptr),
  Tc_type(nullptr), kappa_type(nullptr), cv_type(nullptr), fix_update_special_bonds(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  force_reneighbor = 1;
  next_reneighbor = -1;

  Tc_style = NONE;
  cv_style = NONE;
  conductivity_style = NONE;
  cut_bond = 0;
  comm_forward = 0;

  int ntypes = atom->ntypes;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"conductivity") == 0) {
      // Conductivity arguments
      if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
        conductivity_style = CONSTANT;
        kappa = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        if (kappa < 0.0) error->all(FLERR,"The conductivity must be positive");
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
        conductivity_style = TYPE;
        memory->create(kappa_type,ntypes+1,"rheo_thermal:kappa_type");
        for (int i = 1; i <= ntypes; i++) {
          kappa_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (kappa_type[i] < 0.0) error->all(FLERR,"The conductivity must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg],"specific/heat") == 0) {
      // Cv arguments
      if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for cv option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for cv option");
        cv_style = CONSTANT;
        cv = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        if (cv < 0.0) error->all(FLERR,"The specific heat must be positive");
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for cv option");
        cv_style = TYPE;
        memory->create(cv_type,ntypes + 1,"rheo_thermal:cv_type");
        for (int i = 1; i <= ntypes; i++) {
          cv_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (cv_type[i] < 0.0) error->all(FLERR,"The specific heat must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg],"Tfreeze") == 0) {
      // T freeze arguments
      if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
        Tc_style = CONSTANT;
        Tc = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
        Tc_style = TYPE;
        memory->create(Tc_type, ntypes + 1, "rheo_thermal:Tc_type");
        for (int i = 1; i <= ntypes; i++) {
          Tc_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (Tc_type[i] < 0.0) error->all(FLERR,"The melting temperature must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg],"react") == 0) {
      if (iarg + 2 >= narg) error->all(FLERR, "Insufficient arguments for react option");
      cut_bond = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      btype = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      comm_forward = 1;
      if (cut_bond <= 0.0) error->all(FLERR, "Illegal value for bond lengths");\
      if (btype < 1 || btype > atom->nbondtypes) error->all(FLERR, "Illegal value for bond type");

      cutsq_bond = cut_bond * cut_bond;
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix command, {}", arg[iarg]);
    }
    iarg += 1;
  }

  if (cv_style == NONE || conductivity_style == NONE)
    error->all(FLERR, "Must specify specific heat and conductivity styles\n");
}

/* ---------------------------------------------------------------------- */

FixRHEOThermal::~FixRHEOThermal()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_conductivity", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  memory->destroy(cv_type);
  memory->destroy(Tc_type);
  memory->destroy(kappa_type);
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
  cut_kernel = fix_rheo->h;

  if (cut_bond > cut_kernel)
    error->all(FLERR, "Bonding length exceeds kernel cutoff");

  if (!fix_rheo->thermal_flag)
    error->all(FLERR, "Need to define thermal setting in fix rheo");
  compute_grad = fix_rheo->compute_grad;
  compute_vshift = fix_rheo->compute_vshift;

  dtf = 0.5 * update->dt * force->ftm2v;

  if (atom->temperature_flag != 1)
    error->all(FLERR,"fix rheo/thermal command requires atom property temperature");
  if (atom->heatflow_flag != 1)
    error->all(FLERR,"fix rheo/thermal command requires atom property heatflow");
  if (atom->conductivity_flag != 1)
    error->all(FLERR,"fix rheo/thermal command requires atom property conductivity");


  if (cut_bond > 0.0) {
    if (!force->bond) error->all(FLERR,"Must define a bond style to use reactive bond generation with fix rheo/thermal");
    if (!atom->avec->bonds_allow) error->all(FLERR, "Reactive bond generation in fix rheo/thermal requires atom bonds");

    // all special weights must be 1.0 (no special neighbors) or there must be an instance of fix update/special/bonds
    if (force->special_lj[0] != 1.0 || force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0) {
      auto fixes = modify->get_fix_by_style("UPDATE_SPECIAL_BONDS");
      if (fixes.size() == 0) error->all(FLERR, "Without fix update/special/bonds, reactive bond generation in fix rheo/thermal requires special weights of 1.0");
      fix_update_special_bonds = dynamic_cast<FixUpdateSpecialBonds *>(fixes[0]);
    }

    // need a half neighbor list, built only when particles freeze
    auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    req->set_cutoff(cut_kernel);

    // find instances of bond history to delete data
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

  int *status = atom->status;
  int *mask = atom->mask;
  double *temperature = atom->temperature;
  double **gradt = compute_grad->gradt;
  double **vshift = compute_vshift->array_atom;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (status[i] & STATUS_NO_SHIFT) continue;

    if (mask[i] & groupbit) {
      for (a = 0; a < dim; a++) {
        temperature[i] += dtv * vshift[i][a] * gradt[i][a];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::post_integrate()
{
  int *status = atom->status;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rho = atom->rho;
  int *mask = atom->mask;
  int *type = atom->type;

  double cvi, Tci, Ti;

  int n_melt = 0;
  int n_freeze = 0;

  //Integrate temperature and check status
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_INTEGRATION) continue;

      cvi = calc_cv(i);
      temperature[i] += dtf * heatflow[i] / cvi;

      if (Tc_style != NONE) {
        Ti = temperature[i];
        if (Tc_style == CONSTANT) {
          Tci = Tc;
        } else if (Tc_style == TYPE) {
          Tci = Tc_type[type[i]];
        }

        if (Ti > Tci) {
          // If solid, melt
          if (status[i] & STATUS_SOLID) {
            status[i] &= PHASEMASK;
            status[i] |= STATUS_MELTING;
            n_melt += 1;
          }
        } else {
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
  }

  int n_melt_all, n_freeze_all;
  MPI_Allreduce(&n_melt, &n_melt_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&n_freeze, &n_freeze_all, 1, MPI_INT, MPI_SUM, world);

  if (cut_bond > 0 && (n_melt_all || n_freeze_all)) {
    // Forward status then delete/create bonds
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
  int i;
  int *type = atom->type;
  int *mask = atom->mask;
  double *conductivity = atom->conductivity;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (conductivity_style == CONSTANT) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) conductivity[i] = kappa;
  } else if (conductivity_style == TYPE) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) conductivity[i] = kappa_type[type[i]];
  }
}

/* ----------------------------------------------------------------------
  In the future, update & forward evolving conductivity styles every timestep
------------------------------------------------------------------------- */

void FixRHEOThermal::pre_force(int /*vflag*/)
{
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::final_integrate()
{
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  int *status = atom->status;
  int *mask = atom->mask;

  double cvi;

  //Integrate temperature and check status
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_INTEGRATION) continue;

      cvi = calc_cv(i);
      temperature[i] += dtf * heatflow[i] / cvi;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::break_bonds()
{
  int m, n, nmax, i, j;

  tagint *tag = atom->tag;
  int *status = atom->status;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    for (m = 0; m < num_bond[i]; m++) {
      j = atom->map(bond_atom[i][m]);
      if (!(status[i] & STATUS_MELTING) && !(status[j] & STATUS_MELTING)) continue;

      if (n_histories > 0)
        for (auto &ihistory: histories)
          dynamic_cast<FixBondHistory *>(ihistory)->delete_history(i, num_bond[i] - 1);

      if (fix_update_special_bonds) fix_update_special_bonds->add_broken_bond(i, j);

      // For non-melting neighbors, selectively delete bond if necessary
      if (j >= nlocal || (status[j] & STATUS_MELTING)) continue;
      for (n = 0; n < num_bond[j]; n++) {
        if (bond_atom[j][n] == tag[i]) {
          bond_type[j][n] = 0;
          nmax = num_bond[j] - 1;
          bond_type[j][n] = bond_type[j][nmax];
          bond_atom[j][n] = bond_atom[j][nmax];
          if (n_histories > 0) {
            for (auto &ihistory: histories) {
              dynamic_cast<FixBondHistory *>(ihistory)->shift_history(j, n, nmax);
              dynamic_cast<FixBondHistory *>(ihistory)->delete_history(j, nmax);
            }
          }
          num_bond[j]--;
          break;
        }
      }
    }
    num_bond[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::create_bonds()
{
  int i, j, ii, jj, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int *status = atom->status;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  double **x = atom->x;

  neighbor->build_one(list, 1);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // might be faster to do a full list and just act on the atom that freezes
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(status[i] & STATUS_SOLID)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= SPECIALMASK;

      if (!(status[j] & STATUS_SOLID)) continue;
      if (!(status[i] & STATUS_FREEZING) && !(status[j] & STATUS_FREEZING)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq_bond) continue;

      // Add bonds to owned atoms
      // If newton bond, add to both, otherwise add to whichever has a smaller tag
      if (i < nlocal && (!newton_bond || tag[i] < tag[j])) {
        if (num_bond[i] == atom->bond_per_atom)
          error->one(FLERR,"New bond exceeded bonds per atom in fix rheo/thermal");
        if (fix_update_special_bonds) fix_update_special_bonds->add_created_bond(i, j);
        bond_type[i][num_bond[i]] = btype;
        bond_atom[i][num_bond[i]] = tag[j];
        num_bond[i]++;
      }

      if (j < nlocal && (!newton_bond || tag[j] < tag[i])) {
        if (num_bond[j] == atom->bond_per_atom)
          error->one(FLERR,"New bond exceeded bonds per atom in fix rheo/thermal");
        if (fix_update_special_bonds) fix_update_special_bonds->add_created_bond(i, j);
        bond_type[j][num_bond[j]] = btype;
        bond_atom[j][num_bond[j]] = tag[i];
        num_bond[j]++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_cv(int i)
{
  if (cv_style == CONSTANT) {
    return cv;
  } else if (cv_style == TYPE) {
    return(cv_type[atom->type[i]]);
  }
}


/* ---------------------------------------------------------------------- */

int FixRHEOThermal::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, k, m;
  int *status = atom->status;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(status[j]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  int *status = atom->status;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    status[i] = (int) ubuf(buf[m++]).i;
  }
}
