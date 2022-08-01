// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Dan Bolintineanu (SNL), Ishan Srivastava (SNL), Jeremy Lechman(SNL)
   Leo Silbert (SNL), Gary Grest (SNL)
----------------------------------------------------------------------- */

#include "pair_granular.h"

#include "atom.h"
#include "comm.h"
#include "contact.h"
#include "contact_sub_models.h"
#include "contact_normal_models.h"
#include "contact_tangential_models.h"
#include "contact_damping_models.h"
#include "contact_rolling_models.h"
#include "contact_twisting_models.h"
#include "contact_heat_models.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>
#include <vector>

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

PairGranular::PairGranular(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  finitecutflag = 1;

  single_extra = 12;
  svector = new double[single_extra];

  neighprev = 0;

  nmax = 0;
  mass_rigid = nullptr;

  onerad_dynamic = nullptr;
  onerad_frozen = nullptr;
  maxrad_dynamic = nullptr;
  maxrad_frozen = nullptr;

  dt = update->dt;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  use_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;
  tangential_history_index = 0;
  roll_history_index = 0;
  twist_history_index = 0;
  heat_flag = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>( modify->add_fix("NEIGH_HISTORY_GRANULAR_DUMMY all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairGranular::~PairGranular()
{
  delete[] svector;

  if (!fix_history) modify->delete_fix("NEIGH_HISTORY_GRANULAR_DUMMY");
  else modify->delete_fix("NEIGH_HISTORY_GRANULAR");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }

  memory->destroy(mass_rigid);
}

/* ---------------------------------------------------------------------- */

void PairGranular::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double factor_lj,mi,mj,meff,delx,dely,delz;
  double *forces, *torquesi, *torquesj;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;
  const bool historyupdate = update->setupflag == 0;

  ev_init(eflag,vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *heatflux, *temperature, dq;
  if (heat_flag) {
    heatflux = atom->heatflux;
    temperature = atom->temperature;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      jtype = type[j];

      // Reset model and copy initial geometric data
      models[itype][jtype]->reset_contact();
      models[itype][jtype]->xi = x[i];
      models[itype][jtype]->xj = x[j];
      models[itype][jtype]->radi = radius[i];
      models[itype][jtype]->radj = radius[j];

      touchflag = models[itype][jtype]->check_contact();

      if (!touchflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history*jj];
          for (int k = 0; k < size_history; k++) history[k] = 0.0;
        }
      } else {

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
        models[itype][jtype]->meff = meff;
        models[itype][jtype]->vi = v[i];
        models[itype][jtype]->vj = v[j];
        models[itype][jtype]->omegai = omega[i];
        models[itype][jtype]->omegaj = omega[j];
        models[itype][jtype]->history_update = historyupdate;
        models[itype][jtype]->history = history;
        models[itype][jtype]->prep_contact();

        if (heat_flag) {
          models[itype][jtype]->Ti = temperature[i];
          models[itype][jtype]->Tj = temperature[j];
        }

        if (models[itype][jtype]->beyond_contact) touch[jj] = 1;

        // if any history is needed
        if (use_history) {
          touch[jj] = 1;
          history = &allhistory[size_history*jj];
        }

        models[itype][jtype]->calculate_forces();
        if (heat_flag) dq = models[itype][jtype]->calculate_heat();

        forces = models[itype][jtype]->forces;
        torquesi = models[itype][jtype]->torquesi;
        torquesj = models[itype][jtype]->torquesj;

        // apply forces & torques
        scale3(factor_lj, forces);
        add3(f[i], forces, f[i]);

        scale3(factor_lj, torquesi);
        add3(torque[i], torquesi, torque[i]);
        if (heat_flag) heatflux[i] += dq;

        if (force->newton_pair || j < nlocal) {
          sub3(f[j], forces, f[j]);
          scale3(factor_lj, torquesj);
          add3(torque[j], torquesj, torque[j]);
          if (heat_flag) heatflux[j] -= dq;
        }

        if (evflag) {
          delx = x[i][0] - x[j][0];
          dely = x[i][1] - x[j][1];
          delz = x[i][2] - x[j][2];
          ev_tally_xyz(i,j,nlocal,force->newton_pair,
            0.0,0.0,forces[0],forces[1],forces[2],delx,dely,delz);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranular::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutoff_type,n+1,n+1,"pair:cutoff_type");
  memory->create(models,n+1,n+1,"pair:contact_models");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranular::settings(int narg, char **arg)
{
  if (narg == 1) {
    cutoff_global = utils::numeric(FLERR,arg[0],false,lmp);
  } else {
    cutoff_global = -1; // will be set based on particle sizes, model choice
  }

  normal_history = tangential_history = 0;
  roll_history = twist_history = 0;
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranular::coeff(int narg, char **arg)
{
  double cutoff_one = -1;
  int ncoeff;

  if (narg < 3)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // Create new model stored in vector
  vec_models.push_back(ContactModel());

  //Parse mandatory normal and tangential specifications
  int iarg = 2;
  vec_models.back().init_model(std::string(arg[iarg]), NORMAL);
  ncoeff = vec_models.back().normal_model->num_coeffs;
  iarg += 1;
  if (iarg + ncoeff >= narg)
    error->all(FLERR,"Illegal pair_coeff command"
      "Insufficient arguments provided for normal model.");
  vec_models.back().normal_model->parse_coeffs(arg, iarg);
  iarg += ncoeff;

  if (strcmp(arg[iarg], "tangential") == 0) {
    if (iarg + 1 >= narg)
      error->all(FLERR,"Illegal pair_coeff command, must specify "
          "tangential model after tangential keyword");
    vec_models.back().init_model(std::string(arg[iarg]), TANGENTIAL);
    ncoeff = vec_models.back().tangential_model->num_coeffs;
    iarg += 1;
    if (iarg + ncoeff >= narg)
      error->all(FLERR, "Illegal pair_coeff command"
        "Insufficient arguments provided for tangential model.");
    vec_models.back().tangential_model->parse_coeffs(arg, iarg);
    iarg += ncoeff;
  } else{
    error->all(FLERR, "Illegal pair_coeff command, 'tangential' keyword expected");
  }

  //Parse optional arguments
  while (iarg < narg) {
    if (strcmp(arg[iarg], "damping") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      vec_models.back().init_model(std::string(arg[iarg]), DAMPING);
      ncoeff = vec_models.back().damping_model->num_coeffs;
      iarg += 1;
      if (iarg + ncoeff >= narg)
        error->all(FLERR, "Illegal pair_coeff command"
              "Insufficient arguments provided for damping model.");
      vec_models.back().damping_model->parse_coeffs(arg, iarg);
      iarg += ncoeff;

    } else if (strcmp(arg[iarg], "rolling") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      vec_models.back().init_model(std::string(arg[iarg]), ROLLING);
      ncoeff = vec_models.back().rolling_model->num_coeffs;
      iarg += 1;
      if (iarg + ncoeff >= narg)
        error->all(FLERR, "Illegal pair_coeff command"
              "Insufficient arguments provided for rolling model.");
      vec_models.back().rolling_model->parse_coeffs(arg, iarg);
      iarg += ncoeff;

    } else if (strcmp(arg[iarg], "twisting") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      vec_models.back().init_model(std::string(arg[iarg]), TWISTING);
      ncoeff = vec_models.back().twisting_model->num_coeffs;
      iarg += 1;
      if (iarg + ncoeff >= narg)
        error->all(FLERR, "Illegal pair_coeff command"
              "Insufficient arguments provided for twisting model.");
      vec_models.back().twisting_model->parse_coeffs(arg, iarg);
      iarg += ncoeff;

    } else if (strcmp(arg[iarg], "heat") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      vec_models.back().init_model(std::string(arg[iarg]), HEAT);
      ncoeff = vec_models.back().heat_model->num_coeffs;
      iarg += 1;
      if (iarg + ncoeff >= narg)
        error->all(FLERR, "Illegal pair_coeff command"
              "Insufficient arguments provided for heat model.");
      vec_models.back().heat_model->parse_coeffs(arg, iarg);
      iarg += ncoeff;
      heat_flag = 1;

    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      vec_models.back().cutoff_type = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      vec_models.back().limit_damping = 1;
      iarg += 1;

    } else error->all(FLERR, "Illegal pair_coeff command");
  }

  // Define default damping model if unspecified, takes no args
  if (!vec_models.back().damping_model) {
    vec_models.back().init_model("viscoelastic", DAMPING);
    vec_models.back().damping_model->parse_coeffs(arg, 0);
  }

  if (vec_models.back().limit_damping && !vec_models.back().normal_model->allow_limit_damping)
    error->all(FLERR,"Illegal pair_coeff command, "
        "Cannot limit damping with specified normal contact model");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      models[i][j] = & vec_models.back();
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
  init specific to this pair style
------------------------------------------------------------------------- */

void PairGranular::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, rmass");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  // allocate history and initialize models

  int size_normal_history = 0;
  int size_damping_history = 0;
  int size_tangential_history = 0;
  int size_rolling_history = 0;
  int size_twisting_history = 0;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (models[i][j]->normal_model->size_history != 0 ||
          models[i][j]->damping_model->size_history != 0 ||
          models[i][j]->tangential_model->size_history != 0 ||
          models[i][j]->rolling_model->size_history != 0 ||
          models[i][j]->twisting_model->size_history != 0) use_history = 1;

      if (models[i][j]->normal_model->size_history > size_normal_history)
        size_normal_history = models[i][j]->damping_model->size_history;
      if (models[i][j]->damping_model->size_history > size_damping_history)
        size_damping_history = models[i][j]->normal_model->size_history;
      if (models[i][j]->tangential_model->size_history > size_tangential_history)
        size_tangential_history = models[i][j]->tangential_model->size_history;
      if (models[i][j]->rolling_model->size_history > size_rolling_history)
        size_rolling_history = models[i][j]->rolling_model->size_history;
      if (models[i][j]->twisting_model->size_history > size_twisting_history)
        size_twisting_history = models[i][j]->twisting_model->size_history;
    }
  }

  size_history = size_normal_history + size_damping_history +
      size_tangential_history + size_rolling_history + size_twisting_history;

  damping_history_index = size_normal_history;
  tangential_history_index = size_normal_history + damping_history_index;
  roll_history_index = size_normal_history + damping_history_index + size_tangential_history;
  twist_history_index = size_normal_history + damping_history_index + size_tangential_history + size_rolling_history;

  if (use_history) neighbor->add_request(this, NeighConst::REQ_SIZE|NeighConst::REQ_HISTORY);
  else neighbor->add_request(this, NeighConst::REQ_SIZE);

  dt = update->dt;

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>( modify->replace_fix("NEIGH_HISTORY_GRANULAR_DUMMY",
                                                          "NEIGH_HISTORY_GRANULAR"
                                                          " all NEIGH_HISTORY "
                                                          + std::to_string(size_history),1));
    fix_history->pair = this;
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
      else fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
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

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);

  // set fix which stores history info

  if (size_history > 0) {
    fix_history = dynamic_cast<FixNeighHistory *>( modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranular::init_one(int i, int j)
{
  double cutoff = 0.0;

  if (setflag[i][j] == 0) {
    if ((models[i][i]->normal_model->name != models[j][j]->normal_model->name) ||
        (models[i][i]->damping_model->name != models[j][j]->damping_model->name) ||
        (models[i][i]->tangential_model->name != models[j][j]->tangential_model->name) ||
        (models[i][i]->rolling_model->name != models[j][j]->rolling_model->name) ||
        (models[i][i]->twisting_model->name != models[j][j]->twisting_model->name)) {
      error->all(FLERR,"Granular pair style functional forms are different, "
                 "cannot mix coefficients for types {} and {}. \n"
                 "This combination must be set explicitly via a "
                 "pair_coeff command",i,j);
    }

    vec_models.push_back(ContactModel());
    models[i][j] = models[j][i] = & vec_models.back();
    vec_models.back().init_model(models[i][i]->normal_model->name, NORMAL);
    vec_models.back().init_model(models[i][i]->tangential_model->name, TANGENTIAL);
    vec_models.back().init_model(models[i][i]->damping_model->name, DAMPING);
    vec_models.back().init_model(models[i][i]->rolling_model->name, ROLLING);
    vec_models.back().init_model(models[i][i]->twisting_model->name, TWISTING);
    vec_models.back().init_model(models[i][i]->heat_model->name, HEAT);

    vec_models.back().mix_coeffs(models[i][i], models[j][j]);
  }

  // Check if heat model is defined for all type combinations
  if (heat_flag && !models[i][j]->heat_model)
    error->all(FLERR, "Must specify a heat model for all pair types");

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  double pulloff;
  if (models[i][j]->cutoff_type < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
      pulloff = 0.0;
      if (models[i][j]->beyond_contact) {
        pulloff = models[i][j]->pulloff_distance(maxrad_dynamic[i], maxrad_dynamic[j]);
        cutoff += pulloff;

        pulloff = models[i][j]->pulloff_distance(maxrad_frozen[i], maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j] + pulloff);

        pulloff = models[i][j]->pulloff_distance(maxrad_dynamic[i], maxrad_frozen[j]);
        cutoff = MAX(cutoff,maxrad_dynamic[i] + maxrad_frozen[j] + pulloff);
      }
    } else {

      // radius info about either i or j does not exist
      // (i.e. not present and not about to get poured;
      // set to largest value to not interfere with neighbor list)

      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax,2.0*maxrad_dynamic[k]);
        cutmax = MAX(cutmax,2.0*maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  } else if (models[i][j]->cutoff_type > 0) {
    cutoff = models[i][j]->cutoff_type;
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

  // Copy global options
  models[i][j]->dt = models[j][i]->dt = dt;
  models[i][j]->normal_model->history_index = models[j][i]->normal_model->history_index = normal_history_index;
  models[i][j]->tangential_model->history_index = models[j][i]->tangential_model->history_index = tangential_history_index;
  models[i][j]->rolling_model->history_index = models[j][i]->rolling_model->history_index = roll_history_index;
  models[i][j]->twisting_model->history_index = models[j][i]->twisting_model->history_index = twist_history_index;

  models[i][j]->size_history = models[j][i]->size_history = size_history;
  models[i][j]->init(); // Calculates cumulative properties of sub models
  models[j][i]->init();

  if (models[i][j]->nondefault_history_transfer) nondefault_history_transfer = 1;

  return cutoff;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranular::write_restart(FILE *fp)
{
  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        if (comm->me == 0){
          models[i][j]->write_restart(fp);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranular::read_restart(FILE *fp)
{
  allocate();
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        vec_models.push_back(ContactModel());
        models[i][j] = & vec_models.back();
        models[i][j]->read_restart(fp);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  dt = update->dt;

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      models[i][j]->dt = dt;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
                            double rsq, double /* factor_coul */,
                            double factor_lj, double &fforce)
{

  if (factor_lj == 0) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  int jnum;
  int *jlist;
  double *history,*allhistory;

  int nall = atom->nlocal + atom->nghost;
  if ((i >= nall) || (j >= nall))
    error->all(FLERR,"Not enough atoms for pair granular single function");


  int touchflag;
  double **x = atom->x;
  double *radius = atom->radius;

  // Reset model and copy initial geometric data
  models[itype][jtype]->reset_contact();
  models[itype][jtype]->xi = x[i];
  models[itype][jtype]->xj = x[j];
  models[itype][jtype]->radi = radius[i];
  models[itype][jtype]->radj = radius[j];
  models[i][j]->history_update = 0; // Don't update history

  touchflag = models[itype][jtype]->check_contact();

  if (!touchflag) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle
  double mi, mj, meff;
  double *rmass = atom->rmass;
  int *mask = atom->mask;

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
  double **v = atom->v;
  double **omega = atom->omega;

  models[itype][jtype]->meff = meff;
  models[itype][jtype]->vi = v[i];
  models[itype][jtype]->vj = v[j];
  models[itype][jtype]->omegai = omega[i];
  models[itype][jtype]->omegaj = omega[j];
  models[itype][jtype]->history = history;
  models[itype][jtype]->prep_contact();

  // if any history is needed

  jnum = list->numneigh[i];
  jlist = list->firstneigh[i];

  if (use_history) {
    if ((fix_history == nullptr) || (fix_history->firstvalue == nullptr))
      error->one(FLERR,"Pair granular single computation needs history");
    allhistory = fix_history->firstvalue[i];
    for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
    }
    history = &allhistory[size_history*neighprev];
  }

  double *forces, *torquesi, *torquesj;
  models[itype][jtype]->calculate_forces();
  forces = models[itype][jtype]->forces;
  torquesi = models[itype][jtype]->torquesi;
  torquesj = models[itype][jtype]->torquesj;

  // apply forces & torques

  fforce = MathExtra::len3(forces);

  // set single_extra quantities

  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];

  svector[0] = models[itype][jtype]->fs[0];
  svector[1] = models[itype][jtype]->fs[1];
  svector[2] = models[itype][jtype]->fs[2];
  svector[3] = MathExtra::len3(models[itype][jtype]->fs);
  svector[4] = models[itype][jtype]->fr[0];
  svector[5] = models[itype][jtype]->fr[1];
  svector[6] = models[itype][jtype]->fr[2];
  svector[7] = MathExtra::len3(models[itype][jtype]->fr);
  svector[8] = models[itype][jtype]->magtortwist;
  svector[9] = delx;
  svector[10] = dely;
  svector[11] = delz;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranular::pack_forward_comm(int n, int *list, double *buf,
                                    int /* pbc_flag */, int * /* pbc */)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranular::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGranular::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   transfer history during fix/neigh/history exchange
   only needed if any history entries i-j are not just negative of j-i entries
------------------------------------------------------------------------- */

void PairGranular::transfer_history(double* source, double* target, int itype, int jtype)
{
  if (models[itype][jtype]->nondefault_history_transfer) {
    for (int i = 0; i < size_history; i++) {
      target[i] = models[itype][jtype]->transfer_history_factor[i] * source [i];
    }
  } else {
    for (int i = 0; i < size_history; i++) {
      target[i] = -source[i];
    }
  }
}

/* ----------------------------------------------------------------------
   self-interaction range of particle
------------------------------------------------------------------------- */

double PairGranular::atom2cut(int i)
{
  double cut;

  cut = atom->radius[i] * 2;
  if (beyond_contact) {
    int itype = atom->type[i];
    if (models[itype][itype]->beyond_contact) {
      cut += models[itype][itype]->pulloff_distance(cut, cut);
    }
  }

  return cut;
}

/* ----------------------------------------------------------------------
   maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairGranular::radii2cut(double r1, double r2)
{
  double cut = 0.0;

  if (beyond_contact) {
    int n = atom->ntypes;
    double temp;

    // Check all combinations of i and j to find theoretical maximum pull off distance
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        if (models[i][j]->beyond_contact) {
          temp = models[i][j]->pulloff_distance(r1, r2);
          if (temp > cut) cut = temp;
        }
      }
    }
  }

  cut += r1 + r2;

  return cut;
}
