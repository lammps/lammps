// clang-format off
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
   Dan Bolintineanu (SNL), Joel Clemmer (SNL), Ishan Srivastava (SNL),
   Jeremy Lechman(SNL), Leo Silbert (SNL), Gary Grest (SNL)
----------------------------------------------------------------------- */

#include "pair_granular.h"

#include "atom.h"
#include "comm.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>
#include <vector>

using namespace LAMMPS_NS;
using namespace Granular_NS;
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

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  use_history = 0;
  size_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;
  heat_flag = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>(modify->add_fix("NEIGH_HISTORY_GRANULAR_DUMMY all DUMMY"));
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
    memory->destroy(cutoff_type);
    memory->destroy(types_indices);
    for (int i = 0; i < nmodels; i++) delete models_list[i];
    memory->sfree(models_list);

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
  int i,j,k,ii,jj,inum,jnum,itype,jtype;
  double factor_lj,mi,mj,meff;
  double *forces, *torquesi, *torquesj, dq;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;
  const bool history_update = update->setupflag == 0;

  class GranularModel* model;

  for (int i = 0; i < nmodels; i++)
    models_list[i]->history_update = history_update;

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
  double *heatflow, *temperature;
  if (heat_flag) {
    heatflow = atom->heatflow;
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
      model = models_list[types_indices[itype][jtype]];

      // Reset model and copy initial geometric data
      model->xi = x[i];
      model->xj = x[j];
      model->radi = radius[i];
      model->radj = radius[j];
      if (use_history) model->touch = touch[jj];

      touchflag = model->check_contact();

      if (!touchflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // if any history is needed
      if (use_history) touch[jj] = 1;

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
      model->meff = meff;
      model->vi = v[i];
      model->vj = v[j];
      model->omegai = omega[i];
      model->omegaj = omega[j];
      if (use_history) {
        history = &allhistory[size_history * jj];
        model->history = history;
      }

      if (heat_flag) {
        model->Ti = temperature[i];
        model->Tj = temperature[j];
      }

      model->calculate_forces();

      forces = model->forces;
      torquesi = model->torquesi;
      torquesj = model->torquesj;

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

      if (heat_flag) {
        dq = model->dq;
        heatflow[i] += dq;
        if (force->newton_pair || j < nlocal) heatflow[j] -= dq;
      }

      if (evflag) {
        ev_tally_xyz(i,j,nlocal,force->newton_pair,
          0.0,0.0,forces[0],forces[1],forces[2],model->dx[0],model->dx[1],model->dx[2]);
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
  memory->create(types_indices,n+1,n+1,"pair:types_indices");

  maxmodels = n * n + 1; // should never need any more space
  models_list = (GranularModel **) memory->smalloc(maxmodels * sizeof(GranularModel *), "pair:models_list");
  for (int i = 0; i < maxmodels; i++) models_list[i] = nullptr;
  nmodels = 0;

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
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranular::coeff(int narg, char **arg)
{
  double cutoff_one = -1;

  if (narg < 3)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // Construct new model
  models_list[nmodels] = new GranularModel(lmp);
  class GranularModel* model = models_list[nmodels];
  nmodels += 1;

  //Parse mandatory specification
  int iarg = 2;
  iarg = model->add_sub_model(arg, iarg, narg, NORMAL);

  //Parse optional arguments
  while (iarg < narg) {

    if (strcmp(arg[iarg], "tangential") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
    } else if (strcmp(arg[iarg], "damping") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
    } else if (strcmp(arg[iarg], "rolling") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, ROLLING);
    } else if (strcmp(arg[iarg], "twisting") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, TWISTING);
    } else if (strcmp(arg[iarg], "heat") == 0) {
      iarg = model->add_sub_model(arg, iarg + 1, narg, HEAT);
      heat_flag = 1;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters for cutoff keyword");
      cutoff_one = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      model->limit_damping = 1;
      iarg += 1;
    } else error->all(FLERR, "Illegal pair_coeff command {}", arg[iarg]);
  }

  // Define default damping sub model if unspecified, has no coeffs
  if (!model->damping_model)
    model->construct_sub_model("viscoelastic", DAMPING);
  model->init();

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      types_indices[i][j] = nmodels - 1;
      cutoff_type[i][j] = cutoff_type[j][i] = cutoff_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  // If there are > ntype^2 models, delete unused models
  if (nmodels == maxmodels) prune_models();

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
  init specific to this pair style
------------------------------------------------------------------------- */

void PairGranular::init_style()
{
  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, rmass, omega");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR,"Heat conduction in pair granular requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,"Heat conduction in pair granular requires atom style with heatflow property");
  }

  // allocate history and initialize models
  class GranularModel* model;
  int size_max[NSUBMODELS] = {0};
  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];

    if (model->beyond_contact) {
      beyond_contact = 1;
      use_history = 1; // Need to track if in contact
    }
    if (model->size_history != 0) use_history = 1;

    for (int i = 0; i < NSUBMODELS; i++)
      if (model->sub_models[i]->size_history > size_max[i])
        size_max[i] = model->sub_models[i]->size_history;

    if (model->nondefault_history_transfer) nondefault_history_transfer = 1;
  }

  size_history = 0;
  if (use_history) {
    for (int i = 0; i < NSUBMODELS; i++) size_history += size_max[i];

    // Ensure size history is at least 1 to avoid errors in fix neigh/history
    // This could occur if normal model is beyond_contact but no other entries are required
    // E.g. JKR + linear_nohistory
    size_history = MAX(size_history, 1);
  }

  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];
    int next_index = 0;
    for (int i = 0; i < NSUBMODELS; i++) {
      model->sub_models[i]->history_index = next_index;
      next_index += size_max[i];
    }
  }

  if (use_history) neighbor->add_request(this, NeighConst::REQ_SIZE|NeighConst::REQ_HISTORY);
  else neighbor->add_request(this, NeighConst::REQ_SIZE);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix("NEIGH_HISTORY_GRANULAR_DUMMY",
                                                          "NEIGH_HISTORY_GRANULAR"
                                                          " all NEIGH_HISTORY "
                                                          + std::to_string(size_history),1));
    fix_history->pair = this;
  } else if (use_history) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
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

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranular::init_one(int i, int j)
{
  double cutoff = 0.0;
  class GranularModel* model;

  if (setflag[i][j] == 0) {

    models_list[nmodels] = new GranularModel(lmp);
    types_indices[i][j] = nmodels;
    model = models_list[nmodels];

    nmodels += 1;
    if (nmodels == maxmodels) prune_models();

    class GranularModel* model1 = models_list[types_indices[i][i]];
    class GranularModel* model2 = models_list[types_indices[j][j]];

    int error_code = model->mix_coeffs(model1, model2);
    if (error_code != -1)
      error->all(FLERR,"Granular pair style functional forms are different, "
                 "cannot mix coefficients for types {} and {} \n"
                 "with sub models {} and {}. \n"
                 "This combination must be set explicitly via a "
                 "pair_coeff command",i,j,
                 model1->sub_models[error_code]->name,
                 model2->sub_models[error_code]->name);

    model->init();

    for (int k = 0; k < NSUBMODELS; k++)
      model->sub_models[k]->history_index = model1->sub_models[k]->history_index;

    cutoff_type[i][j] = cutoff_type[j][i] = MAX(cutoff_type[i][i], cutoff_type[j][j]);
  }

  model = models_list[types_indices[i][j]];
  // Check if heat model is defined for all type combinations
  if (heat_flag && !model->heat_model)
    error->all(FLERR, "Must specify a heat model for all pair types");

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  double pulloff;
  if (cutoff_type[i][j] < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
      pulloff = 0.0;
      if (model->beyond_contact) {
        pulloff = model->pulloff_distance(maxrad_dynamic[i], maxrad_dynamic[j]);
        cutoff += pulloff;

        pulloff = model->pulloff_distance(maxrad_frozen[i], maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j] + pulloff);

        pulloff = model->pulloff_distance(maxrad_dynamic[i], maxrad_frozen[j]);
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
  } else if (cutoff_type[i][j] > 0) {
    cutoff = cutoff_type[i][j];
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

  model->dt = update->dt;
  types_indices[j][i] = types_indices[i][j];
  return cutoff;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranular::write_restart(FILE *fp)
{
  int i,j;
  fwrite(&nmodels,sizeof(int),1,fp);
  for (i = 0; i < nmodels; i++) models_list[i]->write_restart(fp);

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cutoff_type[i][j],sizeof(double),1,fp);
        fwrite(&types_indices[i][j],sizeof(int),1,fp);
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

  if (me == 0) utils::sfread(FLERR,&nmodels,sizeof(int),1,fp,nullptr,error);
  MPI_Bcast(&nmodels,1,MPI_INT,0,world);

  for (i = 0; i < nmodels; i++) {
    delete models_list[i];
    models_list[i] = new GranularModel(lmp);
    models_list[i]->read_restart(fp);
    models_list[i]->init();
  }

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cutoff_type[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&types_indices[i][j],sizeof(int),1,fp,nullptr,error);
        }
        MPI_Bcast(&cutoff_type[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&types_indices[i][j],1,MPI_INT,0,world);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  for (int i = 0; i < nmodels; i++) models_list[i]->dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
                            double /*rsq*/, double /* factor_coul */,
                            double factor_lj, double &fforce)
{
  if (factor_lj == 0) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  int nall = atom->nlocal + atom->nghost;
  if ((i >= nall) || (j >= nall))
    error->all(FLERR,"Not enough atoms for pair granular single function");

  class GranularModel* model = models_list[types_indices[itype][jtype]];

  // Reset model and copy initial geometric data
  double **x = atom->x;
  double *radius = atom->radius;

  model->xi = x[i];
  model->xj = x[j];
  model->radi = radius[i];
  model->radj = radius[j];
  model->history_update = 0; // Don't update history

  // If history is needed
  double *history,*allhistory;
  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  if (use_history) {
    if ((fix_history == nullptr) || (fix_history->firstvalue == nullptr))
      error->one(FLERR,"Pair granular single computation needs history");
    allhistory = fix_history->firstvalue[i];
    for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
    }
    history = &allhistory[size_history * neighprev];
    model->touch = fix_history->firstflag[i][neighprev];
  }

  int touchflag = model->check_contact();

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
  double meff = mi * mj / (mi + mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // Copy additional information and calculate forces
  double **v = atom->v;
  double **omega = atom->omega;

  model->meff = meff;
  model->vi = v[i];
  model->vj = v[j];
  model->omegai = omega[i];
  model->omegaj = omega[j];
  model->history = history;

  model->calculate_forces();

  // apply forces & torques
  // Calculate normal component, normalized by r
  fforce = model->Fnormal * model->rinv;

  // set single_extra quantities
  svector[0] = model->fs[0];
  svector[1] = model->fs[1];
  svector[2] = model->fs[2];
  svector[3] = MathExtra::len3(model->fs);
  svector[4] = model->fr[0];
  svector[5] = model->fr[1];
  svector[6] = model->fr[2];
  svector[7] = MathExtra::len3(model->fr);
  svector[8] = model->magtortwist;
  svector[9] = model->dx[0];
  svector[10] = model->dx[1];
  svector[11] = model->dx[2];

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
  class GranularModel* model = models_list[types_indices[itype][jtype]];
  if (model->nondefault_history_transfer) {
    for (int i = 0; i < model->size_history; i++) {
      target[i] = model->transfer_history_factor[i] * source[i];
    }
  } else {
    for (int i = 0; i < model->size_history; i++) {
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
    class GranularModel* model = models_list[types_indices[itype][itype]];
    if (model->beyond_contact) {
      cut += model->pulloff_distance(cut, cut);
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
    class GranularModel* model;
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        model = models_list[types_indices[i][j]];
        if (model->beyond_contact) {
          temp = model->pulloff_distance(r1, r2);
          if (temp > cut) cut = temp;
        }
      }
    }
  }

  cut += r1 + r2;

  return cut;
}

/* ----------------------------------------------------------------------
   remove unused models
------------------------------------------------------------------------- */

void PairGranular::prune_models()
{
  int ntypes = atom->ntypes;
  for (int n = nmodels-1; n >= 0; n--) {

    // Find and delete unused models
    int in_use = 0;
    for (int i = 1; i <= ntypes; i++)
      for (int j = 1; j <= ntypes; j++)
        if (types_indices[i][j] == n) in_use = 1;

    if (in_use) continue;
    delete models_list[n];

    // Shift models if needed
    if (n != nmodels - 1) {
      models_list[n] = models_list[nmodels-1];
      for (int i = 1; i <= ntypes; i++)
        for (int j = 1; j <= ntypes; j++)
          if (types_indices[i][j] == nmodels-1)
            types_indices[i][j] = n;
    }

    models_list[nmodels-1] = nullptr;
    nmodels -= 1;
  }
}
