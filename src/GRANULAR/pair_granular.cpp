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
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;
using namespace Contact;



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

  models = nullptr;

  limit_damping = nullptr;

  history_transfer_factors = nullptr;

  dt = update->dt;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  use_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;
  tangential_history_index = 0;
  roll_history_index = twist_history_index = 0;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>( modify->add_fix("NEIGH_HISTORY_GRANULAR_DUMMY all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairGranular::~PairGranular()
{
  delete[] svector;
  delete[] history_transfer_factors;

  if (!fix_history) modify->delete_fix("NEIGH_HISTORY_GRANULAR_DUMMY");
  else modify->delete_fix("NEIGH_HISTORY_GRANULAR");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutoff_type);

    memory->destroy(models);

    memory->destroy(normal_coeffs);
    memory->destroy(tangential_coeffs);
    memory->destroy(roll_coeffs);
    memory->destroy(twist_coeffs);

    memory->destroy(Emod);
    memory->destroy(poiss);

    memory->destroy(normal_model);
    memory->destroy(damping_model);
    memory->destroy(tangential_model);
    memory->destroy(roll_model);
    memory->destroy(twist_model);
    memory->destroy(limit_damping);

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
  double factor_lj,mi,mj,meff;
  double forces[3], torquesi[3], torquesj[3];

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
      models[itype][jtype].xi = x[i];
      models[itype][jtype].xj = x[j];
      models[itype][jtype].radi = radius[i];
      models[itype][jtype].radj = radius[j];


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
        models[itype][jtype].meff = meff;
        models[itype][jtype].dt = dt;
        models[itype][jtype].history_update = historyupdate;
        models[itype][jtype].roll_history_index = roll_history_index;
        models[itype][jtype].twist_history_index = twist_history_index;
        models[itype][jtype].vi = v[i];
        models[itype][jtype].vj = v[j];
        models[itype][jtype].omegai = omega[i];
        models[itype][jtype].omegaj = omega[j];
        models[itype][jtype] -> prep_contact();


        if (models[itype][jtype].normal_model == JKR) touch[jj] = 1;
        // if any history is needed
        if (use_history) {
          touch[jj] = 1;
          history = &allhistory[size_history*jj];
        }

        models[itype][jtype].calculate_forces(forces, torquesi, torquesj, history);

        // apply forces & torques
        MathExtra::scale3(factor_lj, forces);
        MathExtra::add3(f[i], forces, f[i]);

        MathExtra::scale3(factor_lj, torquesi);
        MathExtra::add3(torques[i], torquesi, torques[i]);

        if (force->newton_pair || j < nlocal) {
          MathExtra::sub3(f[j], forces, f[j]);
          MathExtra::scale3(factor_lj, torquesj);
          MathExtra::add3(torques[j], torquesj, torques[j]);
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,force->newton_pair,
            0.0,0.0,forces[0],forces[1],forces[2],dx[0],dy[1],dx[2]);
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
  memory->create(models,n+1,n+1,"pair:models");

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
  int normal_model_one, damping_model_one;
  int tangential_model_one, roll_model_one, twist_model_one;

  double normal_coeffs_one[4];
  double tangential_coeffs_one[4];
  double roll_coeffs_one[4];
  double twist_coeffs_one[4];

  double cutoff_one = -1;

  if (narg < 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  //Defaults
  normal_model_one = tangential_model_one = -1;
  roll_model_one = ROLL_NONE;
  twist_model_one = TWIST_NONE;
  damping_model_one = VISCOELASTIC;
  int ld_flag = 0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "hooke") == 0) {
      if (iarg + 2 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hooke option");
      normal_model_one = HOOKE;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // kn
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      iarg += 3;
    } else if (strcmp(arg[iarg], "hertz") == 0) {
      if (iarg + 2 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz option");
      normal_model_one = HERTZ;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // kn
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      iarg += 3;
    } else if (strcmp(arg[iarg], "hertz/material") == 0) {
      if (iarg + 3 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz/material option");
      normal_model_one = HERTZ_MATERIAL;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      iarg += 4;
    } else if (strcmp(arg[iarg], "dmt") == 0) {
      if (iarg + 4 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for Hertz option");
      normal_model_one = DMT;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      normal_coeffs_one[3] = utils::numeric(FLERR,arg[iarg+4],false,lmp); // cohesion
      iarg += 5;
    } else if (strcmp(arg[iarg], "jkr") == 0) {
      if (iarg + 4 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, "
                   "not enough parameters provided for JKR option");
      beyond_contact = 1;
      normal_model_one = JKR;
      normal_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp); // E
      normal_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // damping
      normal_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp); // Poisson's ratio
      normal_coeffs_one[3] = utils::numeric(FLERR,arg[iarg+4],false,lmp); // cohesion
      iarg += 5;
    } else if (strcmp(arg[iarg], "damping") == 0) {
      if (iarg+1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, "
                   "not enough parameters provided for damping model");
      if (strcmp(arg[iarg+1], "velocity") == 0) {
        damping_model_one = VELOCITY;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "mass_velocity") == 0) {
        damping_model_one = MASS_VELOCITY;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "viscoelastic") == 0) {
        damping_model_one = VISCOELASTIC;
        iarg += 1;
      } else if (strcmp(arg[iarg+1], "tsuji") == 0) {
        damping_model_one = TSUJI;
        iarg += 1;
      } else error->all(FLERR, "Illegal pair_coeff command, "
                        "unrecognized damping model");
      iarg += 1;
    } else if (strcmp(arg[iarg], "tangential") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR,"Illegal pair_coeff command, must specify "
                   "tangential model after tangential keyword");
      if (strcmp(arg[iarg+1], "linear_nohistory") == 0) {
        if (iarg + 3 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for tangential model");
        tangential_model_one = TANGENTIAL_NOHISTORY;
        tangential_coeffs_one[0] = 0;
        // gammat and friction coeff
        tangential_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        tangential_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        iarg += 4;
      } else if ((strcmp(arg[iarg+1], "linear_history") == 0) ||
               (strcmp(arg[iarg+1], "mindlin") == 0) ||
               (strcmp(arg[iarg+1], "mindlin_rescale") == 0) ||
               (strcmp(arg[iarg+1], "mindlin/force") == 0) ||
               (strcmp(arg[iarg+1], "mindlin_rescale/force") == 0)) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for tangential model");
        if (strcmp(arg[iarg+1], "linear_history") == 0)
          tangential_model_one = TANGENTIAL_HISTORY;
        else if (strcmp(arg[iarg+1], "mindlin") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN;
        else if (strcmp(arg[iarg+1], "mindlin_rescale") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_RESCALE;
        else if (strcmp(arg[iarg+1], "mindlin/force") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_FORCE;
        else if (strcmp(arg[iarg+1], "mindlin_rescale/force") == 0)
          tangential_model_one = TANGENTIAL_MINDLIN_RESCALE_FORCE;
        tangential_history = 1;
        if ((tangential_model_one == TANGENTIAL_MINDLIN ||
             tangential_model_one == TANGENTIAL_MINDLIN_RESCALE ||
             tangential_model_one == TANGENTIAL_MINDLIN_FORCE ||
             tangential_model_one == TANGENTIAL_MINDLIN_RESCALE_FORCE) &&
            (strcmp(arg[iarg+2], "NULL") == 0)) {
          if (normal_model_one == HERTZ || normal_model_one == HOOKE) {
            error->all(FLERR, "NULL setting for Mindlin tangential "
                       "stiffness requires a normal contact model that "
                       "specifies material properties");
          }
          tangential_coeffs_one[0] = -1;
        } else {
          tangential_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp); // kt
        }
        // gammat and friction coeff
        tangential_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        tangential_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "tangential model not recognized");
      }
    } else if (strcmp(arg[iarg], "rolling") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0) {
        roll_model_one = ROLL_NONE;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "sds") == 0) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for rolling model");
        roll_model_one = ROLL_SDS;
        roll_history = 1;
        // kR and gammaR and rolling friction coeff
        roll_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        roll_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        roll_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "rolling friction model not recognized");
      }
    } else if (strcmp(arg[iarg], "twisting") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0) {
        twist_model_one = TWIST_NONE;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "marshall") == 0) {
        twist_model_one = TWIST_MARSHALL;
        twist_history = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg+1], "sds") == 0) {
        if (iarg + 4 >= narg)
          error->all(FLERR,"Illegal pair_coeff command, "
                     "not enough parameters provided for twist model");
        twist_model_one = TWIST_SDS;
        twist_history = 1;
        // kt and gammat and friction coeff
        twist_coeffs_one[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        twist_coeffs_one[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        twist_coeffs_one[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      } else {
        error->all(FLERR, "Illegal pair_coeff command, "
                   "twisting friction model not recognized");
      }
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      cutoff_one = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      ld_flag = 1;
      iarg += 1;
    } else error->all(FLERR, "Illegal pair_coeff command");
  }

  // error not to specify normal or tangential model
  if ((normal_model_one < 0) || (tangential_model_one < 0))
    error->all(FLERR, "Illegal pair_coeff command, "
               "must specify normal or tangential contact model");

  int count = 0;
  double damp;
  if (damping_model_one == TSUJI) {
    double cor;
    cor = normal_coeffs_one[1];
    damp = 1.2728-4.2783*cor+11.087*square(cor)-22.348*cube(cor)+
        27.467*powint(cor,4)-18.022*powint(cor,5)+4.8218*powint(cor,6);
  } else damp = normal_coeffs_one[1];

  if (ld_flag && normal_model_one == JKR)
    error->all(FLERR,"Illegal pair_coeff command, "
        "Cannot limit damping with JKR model");

  if (ld_flag && normal_model_one == DMT)
    error->all(FLERR,"Illegal pair_coeff command, "
        "Cannot limit damping with DMT model");

  double Emod, poisson;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      // Define normal model
      models[i][j].normal_model = models[j][i].normal_model = normal_model_one;
      models[i][j].damping_model = models[j][i].damping_model = damping_model_one;
      models[i][j].gamma_norm = models[j][i].gamma_norm = damp;
      Emod = normal_coeffs_one[0];
      poisson = normal_coeffs_one[2];
      models[i][j].Emod = models[j][i].Emod = Emod;
      models[i][j].poisson = models[j][i].poisson = poisson;

      if (normal_model_one != HERTZ && normal_model_one != HOOKE) {
        models[i][j].k_norm = models[j][i].k_norm =
          FOURTHIRDS*mix_stiffnessE(Emod,Emod,poisson,poisson);
      } else {
        models[i][j].k_norm = models[j][i].k_norm = normal_coeffs_one[0];
      }
      if ((normal_model_one == JKR) || (normal_model_one == DMT))
        models[i][j].cohesion = models[j][i].cohesion = normal_coeffs_one[3];

      // Define tangential model
      models[i][j].tangential_model = models[j][i].tangential_model = tangential_model_one;
      if (tangential_coeffs_one[0] == -1) {
        models[i][j].k_tang = models[j][i].k_tang =
          8*mix_stiffnessG(Emod, Emod, poisson, poisson);
      } else {
        models[i][j].k_tang = models[j][i].k_tang = tangential_coeffs_one[0];
      }
      models[i][j].gamma_tang = models[j][i]].gamma_tang = tangential_coeffs_one[1];
      models[i][j].mu_tang = models[j][i]].mu_tang = tangential_coeffs_one[2];

      // Define rolling model
      model[i][j].roll_model = model[j][i].roll_model = roll_model_one;
      if (roll_model_one != ROLL_NONE) {
        model[i][j].k_roll     = model[j][i].k_roll     = roll_coeffs_one[0];
        model[i][j].gamma_roll = model[j][i].gamma_roll = roll_coeffs_one[1];
        model[i][j].mu_roll    = model[j][i].mu_roll    = roll_coeffs_one[2];
      }

      // Define twisting model
      models[i][j].twist_model = models[j][i].twist_model = twist_model_one;
      if (twist_model_one != TWIST_NONE && twist_model_one != TWIST_MARSHALL) {
        model[i][j].k_twist = model[j][i].k_twist = twist_coeffs_one[0];
        model[i][j].gamma_twist = model[j][i].gamma_twist = twist_coeffs_one[1];
        model[i][j].mu_twist = model[j][i].mu_twist = twist_coeffs_one[2];
      }

      // Define extra options
      model[i][j].cutoff_type = model[j][i].cutoff_type = cutoff_one;
      model[i][j].limit_damping = model[j][i].limit_damping = ld_flag;

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

  // determine whether we need a granular neigh list, how large it needs to be

  use_history = normal_history || tangential_history ||
    roll_history || twist_history;

  // for JKR, will need fix/neigh/history to keep track of touch arrays

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (models[i][j].normal_model == JKR) use_history = 1;

  size_history = 3*tangential_history + 3*roll_history + twist_history;

  // determine location of tangential/roll/twist histories in array

  if (roll_history) {
    if (tangential_history) roll_history_index = 3;
    else roll_history_index = 0;
  }
  if (twist_history) {
    if (tangential_history) {
      if (roll_history) twist_history_index = 6;
      else twist_history_index = 3;
    } else {
      if (roll_history) twist_history_index = 3;
      else twist_history_index = 0;
    }
  }
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (model[i][j].tangential_model == TANGENTIAL_MINDLIN_RESCALE ||
          model[i][j].tangential_model == TANGENTIAL_MINDLIN_RESCALE_FORCE) {
        size_history += 1;
        roll_history_index += 1;
        twist_history_index += 1;
        nondefault_history_transfer = 1;
        history_transfer_factors = new int[size_history];
        for (int ii = 0; ii < size_history; ++ii)
          history_transfer_factors[ii] = -1;
        history_transfer_factors[3] = 1;
        break;
      }

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
    if ((models[i][i].normal_model != models[j][j].normal_model) ||
        (models[i][i].damping_model != models[j][j].damping_model) ||
        (models[i][i].tangential_model != models[j][j].tangential_model) ||
        (models[i][i].roll_model != models[j][j].roll_model) ||
        (models[i][i].twist_model != models[j][j].twist_model)) {
      error->all(FLERR,"Granular pair style functional forms are different, "
                 "cannot mix coefficients for types {} and {}. \n"
                 "This combination must be set explicitly via a "
                 "pair_coeff command",i,j);
    }

    // Mix normal coefficients
    if (models[i][j].normal_model == HERTZ || models[i][j].normal_model == HOOKE)
      models[i][j].k_norm = models[j][i][0].k_norm =
        mix_geom(models[i][i].k_norm, models[j][j].k_norm);
    else
      models[i][j].k_norm = models[j][i].k_norm =
        mix_stiffnessE(models[i][i].Emod, models[j][j].Emod,
                       models[i][i].poisson, models[j][j].poisson);

    models[i][j].gamma_norm = models[j][i].gamma_norm =
      mix_geom(models[i][i].gamma_norm, models[j][j].gamma_norm);
    if ((normal_model[i][j] == JKR) || (normal_model[i][j] == DMT))
      models[i][j].cohesion = models[j][i].cohesion =
        mix_geom(models[i][i].cohesion, models[j][j].cohesion);

    // Mix tangential coefficients
    models[i][j].k_tang = models[j][i].k_tang =
        mix_geom(models[i][i].k_tang, models[j][j].k_tang);
    models[i][j].gamma_tang = models[j][i].gamma_tang =
        mix_geom(models[i][i].gamma_tang, models[j][j].gamma_tang);
    models[i][j].mu_tang = models[j][i].mu_tang =
        mix_geom(models[i][i].mu_tang, models[j][j].mu_tang);

    // Mix rolling coefficients
    if (models.roll_model[i][j] != ROLL_NONE) {
      models[i][j].k_roll = models[j][i].k_roll =
        mix_geom(models[i][i].k_roll, models[j][j].k_roll);
      models[i][j].gamma_roll = models[j][i].gamma_roll =
          mix_geom(models[i][i].gamma_roll, models[j][j].gamma_roll);
      models[i][j].mu_roll = models[j][i].mu_roll =
          mix_geom(models[i][i].mu_roll, models[j][j].mu_roll);
    }

    // Mix twisting coefficients
    if (models[i][j].twist_model != TWIST_NONE && models[i][j].twist_model != TWIST_MARSHALL) {
      models[i][j].k_twist = models[j][i].k_twist =
        mix_geom(models[i][i].k_twist, models[j][j].k_twist);
      models[i][j].gamma_twist = models[j][i].gamma_twist =
        mix_geom(models[i][i].gamma_twist, models[j][j].gamma_twist);
      models[i][j].mu_twist = models[j][i].mu_twist =
        mix_geom(models[i][i].mu_twist, models[j][j].mu_twist);
    }
  }

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  double pulloff;
  if (models[i][j].cutoff_type < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
      pulloff = 0.0;
      if (normal_model[i][j] == JKR) {
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
  } else if (models[i][j].cutoff_type > 0) {
    cutoff = models[i][j].cutoff_type;
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

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
        fwrite(&models[i][j].normal_model,sizeof(int),1,fp);
        fwrite(&models[i][j].damping_model,sizeof(int),1,fp);
        fwrite(&models[i][j].tangential_model,sizeof(int),1,fp);
        fwrite(&models[i][j].roll_model,sizeof(int),1,fp);
        fwrite(&models[i][j].twist_model,sizeof(int),1,fp);
        fwrite(&models[i][j].limit_damping,sizeof(int),1,fp);
        fwrite(&models[i][j].Emod,sizeof(double),1,fp);
        fwrite(&models[i][j].poisson,sizeof(double),1,fp);
        fwrite(&models[i][j].k_norm,sizeof(double),1,fp);
        fwrite(&models[i][j].gamma_norm,sizeof(double),1,fp);
        fwrite(&models[i][j].cohesion,sizeof(double),1,fp);
        fwrite(&models[i][j].k_tang,sizeof(double),1,fp);
        fwrite(&models[i][j].gamma_tang,sizeof(double),1,fp);
        fwrite(&models[i][j].mu_tang,sizeof(double),1,fp);
        fwrite(&models[i][j].k_roll,sizeof(double),1,fp);
        fwrite(&models[i][j].gamma_roll,sizeof(double),1,fp);
        fwrite(&models[i][j].mu_roll,sizeof(double),1,fp);
        fwrite(&models[i][j].k_twist,sizeof(double),1,fp);
        fwrite(&models[i][j].gamma_twist,sizeof(double),1,fp);
        fwrite(&models[i][j].mu_twist,sizeof(double),1,fp);
        fwrite(&models[i][j].cutoff_type,sizeof(double),1,fp);
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
        if (me == 0) {
          utils::sfread(FLERR,&models[i][j].normal_model,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].damping_model,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].tangential_model,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].roll_model,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].twist_model,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].limit_damping,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].k_norm,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].gamma_norm,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].cohesion,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].k_tang,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].gamma_tang,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].mu_tang,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].k_roll,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].gamma_roll,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].mu_roll,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].k_twist,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].gamma_twist,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].mu_twist,sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&models[i][j].cutoff_type,sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&models[i][j].normal_model,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].damping_model,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].tangential_model,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].roll_model,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].twist_model,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].limit_damping,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].k_norm,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].gamma_norm,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].cohesion,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].k_tang,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].gamma_tang,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].mu_tang,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].k_roll,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].gamma_roll,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].mu_roll,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].k_twist,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].gamma_twist,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].mu_twist,1,MPI_INT,0,world);
        MPI_Bcast(&models[i][j].cutoff_type,1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
                            double rsq, double /* factor_coul */,
                            double /* factor_lj */, double &fforce)
{


  fforce = Fntot*rinv;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = fr1;
  svector[5] = fr2;
  svector[6] = fr3;
  svector[7] = fr;
  svector[8] = magtortwist;
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
   mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessE(double Eii, double Ejj,
                                    double poisii, double poisjj)
{
  return 1/((1-poisii*poisii)/Eii+(1-poisjj*poisjj)/Ejj);
}

/* ----------------------------------------------------------------------
   mixing of shear modulus (G)
------------------------------------------------------------------------ */

double PairGranular::mix_stiffnessG(double Eii, double Ejj,
                                    double poisii, double poisjj)
{
  return 1/((2*(2-poisii)*(1+poisii)/Eii) + (2*(2-poisjj)*(1+poisjj)/Ejj));
}

/* ----------------------------------------------------------------------
   mixing of everything else
------------------------------------------------------------------------- */

double PairGranular::mix_geom(double valii, double valjj)
{
  return sqrt(valii*valjj);
}

/* ----------------------------------------------------------------------
   transfer history during fix/neigh/history exchange
   only needed if any history entries i-j are not just negative of j-i entries
------------------------------------------------------------------------- */

void PairGranular::transfer_history(double* source, double* target)
{
  for (int i = 0; i < size_history; i++)
    target[i] = history_transfer_factors[i]*source[i];
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
    if (models[itype][itype].normal_model == JKR) {
      cut += models[itype][itype].pulloff_distance(cut, cut);
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
        if (models[i][j].normal_model == JKR) {
          temp = pulloff_distance(r1, r2);
          if (temp > cut) cut = temp;
        }
      }
    }
  }

  cut += r1 + r2;

  return cut;
}
