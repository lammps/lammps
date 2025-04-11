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
   Charges by type and after option: Agilio Padua (ENS de Lyon & CNRS)
------------------------------------------------------------------------- */

#include "fix_adapt_fep.h"

#include "atom.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PAIR, KSPACE, ATOM};
enum{DIAMETER, CHARGE};

/* ---------------------------------------------------------------------- */

FixAdaptFEP::FixAdaptFEP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) utils::missing_cmd_args(FLERR,"fix adapt/fep", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 0) error->all(FLERR,"Illegal fix adapt/fep every value {}", nevery);

  dynamic_group_allow = 1;
  create_attribute = 1;

  // count # of adaptations

  nadapt = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) utils::missing_cmd_args(FLERR,"fix adapt/fep pair", error);
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix adapt/fep kspace", error);
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"fix adapt/fep atom", error);
      nadapt++;
      iarg += 4;
    } else break;
  }

  if (nadapt == 0) error->all(FLERR,"Nothing to adapt in fix adapt/fep command");
  adapt = new Adapt[nadapt];

  // parse keywords

  nadapt = 0;
  chgflag = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt/fep command");
      adapt[nadapt].which = PAIR;
      adapt[nadapt].pstyle = utils::strdup(arg[iarg+1]);
      adapt[nadapt].pparam = utils::strdup(arg[iarg+2]);
      utils::bounds_typelabel(FLERR, arg[iarg+3], 1, atom->ntypes,
                              adapt[nadapt].ilo, adapt[nadapt].ihi, lmp, Atom::ATOM);
      utils::bounds_typelabel(FLERR, arg[iarg+4], 1, atom->ntypes,
                              adapt[nadapt].jlo, adapt[nadapt].jhi, lmp, Atom::ATOM);

      // switch i,j if i > j, if wildcards were not used

      if ( (adapt[nadapt].ilo == adapt[nadapt].ihi) &&
           (adapt[nadapt].jlo == adapt[nadapt].jhi) &&
           (adapt[nadapt].ilo > adapt[nadapt].jlo) ) {
        adapt[nadapt].jlo = adapt[nadapt].ihi;
        adapt[nadapt].ilo = adapt[nadapt].jhi;
        adapt[nadapt].ihi = adapt[nadapt].ilo;
        adapt[nadapt].jhi = adapt[nadapt].jlo;
      }

      if (utils::strmatch(arg[iarg+5],"^v_")) {
        adapt[nadapt].var = utils::strdup(arg[iarg+5]+2);
      } else error->all(FLERR,"Illegal fix adapt/fep command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/fep command");
      adapt[nadapt].which = KSPACE;
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        adapt[nadapt].var = utils::strdup(arg[iarg+1]+2);
      } else error->all(FLERR,"Illegal fix adapt/fep command");
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt/fep command");
      adapt[nadapt].which = ATOM;
      if (strcmp(arg[iarg+1],"diameter") == 0) {
        adapt[nadapt].aparam = DIAMETER;
        diam_flag = 1;
      } else if (strcmp(arg[iarg+1],"charge") == 0) {
        adapt[nadapt].aparam = CHARGE;
        chgflag = 1;
      } else error->all(FLERR,"Illegal fix adapt/fep command");
      utils::bounds_typelabel(FLERR, arg[iarg+2], 1, atom->ntypes,
                              adapt[nadapt].ilo, adapt[nadapt].ihi, lmp, Atom::ATOM);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        adapt[nadapt].var = utils::strdup(arg[iarg+3]+2);
      } else error->all(FLERR,"Illegal fix adapt/fep command");
      nadapt++;
      iarg += 4;
    } else break;
  }

  // optional keywords

  resetflag = 0;
  scaleflag = 0;
  afterflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"reset") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix adapt/fep reset", error);
      resetflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix adapt/fep scale", error);
      scaleflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"after") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt/fep command");
      afterflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adapt/fep command");
  }

  // allocate pair style arrays

  int n = atom->ntypes;
  for (int m = 0; m < nadapt; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt:array_orig");

  id_fix_diam = id_fix_chg = nullptr;
}

/* ---------------------------------------------------------------------- */

FixAdaptFEP::~FixAdaptFEP()
{
  for (int m = 0; m < nadapt; m++) {
    delete[] adapt[m].var;
    if (adapt[m].which == PAIR) {
      delete[] adapt[m].pstyle;
      delete[] adapt[m].pparam;
      memory->destroy(adapt[m].array_orig);
    }
  }
  delete[] adapt;

  // check nfix in case all fixes have already been deleted

  if (id_fix_diam && modify->nfix) modify->delete_fix(id_fix_diam);
  if (id_fix_chg && modify->nfix) modify->delete_fix(id_fix_chg);
  delete[] id_fix_diam;
  delete[] id_fix_chg;
}

/* ---------------------------------------------------------------------- */

int FixAdaptFEP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_RUN;
  mask |= PRE_FORCE_RESPA;
  return mask;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */

void FixAdaptFEP::post_constructor()
{
  if (!resetflag) return;
  if (!diam_flag && !chgflag) return;

  // new id = fix-ID + FIX_STORE_ATTRIBUTE
  // new fix group = group for this fix

  id_fix_diam = nullptr;
  id_fix_chg = nullptr;

  if (diam_flag && atom->radius_flag) {
    id_fix_diam = utils::strdup(id + std::string("_FIX_STORE_DIAM"));
    fix_diam = dynamic_cast<FixStoreAtom *>(
      modify->add_fix(fmt::format("{} {} STORE/ATOM 1 0 0 1", id_fix_diam,group->names[igroup])));
    if (fix_diam->restart_reset) fix_diam->restart_reset = 0;
    else {
      double *vec = fix_diam->vstore;
      double *radius = atom->radius;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) vec[i] = radius[i];
        else vec[i] = 0.0;
      }
    }
  }

  if (chgflag && atom->q_flag) {
    id_fix_chg = utils::strdup(id + std::string("_FIX_STORE_CHG"));
    fix_chg = dynamic_cast<FixStoreAtom *>(
      modify->add_fix(fmt::format("{} {} STORE/ATOM 1 0 0 1",id_fix_chg,group->names[igroup])));
    if (fix_chg->restart_reset) fix_chg->restart_reset = 0;
    else {
      double *vec = fix_chg->vstore;
      double *q = atom->q;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) vec[i] = q[i];
        else vec[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::init()
{
  int i,j;

  // allow a dynamic group only if ATOM attribute not used

  if (group->dynamic[igroup])
    for (int i = 0; i < nadapt; i++)
      if (adapt[i].which == ATOM)
        error->all(FLERR,"Cannot use dynamic group with fix adapt/fep atom");

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    ad->ivar = input->variable->find(ad->var);
    if (ad->ivar < 0)
      error->all(FLERR,"Variable name {} for fix adapt/fep does not exist", ad->var);
    if (!input->variable->equalstyle(ad->ivar))
      error->all(FLERR,"Variable {} for fix adapt/fep is invalid style", ad->var);

    if (ad->which == PAIR) {
      anypair = 1;
      Pair *pair = nullptr;

      if (lmp->suffix_enable)
        pair = force->pair_match(std::string(ad->pstyle)+"/"+lmp->suffix,1);

      if (pair == nullptr) pair = force->pair_match(ad->pstyle,1);
      if (pair == nullptr)
        error->all(FLERR, "Fix adapt/fep pair style does not exist");
      void *ptr = pair->extract(ad->pparam,ad->pdim);
      if (ptr == nullptr)
        error->all(FLERR,"Fix adapt/fep pair style param not supported");

      if (ad->pdim != 2)
        error->all(FLERR,"Pair style parameter {} is not compatible with fix adapt/fep", ad->pparam);

      if (ad->pdim == 2) ad->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (ad->pdim == 2 && (strcmp(force->pair_style,"hybrid") == 0 ||
                            strcmp(force->pair_style,"hybrid/overlay") == 0)) {
        auto pair = dynamic_cast<PairHybrid *>(force->pair);
        for (i = ad->ilo; i <= ad->ihi; i++)
          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            if (!pair->check_ijtype(i,j,ad->pstyle))
              error->all(FLERR,"Fix adapt/fep type pair range is not valid for "
                         "pair hybrid sub-style");
      }

    } else if (ad->which == KSPACE) {
      if (force->kspace == nullptr)
        error->all(FLERR,"Fix adapt/fep kspace style does not exist");
      kspace_scale = (double *) force->kspace->extract("scale");

    } else if (ad->which == ATOM) {
      if (ad->aparam == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,"Fix adapt/fep requires atom attribute diameter");
      }
      if (ad->aparam == CHARGE) {
        if (!atom->q_flag)
          error->all(FLERR,"Fix adapt/fep requires atom attribute charge");
      }
    }
  }

  // make copy of original pair array values

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR && ad->pdim == 2) {
      for (i = ad->ilo; i <= ad->ihi; i++)
        for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
          ad->array_orig[i][j] = ad->array[i][j];
    }
  }

  // fixes that store initial per-atom values

  if (id_fix_diam) {
    fix_diam = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(id_fix_diam));
    if (!fix_diam) error->all(FLERR,"Could not find fix adapt/fep storage fix ID {}", id_fix_diam);
  }
  if (id_fix_chg) {
    fix_chg = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(id_fix_chg));
    if (!fix_chg) error->all(FLERR,"Could not find fix adapt/fep storage fix ID {}", id_fix_chg);
  }

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::setup_pre_force(int /*vflag*/)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::pre_force(int /*vflag*/)
{
  if (nevery == 0) return;

  if (afterflag) {         // update at n+1 (better with fix ave/time)
    if (nevery == 1 || update->ntimestep == 0)
      change_settings();
    else if (update->ntimestep > 1 && !((update->ntimestep - 1) % nevery))
      change_settings();
  } else {                      // original version: update at n
    if (update->ntimestep % nevery)
      return;
    change_settings();
  }
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::pre_force_respa(int vflag, int ilevel, int)
{
  if (ilevel < nlevels_respa-1) return;
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptFEP::post_run()
{
  if (resetflag) restore_settings();
}

/* ----------------------------------------------------------------------
   change pair,kspace,atom parameters based on variable evaluation
------------------------------------------------------------------------- */

void FixAdaptFEP::change_settings()
{
  int i,j;

  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    double value = input->variable->compute_equal(ad->ivar);

    // set global scalar or type pair array values

    if (ad->which == PAIR) {
      if (ad->pdim == 0) {
        if (scaleflag) *ad->scalar = value * ad->scalar_orig;
        else *ad->scalar = value;
      } else if (ad->pdim == 2) {
        if (scaleflag)
          for (i = ad->ilo; i <= ad->ihi; i++)
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
              ad->array[i][j] = value*ad->array_orig[i][j];
        else
          for (i = ad->ilo; i <= ad->ihi; i++)
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
              ad->array[i][j] = value;
      }

    // set kspace scale factor

    } else if (ad->which == KSPACE) {
      *kspace_scale = value;

    // set per atom values, also make changes for ghost atoms

    } else if (ad->which == ATOM) {

      if (scaleflag)
        error->all(FLERR, "Keyword 'scale yes' is not supported with fix adapt/fep for 'atom'");
      // reset radius from diameter
      // also scale rmass to new value

      if (ad->aparam == DIAMETER) {
        int mflag = 0;
        if (atom->rmass_flag) mflag = 1;
        double density;

        int *atype = atom->type;
        double *radius = atom->radius;
        double *rmass = atom->rmass;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        int nall = nlocal + atom->nghost;

        if (mflag == 0) {
          for (i = 0; i < nall; i++)
            if (atype[i] >= ad->ilo && atype[i] <= ad->ihi)
              if (mask[i] & groupbit)
                radius[i] = 0.5*value;
        } else {
          for (i = 0; i < nall; i++)
            if (atype[i] >= ad->ilo && atype[i] <= ad->ihi)
              if (mask[i] & groupbit) {
                density = rmass[i] / (4.0*MY_PI/3.0 *
                                      radius[i]*radius[i]*radius[i]);
                radius[i] = 0.5*value;
                rmass[i] = 4.0*MY_PI/3.0 *
                  radius[i]*radius[i]*radius[i] * density;
              }
        }
      } else if (ad->aparam == CHARGE) {
        int *atype = atom->type;
        double *q = atom->q;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        int nall = nlocal + atom->nghost;

        for (i = 0; i < nall; i++)
          if (atype[i] >= ad->ilo && atype[i] <= ad->ihi)
            if (mask[i] & groupbit) q[i] = value;
      }
    }
  }

  modify->addstep_compute(update->ntimestep + nevery);

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (anypair) force->pair->reinit();

  // reset KSpace charges if charges have changed

  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ----------------------------------------------------------------------
   restore pair,kspace.atom parameters to original values
------------------------------------------------------------------------- */

void FixAdaptFEP::restore_settings()
{
  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR) {
      if (ad->pdim == 0) *ad->scalar = ad->scalar_orig;
      else if (ad->pdim == 2) {
        for (int i = ad->ilo; i <= ad->ihi; i++)
          for (int j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            ad->array[i][j] = ad->array_orig[i][j];
      }

    } else if (ad->which == KSPACE) {
      *kspace_scale = 1.0;

    } else if (ad->which == ATOM) {
      if (diam_flag) {
        double density;

        double *vec = fix_diam->vstore;
        double *radius = atom->radius;
        double *rmass = atom->rmass;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            density = rmass[i] / (4.0*MY_PI/3.0 *
                                  radius[i]*radius[i]*radius[i]);
            radius[i] = vec[i];
            rmass[i] = 4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i] * density;
          }
      }
      if (chgflag) {
        double *vec = fix_chg->vstore;
        double *q = atom->q;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) q[i] = vec[i];
      }
    }
  }

  if (anypair) force->pair->reinit();
  if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void FixAdaptFEP::set_arrays(int i)
{
  if (fix_diam) fix_diam->vstore[i] = atom->radius[i];
  if (fix_chg) fix_chg->vstore[i] = atom->q[i];
}
