/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_adapt.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PAIR,KSPACE,ATOM};
enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixAdapt::FixAdapt(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix adapt command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");

  // count # of adaptations

  nadapt = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 3;
    } else break;
  }

  if (nadapt == 0) error->all(FLERR,"Illegal fix adapt command");
  adapt = new Adapt[nadapt];

  // parse keywords

  nadapt = 0;
  diamflag = 0;
  chgflag = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      adapt[nadapt].pstyle = new char[n];
      strcpy(adapt[nadapt].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      adapt[nadapt].pparam = new char[n];
      strcpy(adapt[nadapt].pparam,arg[iarg+2]);
      force->bounds(arg[iarg+3],atom->ntypes,
                    adapt[nadapt].ilo,adapt[nadapt].ihi);
      force->bounds(arg[iarg+4],atom->ntypes,
                    adapt[nadapt].jlo,adapt[nadapt].jhi);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        n = strlen(&arg[iarg+5][2]) + 1;
        adapt[nadapt].var = new char[n];
        strcpy(adapt[nadapt].var,&arg[iarg+5][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = KSPACE;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        adapt[nadapt].var = new char[n];
        strcpy(adapt[nadapt].var,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
      adapt[nadapt].which = ATOM;
      if (strcmp(arg[iarg+1],"diameter") == 0) {
        adapt[nadapt].aparam = DIAMETER;
        diamflag = 1;
      } else if (strcmp(arg[iarg+1],"charge") == 0) {
        adapt[nadapt].aparam = CHARGE; 
        chgflag = 1; 
      } else error->all(FLERR,"Illegal fix adapt command");
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        adapt[nadapt].var = new char[n];
        strcpy(adapt[nadapt].var,&arg[iarg+2][2]);
      } else error->all(FLERR,"Illegal fix adapt command");
      nadapt++;
      iarg += 3;
    } else break;
  }

  // optional keywords

  resetflag = 0;
  scaleflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"reset") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      if (strcmp(arg[iarg+1],"no") == 0) resetflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) resetflag = 1;
      else error->all(FLERR,"Illegal fix adapt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
      if (strcmp(arg[iarg+1],"no") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix adapt command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adapt command");
  }

  // allocate pair style arrays

  int n = atom->ntypes;
  for (int m = 0; m < nadapt; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt:array_orig");
}

/* ---------------------------------------------------------------------- */

FixAdapt::~FixAdapt()
{
  for (int m = 0; m < nadapt; m++) {
    delete [] adapt[m].var;
    if (adapt[m].which == PAIR) {
      delete [] adapt[m].pstyle;
      delete [] adapt[m].pparam;
      memory->destroy(adapt[m].array_orig);
    }
  }
  delete [] adapt;
}

/* ---------------------------------------------------------------------- */

int FixAdapt::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdapt::init()
{
  int i,j;

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    ad->ivar = input->variable->find(ad->var);
    if (ad->ivar < 0)
      error->all(FLERR,"Variable name for fix adapt does not exist");
    if (!input->variable->equalstyle(ad->ivar))
      error->all(FLERR,"Variable for fix adapt is invalid style");

    if (ad->which == PAIR) {
      anypair = 1;
      Pair *pair = NULL;

      if (lmp->suffix_enable) {
        char psuffix[128];
        strcpy(psuffix,ad->pstyle);
        strcat(psuffix,"/");
        strcat(psuffix,lmp->suffix);
        pair = force->pair_match(psuffix,1);
      }
      if (pair == NULL) pair = force->pair_match(ad->pstyle,1);
      if (pair == NULL) error->all(FLERR,"Fix adapt pair style does not exist");
      void *ptr = pair->extract(ad->pparam,ad->pdim);
      if (ptr == NULL) 
        error->all(FLERR,"Fix adapt pair style param not supported");

      ad->pdim = 2;
      if (ad->pdim == 0) ad->scalar = (double *) ptr;
      if (ad->pdim == 2) ad->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (ad->pdim == 2 && (strcmp(force->pair_style,"hybrid") == 0 ||
                            strcmp(force->pair_style,"hybrid/overlay") == 0)) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = ad->ilo; i <= ad->ihi; i++)
          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            if (!pair->check_ijtype(i,j,ad->pstyle))
              error->all(FLERR,"Fix adapt type pair range is not valid for "
                         "pair hybrid sub-style");
      }

    } else if (ad->which == KSPACE) {
      if (force->kspace == NULL)
        error->all(FLERR,"Fix adapt kspace style does not exist");
      kspace_scale = (double *) force->kspace->extract("scale");

    } else if (ad->which == ATOM) {
      if (ad->aparam == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,"Fix adapt requires atom attribute diameter");
      }
      if (ad->aparam == CHARGE) {
	if (!atom->q_flag)
	  error->all(FLERR,"Fix adapt requires atom attribute charge");
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
}

/* ---------------------------------------------------------------------- */

void FixAdapt::setup_pre_force(int vflag)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdapt::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdapt::post_run()
{
  if (resetflag) restore_settings();
}

/* ----------------------------------------------------------------------
   change pair,kspace,atom parameters based on variable evaluation
------------------------------------------------------------------------- */

void FixAdapt::change_settings()
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

    } else if (ad->which == ATOM) {

      // set radius from diameter
      // also scale rmass to new value

      if (ad->aparam == DIAMETER) {
        int mflag = 0;
        if (atom->rmass_flag) mflag = 1;
        double density;

        double *radius = atom->radius;
        double *rmass = atom->rmass;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        if (mflag == 0) {
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit)
              radius[i] = 0.5*value;
        } else {
          for (i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) {
              density = rmass[i] / (4.0*MY_PI/3.0 *
                                    radius[i]*radius[i]*radius[i]);
              radius[i] = 0.5*value;
              rmass[i] = 4.0*MY_PI/3.0 *
                radius[i]*radius[i]*radius[i] * density;
            }
        }
      } else if (ad->aparam == CHARGE) {
        double *q = atom->q; 
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	for (i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) q[i] = value; 
      }
    }
  }

  modify->addstep_compute(update->ntimestep + nevery);

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (anypair) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   restore pair,kspace,atom parameters to original values
------------------------------------------------------------------------- */

void FixAdapt::restore_settings()
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

    }
  }

  if (anypair) force->pair->reinit();
}
