// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author:  Samuel Cajahuaringa Macollunco (CCES-UNICAMP/BR)
   lambda scaling forces parameters works with the dynamical
   Clausius-Clapeyron Integration (dcci) method
------------------------------------------------------------------------- */

#include "fix_adapt_dcci.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PAIR, ATOM};

/* ---------------------------------------------------------------------- */

FixAdaptDCCI::FixAdaptDCCI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
nadapt(0), adapt(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix adapt/dcci command");
  lambda = utils::numeric(FLERR,arg[3],false,lmp);
  if (lambda < 0) error->all(FLERR,"Illegal fix adapt/dcci command");

  dynamic_group_allow = 1;
  create_attribute = 1;
  scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 2;
  // count # of adaptations
  nadapt = 0;

  int iarg = 4;  // atribute pair 
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix adapt/dcci command");
      nadapt++;
      iarg += 5;
    } else break;
  }

  if (nadapt == 0) error->all(FLERR,"Illegal fix adapt/dcci command");
  adapt = new Adapt[nadapt];

  // parse keywords

  nadapt = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix adapt/dcci command");
      adapt[nadapt].which = PAIR;
      adapt[nadapt].pstyle = utils::strdup(arg[iarg+1]);
      adapt[nadapt].pparam = utils::strdup(arg[iarg+2]);
      adapt[nadapt].pair = nullptr;
      utils::bounds(FLERR,arg[iarg+3],1,atom->ntypes,
                    adapt[nadapt].ilo,adapt[nadapt].ihi,error);
      utils::bounds(FLERR,arg[iarg+4],1,atom->ntypes,
                    adapt[nadapt].jlo,adapt[nadapt].jhi,error);

      nadapt++;
      iarg += 5;
    } else break;
  }

  // allocate pair style arrays

  int n = atom->ntypes;
  for (int m = 0; m < nadapt; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt/dcci:array_orig");

}

/* ---------------------------------------------------------------------- */

FixAdaptDCCI::~FixAdaptDCCI()
{
  for (int m = 0; m < nadapt; m++) {
    if (adapt[m].which == PAIR) {
      delete [] adapt[m].pstyle;
      delete [] adapt[m].pparam;
      memory->destroy(adapt[m].array_orig);
    }
  }
  delete [] adapt;
}

/* ---------------------------------------------------------------------- */

int FixAdaptDCCI::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_RUN;
  mask |= PRE_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::init()
{
  int i,j;

  // allow a dynamic group only if ATOM attribute not used

  if (group->dynamic[igroup])
    for (int i = 0; i < nadapt; i++)
      if (adapt[i].which == ATOM)
        error->all(FLERR,"Cannot use dynamic group with fix adapt/dcci atom");

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    if (ad->which == PAIR) {
      anypair = 1;
      ad->pair = nullptr;

      // if ad->pstyle has trailing sub-style annotation ":N",
      //   strip it for pstyle arg to pair_match() and set nsub = N
      // this should work for appended suffixes as well

      char *pstyle = utils::strdup(ad->pstyle);
      char *cptr;
      int nsub = 0;
      if ((cptr = strchr(pstyle,':'))) {
        *cptr = '\0';
        nsub = utils::inumeric(FLERR,cptr+1,false,lmp);
      }

      if (lmp->suffix_enable) {
        if (lmp->suffix)
          ad->pair = force->pair_match(fmt::format("{}/{}",pstyle,lmp->suffix),1,nsub);
        if ((ad->pair == nullptr) && lmp->suffix2)
          ad->pair = force->pair_match(fmt::format("{}/{}",pstyle,lmp->suffix2),1,nsub);
      }

      if (ad->pair == nullptr) ad->pair = force->pair_match(pstyle,1,nsub);
      if (ad->pair == nullptr) error->all(FLERR,"Fix adapt pair style {} not found", pstyle);

      void *ptr = ad->pair->extract(ad->pparam,ad->pdim);
      if (ptr == nullptr)
        error->all(FLERR,"Fix adapt/dcci pair style param not supported");

      // for pair styles only parameters that are 2-d arrays in atom types or
      // scalars are supported

      if (ad->pdim != 2 && ad->pdim != 0)
        error->all(FLERR,"Fix adapt/dcci pair style param is not compatible");

      if (ad->pdim == 2) ad->array = (double **) ptr;
      if (ad->pdim == 0) ad->scalar = (double *) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (utils::strmatch(force->pair_style,"^hybrid")) {
        auto pair = dynamic_cast<PairHybrid *>(force->pair);
        for (i = ad->ilo; i <= ad->ihi; i++)
          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            if (!pair->check_ijtype(i,j,pstyle))
              error->all(FLERR,"Fix adapt/dcci type pair range is not valid "
                         "for pair hybrid sub-style {}", pstyle);
      }

      delete [] pstyle;
    }
  }

  // make copy of original pair array values

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR && ad->pdim == 2) {
      for (i = ad->ilo; i <= ad->ihi; i++)
        for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
          ad->array_orig[i][j] = ad->array[i][j];
    } else if (ad->which == PAIR && ad->pdim == 0){
      ad->scalar_orig = *ad->scalar;
    }
  }

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::setup_pre_force(int /*vflag*/)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::pre_force(int /*vflag*/)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::pre_force_respa(int vflag, int ilevel, int)
{
  if (ilevel < nlevels_respa-1) return;
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::post_run()
{
  //if (resetflag) restore_settings();
}

/* ----------------------------------------------------------------------
   change pair,kspace parameters based on variable evaluation
------------------------------------------------------------------------- */

void FixAdaptDCCI::change_settings()
{
  int i,j;

  // variable evaluation may invoke computes so wrap with clear/add

  //modify->clearstep_compute();

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    double value = lambda; 

    // set global scalar or type pair array values

    if (ad->which == PAIR) {
      if (ad->pdim == 0) {
        *ad->scalar = value;
      } else if (ad->pdim == 2) {
          for (i = ad->ilo; i <= ad->ihi; i++)
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
              ad->array[i][j] = value;
      }
    }
  }

  if (anypair) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   restore pair,kspace,atom parameters to original values
------------------------------------------------------------------------- */

void FixAdaptDCCI::restore_settings()
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
    }
  }

  if (anypair) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   lambda scaling parameter to force calculations
------------------------------------------------------------------------- */

double FixAdaptDCCI::compute_scalar()
{
  return lambda;
}

