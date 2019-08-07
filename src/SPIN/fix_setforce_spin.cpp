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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_setforce_spin.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixSetForceSpin::FixSetForceSpin(LAMMPS *lmp, int narg, char **arg) :
  FixSetForce(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixSetForceSpin::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **fm = atom->fm;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += fm[i][0];
        foriginal[1] += fm[i][1];
        foriginal[2] += fm[i][2];
        if (xstyle) fm[i][0] = xvalue;
        if (ystyle) fm[i][1] = yvalue;
        if (zstyle) fm[i][2] = zvalue;
      }

  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += fm[i][0];
        foriginal[1] += fm[i][1];
        foriginal[2] += fm[i][2];
        if (xstyle == ATOM) fm[i][0] = sforce[i][0];
        else if (xstyle) fm[i][0] = xvalue;
        if (ystyle == ATOM) fm[i][1] = sforce[i][1];
        else if (ystyle) fm[i][1] = yvalue;
        if (zstyle == ATOM) fm[i][2] = sforce[i][2];
        else if (zstyle) fm[i][2] = zvalue;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForceSpin::single_setforce_spin(int i, double fmi[3])
{
  double **x = atom->x;
  int *mask = atom->mask;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;
 
  // constant force

  if (varflag == CONSTANT) {
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) return;
      foriginal[0] += fmi[0];
      foriginal[1] += fmi[1];
      foriginal[2] += fmi[2];
      if (xstyle) fmi[0] = xvalue;
      if (ystyle) fmi[1] = yvalue;
      if (zstyle) fmi[2] = zvalue;
    }

  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) return;
      foriginal[0] += fmi[0];
      foriginal[1] += fmi[1];
      foriginal[2] += fmi[2];
      if (xstyle == ATOM) fmi[0] = sforce[i][0];
      else if (xstyle) fmi[0] = xvalue;
      if (ystyle == ATOM) fmi[1] = sforce[i][1];
      else if (ystyle) fmi[1] = yvalue;
      if (zstyle == ATOM) fmi[2] = sforce[i][2];
      else if (zstyle) fmi[2] = zvalue;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForceSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // set force to desired value on requested level, 0.0 on other levels

  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = NULL;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }

    double **x = atom->x;
    double **fm = atom->fm;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if (xstyle) fm[i][0] = 0.0;
        if (ystyle) fm[i][1] = 0.0;
        if (zstyle) fm[i][2] = 0.0;
      }
  }
}
