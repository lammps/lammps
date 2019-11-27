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

#include "atom_vec_hybrid.h"
#include <cstring>
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybrid::AtomVecHybrid(LAMMPS *lmp) : AtomVec(lmp) {}

/* ---------------------------------------------------------------------- */

AtomVecHybrid::~AtomVecHybrid()
{
  for (int k = 0; k < nstyles; k++) delete styles[k];
  delete [] styles;
  for (int k = 0; k < nstyles; k++) delete [] keywords[k];
  delete [] keywords;

  // these strings will be concatenated from sub-style strings
  // fields_data_atom must start with fields common to all styles

  fields_grow = fields_copy = fields_comm = fields_comm_vel = NULL;
  fields_reverse = fields_border = fields_border_vel = NULL;
  fields_exchange = fields_restart = fields_create = NULL;
  fields_data_atom = (char *) "id type x";
  fields_data_vel = NULL;
}

/* ----------------------------------------------------------------------
   process sub-style args
------------------------------------------------------------------------- */

void AtomVecHybrid::process_args(int narg, char **arg)
{
  // build list of all known atom styles

  build_styles();

  // allocate list of sub-styles as big as possibly needed if no extra args

  styles = new AtomVec*[narg];
  keywords = new char*[narg];

  // allocate each sub-style
  // call process_args() with set of args that are not atom style names
  // use known_style() to determine which args these are

  int i,jarg,dummy;

  int iarg = 0;
  nstyles = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"hybrid") == 0)
      error->all(FLERR,"Atom style hybrid cannot have hybrid as an argument");
    for (i = 0; i < nstyles; i++)
      if (strcmp(arg[iarg],keywords[i]) == 0)
        error->all(FLERR,"Atom style hybrid cannot use same atom style twice");
    styles[nstyles] = atom->new_avec(arg[iarg],1,dummy);
    keywords[nstyles] = new char[strlen(arg[iarg])+1];
    strcpy(keywords[nstyles],arg[iarg]);
    jarg = iarg + 1;
    while (jarg < narg && !known_style(arg[jarg])) jarg++;
    styles[nstyles]->process_args(jarg-iarg-1,&arg[iarg+1]);
    iarg = jarg;
    nstyles++;
  }

  // free allstyles created by build_styles()

  for (int i = 0; i < nallstyles; i++) delete [] allstyles[i];
  delete [] allstyles;

  // concatenate field strings from all sub-styles

  concatenate_fields();

  // parent AtomVec will now operate on concatenated fields

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::concatenate_fields()
{
  for (int k = 0; k < nstyles; k++) {
    concatenate(fields_grow,styles[k]->fields_grow);
    concatenate(fields_copy,styles[k]->fields_copy);
    concatenate(fields_comm,styles[k]->fields_comm);
    concatenate(fields_comm_vel,styles[k]->fields_comm_vel);
    concatenate(fields_reverse,styles[k]->fields_reverse);
    concatenate(fields_border,styles[k]->fields_border);
    concatenate(fields_border_vel,styles[k]->fields_border_vel);
    concatenate(fields_exchange,styles[k]->fields_exchange);
    concatenate(fields_restart,styles[k]->fields_restart);
    concatenate(fields_create,styles[k]->fields_create);
    concatenate(fields_data_atom,styles[k]->fields_data_atom);
    concatenate(fields_data_vel,styles[k]->fields_data_vel);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::concatenate(char *&root, char *add)
{
  /*
  char **rootwords,**addwords;
  int nroot = parse(root,rootwords);
  int nadd = parse(add,addwords);

  for (int iadd = 0; iadd < nadd; iadd++) {
    if (check(addwords[iadd],nroot,rootwords)) continue;
    addone(addwords[iadd],nroot,rootwords);
  }
  */
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::init()
{
  AtomVec::init();
  for (int k = 0; k < nstyles; k++) styles[k]->init();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybrid::force_clear(int n, size_t nbytes)
{
  for (int k = 0; k < nstyles; k++)
    if (styles[k]->forceclearflag) styles[k]->force_clear(n,nbytes);
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   returned value encodes which sub-style and index returned by sub-style
   return -1 if name is unknown to any sub-styles
------------------------------------------------------------------------- */

int AtomVecHybrid::property_atom(char *name)
{
  for (int k = 0; k < nstyles; k++) {
    int index = styles[k]->property_atom(name);
    if (index >= 0) return index*nstyles + k;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecHybrid::pack_property_atom(int multiindex, double *buf,
                                       int nvalues, int groupbit)
{
  int k = multiindex % nstyles;
  int index = multiindex/nstyles;
  styles[k]->pack_property_atom(index,buf,nvalues,groupbit);
}

/* ----------------------------------------------------------------------
   allstyles = list of all atom styles in this LAMMPS executable
------------------------------------------------------------------------- */

void AtomVecHybrid::build_styles()
{
  nallstyles = 0;
#define ATOM_CLASS
#define AtomStyle(key,Class) nallstyles++;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS

  allstyles = new char*[nallstyles];

  int n;
  nallstyles = 0;
#define ATOM_CLASS
#define AtomStyle(key,Class)                \
  n = strlen(#key) + 1;                     \
  allstyles[nallstyles] = new char[n];      \
  strcpy(allstyles[nallstyles],#key);       \
  nallstyles++;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
}

/* ----------------------------------------------------------------------
   allstyles = list of all known atom styles
------------------------------------------------------------------------- */

int AtomVecHybrid::known_style(char *str)
{
  for (int i = 0; i < nallstyles; i++)
    if (strcmp(str,allstyles[i]) == 0) return 1;
  return 0;
}
