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

#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "force.h"
#include "atom.h"
#include "comm.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "bond.h"
#include "bond_hybrid.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "group.h"
#include "memory.h"
#include "error.h"

#define BondInclude
#define AngleInclude
#define DihedralInclude
#define ImproperInclude
#define PairInclude
#define KSpaceInclude
#include "style.h"
#undef BondInclude
#undef AngleInclude
#undef DihedralInclude
#undef ImproperInclude
#undef PairInclude
#undef KSpaceInclude

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Force::Force(LAMMPS *lmp) : Pointers(lmp)
{
  newton = newton_pair = newton_bond = 1;

  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_dihedral = 0;
  special_extra = 0;

  dielectric = 1.0;

  pair = NULL;
  bond = NULL;
  angle = NULL;
  dihedral = NULL;
  improper = NULL;
  kspace = NULL;

  char *str = (char *) "none";
  int n = strlen(str) + 1;
  pair_style = new char[n];
  strcpy(pair_style,str);
  bond_style = new char[n];
  strcpy(bond_style,str);
  angle_style = new char[n];
  strcpy(angle_style,str);
  dihedral_style = new char[n];
  strcpy(dihedral_style,str);
  improper_style = new char[n];
  strcpy(improper_style,str);
  kspace_style = new char[n];
  strcpy(kspace_style,str);
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
  delete [] pair_style;
  delete [] bond_style;
  delete [] angle_style;
  delete [] dihedral_style;
  delete [] improper_style;
  delete [] kspace_style;

  if (pair) delete pair;
  if (bond) delete bond;
  if (angle) delete angle;
  if (dihedral) delete dihedral;
  if (improper) delete improper;
  if (kspace) delete kspace;
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
  qqrd2e = qqr2e/dielectric;

  comm->maxforward_pair = comm->maxreverse_pair = 0;

  if (kspace) kspace->init();         // kspace must come before pair
  if (pair) pair->init();             // so g_ewald is defined
  if (bond) bond->init();
  if (angle) angle->init();
  if (dihedral) dihedral->init();
  if (improper) improper->init();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_pair(const char *style)
{
  delete [] pair_style;
  if (pair) delete pair;

  pair = new_pair(style);
  int n = strlen(style) + 1;
  pair_style = new char[n];
  strcpy(pair_style,style);
}

/* ----------------------------------------------------------------------
   generate a pair class
------------------------------------------------------------------------- */

Pair *Force::new_pair(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define PairClass
#define PairStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style.h"
#undef PairClass

  else error->all("Invalid pair style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to current pair class or hybrid sub-class
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   return NULL if no match
   return NULL if not exact and multiple hybrid sub-styles match
------------------------------------------------------------------------- */

Pair *Force::pair_match(const char *word, int exact)
{
  int iwhich,count;

  if (exact && strcmp(pair_style,word) == 0) return pair;
  else if (!exact && strstr(pair_style,word)) return pair;

  else if (strcmp(pair_style,"hybrid") == 0) {
    PairHybrid *hybrid = (PairHybrid *) pair;
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++) {
      if (exact && strcmp(hybrid->keywords[i],word) == 0)
	return hybrid->styles[i];
      if (!exact && strstr(hybrid->keywords[i],word)) {
	iwhich = i;
	count++;
      }
    }
    if (!exact && count == 1) return hybrid->styles[iwhich];

  } else if (strcmp(pair_style,"hybrid/overlay") == 0) {
    PairHybridOverlay *hybrid = (PairHybridOverlay *) pair;
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++) {
      if (exact && strcmp(hybrid->keywords[i],word) == 0)
	return hybrid->styles[i];
      else if (!exact && strstr(hybrid->keywords[i],word)) {
	iwhich = i;
	count++;
      }
    }
    if (!exact && count == 1) return hybrid->styles[iwhich];
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   create a bond style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_bond(const char *style)
{
  delete [] bond_style;
  if (bond) delete bond;

  bond = new_bond(style);
  int n = strlen(style) + 1;
  bond_style = new char[n];
  strcpy(bond_style,style);
}

/* ----------------------------------------------------------------------
   generate a bond class
------------------------------------------------------------------------- */

Bond *Force::new_bond(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define BondClass
#define BondStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style.h"
#undef BondClass

  else error->all("Invalid bond style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to current bond class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Bond *Force::bond_match(const char *style)
{
  if (strcmp(bond_style,style) == 0) return bond;
  else if (strcmp(bond_style,"hybrid") == 0) {
    BondHybrid *hybrid = (BondHybrid *) bond;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strcmp(hybrid->keywords[i],style) == 0) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create an angle style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_angle(const char *style)
{
  delete [] angle_style;
  if (angle) delete angle;

  angle = new_angle(style);
  int n = strlen(style) + 1;
  angle_style = new char[n];
  strcpy(angle_style,style);
}

/* ----------------------------------------------------------------------
   generate an angle class
------------------------------------------------------------------------- */

Angle *Force::new_angle(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define AngleClass
#define AngleStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style.h"
#undef AngleClass

  else error->all("Invalid angle style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create a dihedral style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_dihedral(const char *style)
{
  delete [] dihedral_style;
  if (dihedral) delete dihedral;

  dihedral = new_dihedral(style);
  int n = strlen(style) + 1;
  dihedral_style = new char[n];
  strcpy(dihedral_style,style);
}

/* ----------------------------------------------------------------------
   generate a dihedral class
------------------------------------------------------------------------- */

Dihedral *Force::new_dihedral(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define DihedralClass
#define DihedralStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style.h"
#undef DihedralClass

  else error->all("Invalid dihedral style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create an improper style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_improper(const char *style)
{
  delete [] improper_style;
  if (improper) delete improper;

  improper = new_improper(style);
  int n = strlen(style) + 1;
  improper_style = new char[n];
  strcpy(improper_style,style);
}

/* ----------------------------------------------------------------------
   generate a improper class
------------------------------------------------------------------------- */

Improper *Force::new_improper(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define ImproperClass
#define ImproperStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style.h"
#undef ImproperClass

  else error->all("Invalid improper style");
  return NULL;
}

/* ----------------------------------------------------------------------
   new kspace style 
------------------------------------------------------------------------- */

void Force::create_kspace(int narg, char **arg)
{
  delete [] kspace_style;
  if (kspace) delete kspace;

  if (strcmp(arg[0],"none") == 0) kspace = NULL;

#define KSpaceClass
#define KSpaceStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) kspace = new Class(lmp,narg-1,&arg[1]);
#include "style.h"
#undef KSpaceClass

  else error->all("Invalid kspace style");

  int n = strlen(arg[0]) + 1;
  kspace_style = new char[n];
  strcpy(kspace_style,arg[0]);
}

/* ----------------------------------------------------------------------
   set special bond values 
------------------------------------------------------------------------- */

void Force::set_special(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal special_bonds command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"amber") == 0) {
      if (iarg+1 > narg) error->all("Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.5;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 5.0/6.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"charmm") == 0) {
      if (iarg+1 > narg) error->all("Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 0.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"fene") == 0) {
      if (iarg+1 > narg) error->all("Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 1.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 1.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"lj/coul") == 0) {
      if (iarg+4 > narg) error->all("Illegal special_bonds command");
      special_lj[1] = special_coul[1] = atof(arg[iarg+1]);
      special_lj[2] = special_coul[2] = atof(arg[iarg+2]);
      special_lj[3] = special_coul[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"lj") == 0) {
      if (iarg+4 > narg) error->all("Illegal special_bonds command");
      special_lj[1] = atof(arg[iarg+1]);
      special_lj[2] = atof(arg[iarg+2]);
      special_lj[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"coul") == 0) {
      if (iarg+4 > narg) error->all("Illegal special_bonds command");
      special_coul[1] = atof(arg[iarg+1]);
      special_coul[2] = atof(arg[iarg+2]);
      special_coul[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all("Illegal special_bonds command");
      if (strcmp(arg[iarg+1],"no") == 0) special_dihedral = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) special_dihedral = 1;
      else error->all("Illegal special_bonds command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all("Illegal special_bonds command");
      special_extra = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal special_bonds command");
  }

  for (int i = 1; i <= 3; i++)
    if (special_lj[i] < 0.0 || special_lj[i] > 1.0 ||
	special_coul[i] < 0.0 || special_coul[i] > 1.0)
      error->all("Illegal special_bonds command");

  if (special_extra < 0) error->all("Illegal special_bonds command");
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to nmax,
     (3) i* = i to nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Force::bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),1);
    nhi = MIN(atoi(str),nmax);
  } else if (strlen(str) == 1) {
    nlo = 1;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 1;
    nhi = MIN(atoi(ptr+1),nmax);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),1);
    nhi = nmax;
  } else {
    nlo = MAX(atoi(str),1);
    nhi = MIN(atoi(ptr+1),nmax);
  }
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by force fields to check validity of their arguments
------------------------------------------------------------------------- */

double Force::numeric(char *str)
{
  int n = strlen(str);
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    error->all("Expecting floating point argument in input script");
  }

  return atof(str);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by force fields to check validity of their arguments
------------------------------------------------------------------------- */

int Force::inumeric(char *str)
{
  int n = strlen(str);
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all("Expecting integer argument in input script");
  }

  return atoi(str);
}

/* ----------------------------------------------------------------------
   memory usage of force classes
------------------------------------------------------------------------- */

double Force::memory_usage()
{
  double bytes = 0.0;
  if (pair) bytes += pair->memory_usage();
  if (bond) bytes += bond->memory_usage();
  if (angle) bytes += angle->memory_usage();
  if (dihedral) bytes += dihedral->memory_usage();
  if (improper) bytes += improper->memory_usage();
  if (kspace) bytes += kspace->memory_usage();
  return bytes;
}
