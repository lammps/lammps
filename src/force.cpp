/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
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
#include "temperature.h"
#include "pressure.h"
#include "group.h"
#include "memory.h"
#include "error.h"

#define BondInclude
#define AngleInclude
#define DihedralInclude
#define ImproperInclude
#define PairInclude
#define KSpaceInclude
#define TempInclude
#include "style.h"
#undef BondInclude
#undef AngleInclude
#undef DihedralInclude
#undef ImproperInclude
#undef PairInclude
#undef KSpaceInclude
#undef TempInclude

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 1

/* ---------------------------------------------------------------------- */

Force::Force()
{
  dimension = 3;
  newton = newton_pair = newton_bond = 1;

  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;

  dielectric = 1.0;

  pair = NULL;
  bond = NULL;
  angle = NULL;
  dihedral = NULL;
  improper = NULL;
  kspace = NULL;

  char *str = "none";
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

  // default temperature = group all and style full

  ntemp = maxtemp = 0;
  templist = NULL;

  char **arg = new char*[3];
  arg[0] = "default";
  arg[1] = "all";
  arg[2] = "full";
  add_temp(3,arg,0);
  delete [] arg;

  pressure = new Pressure;
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

  for (int i = 0; i < ntemp; i++) delete templist[i];
  memory->sfree(templist);
  delete pressure;
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

  for (int itemp = 0; itemp < ntemp; itemp++) templist[itemp]->init();
  pressure->init();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_pair(char *style)
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

Pair *Force::new_pair(char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define PairClass
#define PairStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class();
#include "style.h"
#undef PairClass

  else error->all("Invalid pair style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to current pair class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Pair *Force::pair_match(char *style)
{
  if (strcmp(pair_style,style) == 0) return pair;
  else if (strcmp(pair_style,"hybrid") == 0) {
    PairHybrid *hpair = (PairHybrid *) pair;
    for (int i = 0; i < hpair->nstyles; i++)
      if (strcmp(hpair->keywords[i],style) == 0) return hpair->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create a bond style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_bond(char *style)
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

Bond *Force::new_bond(char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define BondClass
#define BondStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class();
#include "style.h"
#undef BondClass

  else error->all("Invalid bond style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to current bond class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Bond *Force::bond_match(char *style)
{
  if (strcmp(bond_style,style) == 0) return bond;
  else if (strcmp(bond_style,"hybrid") == 0) {
    BondHybrid *hbond = (BondHybrid *) bond;
    for (int i = 0; i < hbond->nstyles; i++)
      if (strcmp(hbond->keywords[i],style) == 0) return hbond->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create an angle style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_angle(char *style)
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

Angle *Force::new_angle(char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define AngleClass
#define AngleStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class();
#include "style.h"
#undef AngleClass

  else error->all("Invalid angle style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create a dihedral style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_dihedral(char *style)
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

Dihedral *Force::new_dihedral(char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define DihedralClass
#define DihedralStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class();
#include "style.h"
#undef DihedralClass

  else error->all("Invalid dihedral style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create an improper style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_improper(char *style)
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

Improper *Force::new_improper(char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define ImproperClass
#define ImproperStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class();
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
  else if (strcmp(arg[0],#key) == 0) kspace = new Class(narg-1,&arg[1]);
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

  if (narg == 1 && strcmp(arg[0],"charmm") == 0) {
    special_lj[1] = 0.0;
    special_lj[2] = 0.0;
    special_lj[3] = 0.0;
    special_coul[1] = 0.0;
    special_coul[2] = 0.0;
    special_coul[3] = 0.0;
  } else if (narg == 1 && strcmp(arg[0],"amber") == 0) {
    special_lj[1] = 0.0;
    special_lj[2] = 0.0;
    special_lj[3] = 0.5;
    special_coul[1] = 0.0;
    special_coul[2] = 0.0;
    special_coul[3] = 5.0/6.0;
  } else if (narg == 3) {
    special_lj[1] = special_coul[1] = atof(arg[0]);
    special_lj[2] = special_coul[2] = atof(arg[1]);
    special_lj[3] = special_coul[3] = atof(arg[2]);
  } else if (narg == 6) {
    special_lj[1] = atof(arg[0]);
    special_lj[2] = atof(arg[1]);
    special_lj[3] = atof(arg[2]);
    special_coul[1] = atof(arg[3]);
    special_coul[2] = atof(arg[4]);
    special_coul[3] = atof(arg[5]);
  } else error->all("Illegal special_bonds command");
}

/* ----------------------------------------------------------------------
   create a new temperature
   called from input script or fixes that create temperature
   if flag, then allow overwrite of existing temperature with same ID
------------------------------------------------------------------------- */

void Force::add_temp(int narg, char **arg, int flag)
{
  if (narg < 3) error->all("Illegal temperature command");

  // error checks

  int itemp;
  for (itemp = 0; itemp < ntemp; itemp++)
    if (strcmp(arg[0],templist[itemp]->id) == 0) break;
  if (flag == 0 && itemp < ntemp) error->all("Reuse of temperature ID");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find temperature group ID");

  // if resetting existing temperature, delete it first

  int newflag = 1;
  if (itemp < ntemp) {
    newflag = 0;
    delete templist[itemp];
  }

  // extend Temperature list if necessary

  if (itemp == maxtemp) {
    maxtemp += DELTA;
    templist = (Temperature **) 
      memory->srealloc(templist,maxtemp*sizeof(Temperature *),
		       "modify:templist");
  }

  // create the Temperature class

  if (0) return;         // dummy line to enable else-if macro expansion

#define TempClass
#define TempStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) templist[itemp] = new Class(narg,arg);
#include "style.h"
#undef TempClass

  else error->all("Invalid temperature style");

  if (newflag) ntemp++;
}

/* ----------------------------------------------------------------------
   return ptr to temperature that matches ID
   else return NULL
------------------------------------------------------------------------- */

Temperature *Force::find_temp(char *id)
{
  for (int itemp = 0; itemp < ntemp; itemp++)
    if (strcmp(id,templist[itemp]->id) == 0) return templist[itemp];
  return NULL;
}

/* ----------------------------------------------------------------------
   modify temperature parameters
------------------------------------------------------------------------- */

void Force::modify_temp(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal temperature_modify command");

  // lookup Temperature ID

  int itemp;
  for (itemp = 0; itemp < ntemp; itemp++)
    if (strcmp(arg[0],templist[itemp]->id) == 0) break;
  if (itemp == ntemp) error->all("Could not find temp_modify ID");
  
  templist[itemp]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   nmax = upper bound
   5 possibilities: (1) i = i to i, (2) * = 1 to nmax,
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
   memory usage of force classes
------------------------------------------------------------------------- */

int Force::memory_usage()
{
  int bytes = 0;
  if (pair) bytes += pair->memory_usage();
  if (bond) bytes += bond->memory_usage();
  if (angle) bytes += angle->memory_usage();
  if (dihedral) bytes += dihedral->memory_usage();
  if (improper) bytes += improper->memory_usage();
  if (kspace) bytes += kspace->memory_usage();
  return bytes;
}
