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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "string.h"
#include "info.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void Info::command(int narg, char **arg)
{
  if ((narg != 1) && (narg !=3)) error->all(FLERR,"Illegal info command");

  char *varname = NULL;
  if (narg == 3) {
    if (strcmp(arg[1],"variable") == 0) {
      varname = arg[2];
    } else error->all(FLERR,"Illegal info command");
  }

  if (strcmp(arg[0],"groups") == 0) {
    int ngroup = group->ngroup;
    char **names = group->names;
    for (int i=0; i < ngroup; ++i) {
      if (screen) fprintf(screen,"group[%d]: %s\n",i,names[i]);
      if (logfile) fprintf(logfile,"group[%d]: %s\n",i,names[i]);
    }

#if 0
    if (varname) {
      if (screen) fprintf(screen,"storing list of groups in index variable %s\n",varname);
      char **varcmd = new char*[ngroup+2];
      varcmd[0] = varname;
      varcmd[1] = (char *)"index";
      for (int i=1; i < ngroup-1; ++i) {
        varcmd[i+1] = names[i];
      }

      input->variable->set(ngroup+1,varcmd);
      delete[] varcmd;
    }
#endif
 
  } else if (strcmp(arg[0],"units") == 0) {
    if (screen) fprintf(screen,"units %s\n",update->unit_style);
    if (logfile) fprintf(logfile,"units %s\n",update->unit_style);

  } else {
    error->all(FLERR,"Unknown info command style");
  }
}
