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
#include "atom.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "python", "(unknown)"};

enum{INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
     SCALARFILE,ATOMFILE,FORMAT,EQUAL,ATOM,PYTHON};

/* ---------------------------------------------------------------------- */

void Info::command(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal info command");

  if (!screen) return;

  if (strcmp(arg[0],"groups") == 0) {
    int ngroup = group->ngroup;
    char **names = group->names;
    fprintf(screen,"Group information:\n");
    for (int i=0; i < ngroup; ++i) {
      fprintf(screen,"Group[%2d]: %s\n",i,names[i]);
    }

  } else if (strcmp(arg[0],"variables") == 0) {
    int nvar = input->variable->nvar;
    int *style = input->variable->style;
    char **names = input->variable->names;
    char ***data = input->variable->data;
    fprintf(screen,"Variable information:\n");
    for (int i=0; i < nvar; ++i) {
      fprintf(screen,"Variable[%3d]: %-10s  style = %-10s  def = %s\n",
             i,names[i],varstyles[style[i]],data[i][0]);
   
    }

  } else if (strcmp(arg[0],"system") == 0) {
    fprintf(screen,"System information:\n");
    fprintf(screen,"Units =  %s\n",update->unit_style);
    fprintf(screen,"Atom style = %s\n", atom->atom_style);
    fprintf(screen,"Natoms     = " BIGINT_FORMAT "\n", atom->natoms);
    fprintf(screen,"Nbonds     = " BIGINT_FORMAT "\n", atom->nbonds);
    fprintf(screen,"Nangles    = " BIGINT_FORMAT "\n", atom->nangles);
    fprintf(screen,"Ndihedrals = " BIGINT_FORMAT "\n", atom->ndihedrals);
    fprintf(screen,"Nimpropers = " BIGINT_FORMAT "\n", atom->nimpropers);

  } else {
    error->all(FLERR,"Unknown info command style");
  }
}
