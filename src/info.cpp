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
#include "force.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "python", "(unknown)"};

static const char *mapstyles[] = { "none", "array", "hash" };

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
      int ndata = 1;
      fprintf(screen,"Variable[%3d]: %-10s  style = %-10s  def =",
             i,names[i],varstyles[style[i]]);
      if ((style[i] != LOOP) && (style[i] != ULOOP))
        ndata = input->variable->num[i];
      for (int j=0; j < ndata; ++j)
        fprintf(screen," %s",data[i][j]);
      fputs("\n",screen);
    }

  } else if (strcmp(arg[0],"system") == 0) {
    fprintf(screen,"System information:\n");
    fprintf(screen,"Units      = %s\n",update->unit_style);
    fprintf(screen,"Atom style = %s\n", atom->atom_style);
    fprintf(screen,"Atom map   = %s\n", mapstyles[atom->map_style]);
    fprintf(screen,"Atoms     / types = " BIGINT_FORMAT " / %d\n",
            atom->natoms, atom->ntypes);

    if (atom->molecular > 0) {

      fprintf(screen,"Molecule type = %s\n",
              (atom->molecular == 2) ? "template" : "standard");
      fprintf(screen,"Bonds     / types = " BIGINT_FORMAT " / %d\n",
              atom->nbonds, atom->nbondtypes);
      fprintf(screen,"Angles    / types = " BIGINT_FORMAT " / %d\n",
              atom->nangles, atom->nangletypes);
      fprintf(screen,"Dihedrals / types = " BIGINT_FORMAT " / %d\n",
              atom->ndihedrals, atom->ndihedraltypes);
      fprintf(screen,"impropers / types = " BIGINT_FORMAT " / %d\n",
              atom->nimpropers, atom->nimpropertypes);

      const double * const special_lj   = force->special_lj;
      const double * const special_coul = force->special_coul;

      fprintf(screen,"Special bond factors lj:   %-10g %-10g %-10g\n"
                     "Special bond factors coul: %-10g %-10g %-10g\n",
                     special_lj[1],special_lj[2],special_lj[3],
                     special_coul[1],special_coul[2],special_coul[3]);
    }

  } else {
    error->all(FLERR,"Unknown info command style");
  }
}
