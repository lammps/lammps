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

#include "write_dump.h"
#include "style_dump.h"
#include "dump.h"
#include "atom.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void WriteDump::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal write_dump command");

  if (atom->tag_enable == 0)
      error->all(FLERR,"Must have atom IDs for write_dump command");

  // modindex = index in args of "modify", narg if doesn't exist

  int modindex;
  for (modindex = 0; modindex < narg; modindex++)
    if (strcmp(arg[modindex],"modify") == 0) break;

  // create the Dump instance

  Dump *dump;

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) dump = new Class(lmp,modindex,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  // pass additional args to modify_params

  if (modindex < narg) dump->modify_params(narg-modindex-1,&arg[modindex+1]);

  // write the current snapshot to dump file

  dump->init();
  dump->write();

  // delete the Dump instance

  delete dump;
}

