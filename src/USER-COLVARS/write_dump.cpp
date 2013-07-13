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
#include "dump.h"
#include "style_dump.h"
#include "atom.h"
#include "group.h"
#include "error.h"

#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

// syntax: write_dump <group-ID> <dump-style> <filename> [<flags>]
//  flags are keywords supported by the dump command

void WriteDump::command(int narg, char **arg)
{
  Dump *dump;
  int lnarg;
  char **larg;

  if (narg < 3) error->all(FLERR,"Illegal write_dump command");

  if (atom->tag_enable == 0)
      error->all(FLERR,"Must have atom IDs for write_dump command");

  // create the Dump class instance
  lnarg = narg + 2;
  larg = new char*[lnarg];
  larg[0] = (char *) "WRITE_DUMP";
  larg[1] = arg[0]; // group
  larg[2] = arg[1]; // dump style
  larg[3] = (char *) "0";
  for (int i=2; i < narg; ++i)
    larg[i+2] = arg[i];

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) dump = new Class(lmp,lnarg,larg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  // tweak dump settings. write out sorted dump
  lnarg = 2;
  larg[0] = (char *) "sort";
  larg[1] = (char *) "id";
  dump->modify_params(lnarg,larg);

  // write out one frame and then delete the dump again
  dump->init();
  dump->write();
  delete dump;
  delete[] larg;
}

