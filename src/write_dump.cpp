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
#include "write_dump.h"
#include "style_dump.h"
#include "dump.h"
#include "dump_image.h"
#include "atom.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void WriteDump::command(int narg, char **arg)
{

  if (narg < 3) error->all(FLERR,"Illegal write_dump command");

  // modindex = index in args of "modify" keyword
  // will be narg if "modify" is not present

  int modindex;
  for (modindex = 0; modindex < narg; modindex++)
    if (strcmp(arg[modindex],"modify") == 0) break;

  // create the Dump instance
  // create dump command line with extra required args

  Dump *dump;

  char **dumpargs = new char*[modindex+2];
  dumpargs[0] = (char *) "WRITE_DUMP"; // dump id
  dumpargs[1] = arg[0];                // group
  dumpargs[2] = arg[1];                // dump style
  dumpargs[3] = (char *) "0";          // dump frequency

  for (int i = 2; i < modindex; ++i)
    dumpargs[i+2] = arg[i];

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) dump = new Class(lmp,modindex+2,dumpargs);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  if (modindex < narg) dump->modify_params(narg-modindex-1,&arg[modindex+1]);

  // write out one frame and then delete the dump again
  // set multifile_override for DumpImage so that filename needs no "*"

  if (strcmp(arg[1],"image") == 0) 
    ((DumpImage *) dump)->multifile_override = 1;

  if (strcmp(arg[1],"cfg") == 0)
    ((DumpCFG *) dump)->multifile_override = 1;

  dump->init();
  dump->write();

  // delete the Dump instance and local storage

  delete dump;
  delete[] dumpargs;
}

