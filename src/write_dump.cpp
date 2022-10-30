/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "comm.h"
#include "dump.h"
#include "dump_cfg.h"
#include "dump_image.h"
#include "error.h"
#include "output.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void WriteDump::command(int narg, char **arg)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "write_dump", error);

  // modindex = index in args of "modify" keyword
  // will be narg if "modify" is not present

  int modindex;
  for (modindex = 0; modindex < narg; modindex++)
    if (strcmp(arg[modindex], "modify") == 0) break;

  // create the Dump instance
  // create dump command line with extra required args

  auto dumpargs = new char *[modindex + 2];
  dumpargs[0] = (char *) "WRITE_DUMP";                                       // dump id
  dumpargs[1] = arg[0];                                                      // group
  dumpargs[2] = arg[1];                                                      // dump style
  dumpargs[3] = utils::strdup(std::to_string(MAX(update->ntimestep, 1)));    // dump frequency

  for (int i = 2; i < modindex; ++i) dumpargs[i + 2] = arg[i];

  Dump *dump = output->add_dump(modindex + 2, dumpargs);
  if (modindex < narg) dump->modify_params(narg - modindex - 1, &arg[modindex + 1]);

  // write out one frame and then delete the dump again
  // set multifile_override for DumpImage so that filename needs no "*"

  if (strcmp(arg[1], "image") == 0) (dynamic_cast<DumpImage *>(dump))->multifile_override = 1;
  if (strcmp(arg[1], "cfg") == 0) (dynamic_cast<DumpCFG *>(dump))->multifile_override = 1;
  if ((update->first_update == 0) && (comm->me == 0))
    error->warning(FLERR, "Calling write_dump before a full system init.");

  dump->init();
  dump->write();

  // delete the Dump instance and local storage

  output->delete_dump(dumpargs[0]);
  delete[] dumpargs[3];
  delete[] dumpargs;
}
