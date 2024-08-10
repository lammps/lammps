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

  // work around "fix not computed at compatible times" errors.

  int dumpfreq = MAX(1, update->nsteps);
  dumpfreq += update->ntimestep % dumpfreq;

  std::string dump_id = "WRITE_DUMP";
  auto dumpargs = new char *[modindex + 2];
  dumpargs[0] = (char *) dump_id.c_str();                   // dump id
  dumpargs[1] = arg[0];                                     // group
  dumpargs[2] = arg[1];                                     // dump style
  dumpargs[3] = utils::strdup(std::to_string(dumpfreq));    // dump frequency

  // copy arguments up to modify, but skip over "noinit" if present

  int noinitwarn = 0;
  for (int i = 2; i < modindex; ++i) {
    if (strcmp(arg[i], "noinit") == 0) {
      noinitwarn = 1;
    } else {
      dumpargs[i + 2 - noinitwarn] = arg[i];
    }
  }

  auto *dump = output->add_dump(modindex + 2 - noinitwarn, dumpargs);

  try {
    if (modindex < narg) dump->modify_params(narg - modindex - 1, &arg[modindex + 1]);
  } catch (LAMMPSException &e) {
    // delete dump after error and then rethrow the exception to avoid re-use of dump-ID error
    dump = output->get_dump_by_id(dump_id);
    if (dump) output->delete_dump(dump_id);
    throw e;
  }

  // write out one frame and then delete the dump again
  // set multifile_override for DumpImage so that filename needs no "*"

  if (strcmp(arg[1], "image") == 0) (dynamic_cast<DumpImage *>(dump))->multifile_override = 1;
  if (strcmp(arg[1], "cfg") == 0) (dynamic_cast<DumpCFG *>(dump))->multifile_override = 1;
  if ((update->first_update == 0) && (comm->me == 0) && (noinitwarn == 0))
    error->warning(FLERR, "Calling write_dump before a full system init.");

  dump->init();
  dump->write();

  // delete the Dump instance and local storage

  output->delete_dump(dump_id);
  delete[] dumpargs[3];
  delete[] dumpargs;
}
