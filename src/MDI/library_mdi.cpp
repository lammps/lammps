/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// ----------------------------------------------------------------------
// MolSSI Driver Interface functions
// ----------------------------------------------------------------------

#include "library_mdi.h"

// needed to enable MPI support
#define LAMMPS_LIB_MPI 1
#include "library.h"

#include "fix_mdi_engine.h"

#include <cstring>

using namespace LAMMPS_NS;

/** Initialize an instance of LAMMPS as an MDI plugin
 *
\verbatim embed:rst

This function is called by the MolSSI Driver Interface library (MDI)
when LAMMPS is run as a plugin, and should not otherwise be used.

The function initializes MDI, then creates and initializes an instance
of LAMMPS.  The command-line arguments ``argc`` and ``argv`` used to
initialize LAMMPS are recieved from MDI.  The LAMMPS instance runs an
input file, which must include the ``mdi/engine`` command; when LAMMPS
executes this command, it will begin listening for commands from the
driver.  The name of the input file is obtained from the ``-in``
command-line argument, which must be provided by the MDI driver.

\endverbatim
 * \param  command    string buffer corresponding to the command to be executed
 * \param  comm       MDI communicator that can be used to communicated with the driver.
 * \param  class_obj  pointer to an instance of an mdi/engine fix cast to ``void *``.
 * \return 0 on no error. */

int MDI_Plugin_init_lammps()
{
  // initialize MDI
  int mdi_argc;
  char **mdi_argv;
  if (MDI_Plugin_get_argc(&mdi_argc)) MPI_Abort(MPI_COMM_WORLD, 1);
  if (MDI_Plugin_get_argv(&mdi_argv)) MPI_Abort(MPI_COMM_WORLD, 1);
  if (MDI_Init(&mdi_argc, &mdi_argv)) MPI_Abort(MPI_COMM_WORLD, 1);

  // get the MPI intra-communicator for this code
  MPI_Comm mpi_world_comm = MPI_COMM_WORLD;
  if (MDI_MPI_get_world_comm(&mpi_world_comm)) MPI_Abort(MPI_COMM_WORLD, 1);

  // find the -in argument
  int iarg = 0;
  char *filename;
  bool found_filename = false;
  while (iarg < mdi_argc && !found_filename) {

    if ((strcmp(mdi_argv[iarg], "-in") == 0) || (strcmp(mdi_argv[iarg], "-i") == 0)) {

      if (iarg + 2 > mdi_argc) MPI_Abort(MPI_COMM_WORLD, 1);
      filename = mdi_argv[iarg + 1];
      found_filename = true;

      // remove -in argument from the command list
      mdi_argc -= 2;
      for (int jarg = iarg; jarg < mdi_argc; jarg++) mdi_argv[jarg] = mdi_argv[jarg + 2];
    }
    iarg++;
  }
  if (!found_filename) MPI_Abort(MPI_COMM_WORLD, 1);

  // create and run a LAMMPS instance
  void *lmp = nullptr;
  if (lammps_config_has_mpi_support() > 0)
    lmp = lammps_open(mdi_argc, mdi_argv, mpi_world_comm, nullptr);
  else
    lmp = lammps_open_no_mpi(mdi_argc, mdi_argv, nullptr);
  lammps_file(lmp, filename);
  lammps_close(lmp);

  return 0;
}

/** Execute an MDI command
 *
\verbatim embed:rst

This function is called by the MolSSI Driver Interface Library (MDI)
when LAMMPS is run as a plugin, and should not otherwise be used.
The function executes a single command from an external MDI driver.

\endverbatim
 * \param  command   string buffer corresponding to the command to be executed
 * \param  comm      MDI communicator that can be used to communicated with the driver.
 * \param  class_obj pointer to an instance of an mdi/engine fix cast to ``void *``.
 * \return 0 on no error, 1 on error. */
int lammps_execute_mdi_command(const char *command, MDI_Comm comm, void *class_obj)
{
  FixMDIEngine *mdi_fix = (FixMDIEngine *) class_obj;
  return mdi_fix->execute_command(command, comm);
}
