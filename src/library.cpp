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

// C or Fortran style library interface to LAMMPS
// customize by adding new LAMMPS-specific functions

#include "library.h"
#include <mpi.h>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include "universe.h"
#include "atom_vec.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"
#include "compute.h"
#include "fix.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "fix_external.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#if defined(LAMMPS_EXCEPTIONS)
#include "exceptions.h"
#endif

using namespace LAMMPS_NS;

// ----------------------------------------------------------------------
// utility macros
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   macros for optional code path which captures all exceptions
   and stores the last error message. These assume there is a variable lmp
   which is a pointer to the current LAMMPS instance.

   Usage:

   BEGIN_CAPTURE
   {
     // code paths which might throw exception
     ...
   }
   END_CAPTURE
------------------------------------------------------------------------- */

#ifdef LAMMPS_EXCEPTIONS
#define BEGIN_CAPTURE \
  Error * error = lmp->error; \
  try

#define END_CAPTURE \
  catch(LAMMPSAbortException & ae) { \
    int nprocs = 0; \
    MPI_Comm_size(ae.universe, &nprocs ); \
    \
    if (nprocs > 1) { \
      error->set_last_error(ae.message.c_str(), ERROR_ABORT); \
    } else { \
      error->set_last_error(ae.message.c_str(), ERROR_NORMAL); \
    } \
  } catch(LAMMPSException & e) { \
    error->set_last_error(e.message.c_str(), ERROR_NORMAL); \
  }
#else
#define BEGIN_CAPTURE
#define END_CAPTURE
#endif

// ----------------------------------------------------------------------
// library API functions to create/destroy an instance of LAMMPS
//   and communicate commands to it
// ----------------------------------------------------------------------

/** \brief Create an instance of the LAMMPS class and store reference in ptr
 *
\verbatim embed:rst
The :cpp:func:`lammps_open` function will create a new
LAMMPS instance while passing in a list of strings as if they were
:doc:`command-line arguments <Run_options>` when LAMMPS is run in
stand-alone mode from the command line, and an MPI communicator for
LAMMPS to run under.  Since the list of arguments is **exactly** as
when called from the command line, the first argument would be the
name of the executable and thus is ignored.  However ``argc`` may
be set to 0 and then ``argv`` may be a ``NULL`` pointer.

If for some reason the initialization of the LAMMPS instance fails
the ``ptr`` handle will be set to a ``NULL`` pointer.
\endverbatim
 *
 * \param argc number of command line arguments
 * \param argv list of command line argument strings
 * \param communicator MPI communicator for this LAMMPS instance.
 * \param ptr pointer to a void pointer variable which serves as a handle
 */
void lammps_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
#ifdef LAMMPS_EXCEPTIONS
  try
  {
    LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
    *ptr = (void *) lmp;
  }
  catch(LAMMPSException & e) {
    fprintf(stderr, "LAMMPS Exception: %s", e.message.c_str());
    *ptr = (void *) NULL;
  }
#else
  LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lmp;
#endif
}

/** \brief Variant of lammps_open() that implicitly uses ``MPI_COMM_WORLD``
 *
\verbatim embed:rst
This function works like :cpp:func:`lammps_open`, but
does not require to pass an MPI communicator.  Instead it will use
``MPI_COMM_WORLD`` when creating a LAMMPS instance.  It first calls
``MPI_Init()`` if the MPI environment is not yet initialized.

Outside of the convenience, this function is useful, when the LAMMPS
library was compiled in serial mode with the MPI ``STUBS`` library,
but is called from a parallel MPI program.
\endverbatim
 *
 * \param argc number of command line arguments
 * \param argv list of command line argument strings
 * \param ptr pointer to a void pointer variable which serves as a handle
 */
void lammps_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

#ifdef LAMMPS_EXCEPTIONS
  try
  {
    LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
    *ptr = (void *) lmp;
  }
  catch(LAMMPSException & e) {
    fprintf(stderr, "LAMMPS Exception: %s", e.message.c_str());
    *ptr = (void*) NULL;
  }
#else
  LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lmp;
#endif
}

/** \brief Delete a LAMMPS instance created by lammps_open() or
 *   lammps_open_no_mpi()
 *
 * \param ptr pointer to a previously created LAMMPS instance.
 *
\verbatim embed:rst
This function deletes a LAMMPS class instance pointed to by ``ptr``.
It does **not** call ``MPI_Finalize()`` to allow creating and deleting
multiple LAMMPS instances.  See :cpp:func:`lammps_finalize` for a
function to call at the end of the program in order to close down the
MPI infrastructure in case the calling program is not using MPI and was
creating the LAMMPS instance through :c:func:`lammps_open_no_mpi`.
\endverbatim
*/
void lammps_close(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  delete lmp;
}

/** \brief Get the numerical representation of the current LAMMPS version.
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \return an integer representing the version data in the format YYYYMMDD

\verbatim embed:rst

The :cpp:func:`lammps_version` function can be used to
determine the specific version of the underlying LAMMPS code. This is
particularly useful when loading LAMMPS as a shared library via dlopen()
so that the calling code can handle changes in the C API, or differences
in behavior, or available features (exported functions).

The returned LAMMPS version code is an integer in the format YYYYMMDD
(e.g. a LAMMPS version string of 2 Sep 2015 results in 20150902) that
will grow with every new LAMMPS release and thus is suitable for simple
numerical comparisons, e. g. when a program calling the library interface
is meant to be compatible with multiple LAMMPS versions with changes
in the command syntax.
\endverbatim
 */

int lammps_version(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return atoi(lmp->universe->num_ver);
}

/** \brief Shut down the MPI infrastructure

\verbatim embed:rst
The MPI standard requires that any MPI application calls
``MPI_Finalize()`` before exiting.  If a calling program does
not do any MPI calls and creates the LAMMPS instance through
:cpp:func:`lammps_open_no_mpi`,
MPI is initialized implicitly using the ``MPI_COMM_WORLD``
communicator, then this function should then be called
before exiting the program to wait until all parallel
tasks are completed and then cleanly shut down MPI.
\endverbatim
*/
void lammps_finalize()
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

/** \brief free memory buffer allocated by LAMMPS

\verbatim embed:rst
Some of the library interface functions return data as pointer to a
buffer that has been allocated by LAMMPS or the library interface.
This function can be used to delete those in order to avoid memory leaks.
\endverbatim
 *
 * \param ptr pointer to data allocated by LAMMPS
*/
void lammps_free(void *ptr)
{
  free(ptr);
}

// ----------------------------------------------------------------------
// library API functions to process commands
// ----------------------------------------------------------------------

/** \brief Process LAMMPS input from a file
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param filename name of a file with LAMMPS input

\verbatim embed:rst
This function will process the commands in the file pointed to
by ``filename`` line by line and thus functions very much like
the :doc:`include <include>` command.  The function returns when
the end of the file is reached and the commands have completed.
This may take considerable amounts of time, if the input file
contains :doc:`run <run>` or :doc:`minimize <minimize>` commands.
\endverbatim
  */
void lammps_file(void *ptr, char *filename)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0)
      lmp->error->all(FLERR,"Library error: issuing LAMMPS command during run");
    else
      lmp->input->file(filename);
  }
  END_CAPTURE
}

/** \brief Process a single LAMMPS input command from a string
 *
\verbatim embed:rst
This function will process a single command in the string ``str``.

The string may but need not have a final newline character.
Newline characters in the body of the string will be treated as
part of the command and will not start a second command.  You may
use the function :cpp:func:`lammps_commands_string`
to process a string with multiple command lines.

The function returns the name of the command on success or ``NULL``
when passing a string without a command.
\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param str string with the LAMMPS command
 * \return string with command name or ``NULL``
 */
char *lammps_command(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  char *result = NULL;

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0)
      lmp->error->all(FLERR,"Library error: issuing LAMMPS command during run");
    else
      result = lmp->input->one(str);
  }
  END_CAPTURE

  return result;
}

/** \brief Process multiple LAMMPS input commands from list of strings
 *
\verbatim embed:rst
This function allows to process multiple commands from a list of strings.
It does so by first concatenating the individual strings in ``cmds``
into a single string,  while inserting newline characters, if needed.
The combined string would be equivalent to a multi-line chunk of an
input file and is then internally passed to
:cpp:func:`lammps_commands_string` for processing.
\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param ncmd number of strings in the list
 * \param cmds list of strings with the LAMMPS commands
 */
void lammps_commands_list(void *ptr, int ncmd, char **cmds)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  int n = ncmd+1;
  for (int i = 0; i < ncmd; i++) n += strlen(cmds[i]);

  char *str = (char *) lmp->memory->smalloc(n,"lib/commands/list:str");
  str[0] = '\0';
  n = 0;

  for (int i = 0; i < ncmd; i++) {
    strcpy(&str[n],cmds[i]);
    n += strlen(cmds[i]);
    if (str[n-1] != '\n') {
      str[n] = '\n';
      str[n+1] = '\0';
      n++;
    }
  }

  lammps_commands_string(ptr,str);
  lmp->memory->sfree(str);
}

/** \brief Process a block of LAMMPS input commands from a single string
 *
\verbatim embed:rst
This function will process a string similar to a block of commands from
an input file.  The string may have multiple lines (separated by newline
characters) and also single commands may be distributed over multiple
lines with continuation characters ('&').  Those lines are combined by
removing the '&' and the following newline character.  After this processing
the string is handed to LAMMPS for parsing and executing.

.. note::

   Multi-line commands enabled by triple quotes will NOT work with
   this function.

\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param str string with block of LAMMPS input
 */
void lammps_commands_string(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  // make copy of str so can strtok() it

  int n = strlen(str) + 1;
  char *copy = new char[n];
  strcpy(copy,str);

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0) {
      lmp->error->all(FLERR,"Library error: issuing LAMMPS command during run");
    }

    char *ptr = copy;
    for (int i=0; i < n-1; ++i) {

      // handle continuation character as last character in line or string
      if ((copy[i] == '&') && (copy[i+1] == '\n'))
        copy[i+1] = copy[i] = ' ';
      else if ((copy[i] == '&') && (copy[i+1] == '\0'))
        copy[i] = ' ';

      if (copy[i] == '\n') {
        copy[i] = '\0';
        lmp->input->one(ptr);
        ptr = copy + i+1;
      } else if (copy[i+1] == '\0')
        lmp->input->one(ptr);
    }
  }
  END_CAPTURE

  delete [] copy;
}

// ----------------------------------------------------------------------
// library API functions to extract info from LAMMPS or set info in LAMMPS
// ----------------------------------------------------------------------

/** \brief Query LAMMPS about global settings that can be expressed as an integer
 *
\verbatim embed:rst

This function will retrieve or compute global properties and return them
as an integer.
The following query keywords are currently supported. If a keyword is
not recognized, the function returns -1.

Please also see :cpp:func:`lammps_extract_global`.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Description / Return value
   * - bigint
     - size of the ``bigint`` integer type, 4 or 8 bytes.
       Set at :ref:`compile time <size>`.
   * - tagint
     - size of the ``tagint`` integer type, 4 or 8 bytes.
       Set at :ref:`compile time <size>`.
   * - imageint
     - size of the ``imageint`` integer type, 4 or 8 bytes.
       Set at :ref:`compile time <size>`.
   * - dimension
     - Number of dimensions: 2 or 3. See :doc:`dimension`.
   * - box_exist
     - 1 if the simulation box is defined, 0 if not. See :doc:`create_box`.
   * - triclinic
     - 1 if the the simulation box is triclinic, 0 if orthogonal. See :doc:`change_box`.
   * - nlocal
     - number of "owned" atoms of the current MPI rank.
   * - nghost
     - number of "ghost" atoms of the current MPI rank.
   * - nall
     - number of all "owned" and "ghost" atoms of the current MPI rank.
   * - nmax
     - maximum of nlocal+nghost across all MPI ranks (for per-atom data array size).
   * - ntypes
     - number of atom types
   * - molecule_flag
     - 1 if the atom style includes molecular topology data. See :doc:`atom_style`.
   * - q_flag
     - 1 if the atom style includes point charges. See :doc:`atom_style`.
   * - mu_flag
     - 1 if the atom style includes point dipoles. See :doc:`atom_style`.
   * - rmass_flag
     - 1 if the atom style includes per-atom masses, 0 if there are per-type masses. See :doc:`atom_style`.

\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param name string with keyword of the setting
 * \return the value of the queried setting as a signed integer or -1 if unknown
 */
int lammps_extract_setting(void * ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

// This can be customized by adding keywords and documenting them in the section above.
  if (strcmp(name,"bigint") == 0) return sizeof(bigint);
  if (strcmp(name,"tagint") == 0) return sizeof(tagint);
  if (strcmp(name,"imageint") == 0) return sizeof(imageint);

  if (strcmp(name,"dimension") == 0) return lmp->domain->dimension;
  if (strcmp(name,"box_exist") == 0) return lmp->domain->box_exist;
  if (strcmp(name,"triclinic") == 0) return lmp->domain->triclinic;

  if (strcmp(name,"nlocal") == 0) return lmp->atom->nlocal;
  if (strcmp(name,"nghost") == 0) return lmp->atom->nghost;
  if (strcmp(name,"nall") == 0) return lmp->atom->nlocal+lmp->atom->nghost;
  if (strcmp(name,"nmax") == 0) return lmp->atom->nmax;
  if (strcmp(name,"ntypes") == 0) return lmp->atom->ntypes;

  if (strcmp(name,"molecule_flag") == 0) return lmp->atom->molecule_flag;
  if (strcmp(name,"q_flag") == 0) return lmp->atom->q_flag;
  if (strcmp(name,"mu_flag") == 0) return lmp->atom->mu_flag;
  if (strcmp(name,"rmass_flag") == 0) return lmp->atom->rmass_flag;

  return -1;
}

/** \brief Get pointer to internal global LAMMPS variables or arrays.
 *
\verbatim embed:rst
This function returns a pointer to the location of some global
property stored in one of the constituent classes of a LAMMPS instance.
The returned pointer is cast to ``void *`` and needs to be cast to a
pointer of type that the entity represents.  The pointers returned
by this function are generally persistent unless, i.e. it is not
necessary to call the function again, unless a :doc:`clear` command
is issued and the contents of the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
class are wiped out and recreated.

Please also see :cpp:func:`lammps_extract_setting` and
:cpp:func:`lammps_extract_box`.

.. warning::

   Modifying the data in the location pointed to by the returned
   pointer may lead to inconsistent internal data and thus may
   cause failures or crashes or bogus simulations.  In general it
   is thus usually better to use a LAMMPS input command that sets
   or changes these parameters.  Those will takes care of all side
   effects and necessary updates of settings derived from such
   settings.  Where possible a reference to such a command or a
   relevant section of the manual is given below.

This table lists the supported names, their data types, length
of the data area, and a short description.  The ``bigint`` type
may be defined to be either an ``int`` or an ``int64_t``.  This
is selected at :ref:`compile time <size>`.


.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Length
     - Description
   * - units
     - char \*
     - 1
     - string with the current unit style. See :doc:`units`.
   * - dt
     - double
     - 1
     - length of the time step. See :doc:`timestep`.
   * - ntimestep
     - bigint
     - 1
     - current time step number. See :doc:`reset_timestep`.
   * - boxlo
     - double
     - 3
     - lower box boundaries. See :doc:`create_box`.
   * - boxhi
     - double
     - 3
     - upper box boundaries. See :doc:`create_box`.
   * - boxxlo
     - double
     - 1
     - lower box boundary in x-direction. See :doc:`create_box`.
   * - boxxhi
     - double
     - 1
     - upper box boundary in x-direction. See :doc:`create_box`.
   * - boxylo
     - double
     - 1
     - lower box boundary in y-direction. See :doc:`create_box`.
   * - boxyhi
     - double
     - 1
     - upper box boundary in y-direction. See :doc:`create_box`.
   * - boxzlo
     - double
     - 1
     - lower box boundary in z-direction. See :doc:`create_box`.
   * - boxzhi
     - double
     - 1
     - upper box boundary in z-direction. See :doc:`create_box`.
   * - periodicity
     - int
     - 3
     - 0 if non-periodic, 1 if periodic for x, y, and z. See :doc:`boundary`.
   * - triclinic
     - int
     - 1
     - 1 if the the simulation box is triclinic, 0 if orthogonal. See :doc:`change_box`.
   * - xy
     - double
     - 1
     - triclinic tilt factor. See :doc:`Howto_triclinic`.
   * - yz
     - double
     - 1
     - triclinic tilt factor. See :doc:`Howto_triclinic`.
   * - xz
     - double
     - 1
     - triclinic tilt factor. See :doc:`Howto_triclinic`.
   * - natoms
     - bigint
     - 1
     - total number of atoms in the simulation.
   * - nbonds
     - bigint
     - 1
     - total number of bonds in the simulation.
   * - nangles
     - bigint
     - 1
     - total number of angles in the simulation.
   * - ndihedrals
     - bigint
     - 1
     - total number of dihedrals in the simulation.
   * - nimpropers
     - bigint
     - 1
     - total number of impropers in the simulation.
   * - nlocal
     - int
     - 1
     - number of "owned" atoms of the current MPI rank.
   * - nghost
     - int
     - 1
     - number of "ghost" atoms of the current MPI rank.
   * - nmax
     - int
     - 1
     - maximum of nlocal+nghost across all MPI ranks (for per-atom data array size).
   * - ntypes
     - int
     - 1
     - number of atom types
   * - q_flag
     - int
     - 1
     - 1 if the atom style includes point charges. See :doc:`atom_style`.
   * - atime
     - double
     - 1
     - accumulated simulation time in time units.
   * - atimestep
     - bigint
     - 1
     - the number of the timestep when "atime" was last updated.
   * - boltz
     - double
     - 1
     - value of the "boltz" constant. See :doc:`units`.
   * - hplanck
     - double
     - 1
     - value of the "hplanck" constant. See :doc:`units`.
   * - mvv2e
     - double
     - 1
     - factor to convert :math:`\frac{1}{2}mv^2` for a particle to the current energy unit. See :doc:`units`.
   * - ftm2v
     - double
     - 1
     -  . See :doc:`units`.
   * - mv2d
     - double
     - 1
     -  . See :doc:`units`.
   * - nktv2p
     - double
     - 1
     -  . See :doc:`units`.
   * - qqr2e
     - double
     - 1
     - factor to convert :math:`\frac{q_i q_j}{r}` to energy units. See :doc:`units`.
   * - qe2f
     - double
     - 1
     -  . See :doc:`units`.
   * - vxmu2f
     - double
     - 1
     -  . See :doc:`units`.
   * - xxt2kmu
     - double
     - 1
     -  . See :doc:`units`.
   * - dielectric
     - double
     - 1
     - value of the dielectric constant. See :doc:`dielectric`.
   * - qqrd2e
     - double
     - 1
     -  . See :doc:`units`.
   * - e_mass
     - double
     - 1
     -  . See :doc:`units`.
   * - hhmrr2e
     - double
     - 1
     -  . See :doc:`units`.
   * - mvh2r
     - double
     - 1
     -  . See :doc:`units`.
   * - angstrom
     - double
     - 1
     - constant to convert current length unit to angstroms. 1.0 for reduced (aka "lj") units. See :doc:`units`.
   * - femtosecond
     - double
     - 1
     - constant to convert current time unit to femtoseconds. 1.0 for reduced (aka "lj") units
   * - qelectron
     - double
     - 1
     -  . See :doc:`units`.

\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param name string with name of the entity
 * \return pointer cast to ``void *`` to the location of the requested property. NULL if unknown keyword name
 */

void *lammps_extract_global(void *ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  if (strcmp(name,"units") == 0) return (void *) lmp->update->unit_style;
  if (strcmp(name,"dt") == 0) return (void *) &lmp->update->dt;
  if (strcmp(name,"ntimestep") == 0) return (void *) &lmp->update->ntimestep;
  if (strcmp(name,"boxlo") == 0) return (void *) lmp->domain->boxlo;
  if (strcmp(name,"boxhi") == 0) return (void *) lmp->domain->boxhi;
  if (strcmp(name,"boxxlo") == 0) return (void *) &lmp->domain->boxlo[0];
  if (strcmp(name,"boxxhi") == 0) return (void *) &lmp->domain->boxhi[0];
  if (strcmp(name,"boxylo") == 0) return (void *) &lmp->domain->boxlo[1];
  if (strcmp(name,"boxyhi") == 0) return (void *) &lmp->domain->boxhi[1];
  if (strcmp(name,"boxzlo") == 0) return (void *) &lmp->domain->boxlo[2];
  if (strcmp(name,"boxzhi") == 0) return (void *) &lmp->domain->boxhi[2];
  if (strcmp(name,"periodicity") == 0) return (void *) lmp->domain->periodicity;
  if (strcmp(name,"triclinic") == 0) return (void *) &lmp->domain->triclinic;
  if (strcmp(name,"xy") == 0) return (void *) &lmp->domain->xy;
  if (strcmp(name,"xz") == 0) return (void *) &lmp->domain->xz;
  if (strcmp(name,"yz") == 0) return (void *) &lmp->domain->yz;
  if (strcmp(name,"natoms") == 0) return (void *) &lmp->atom->natoms;
  if (strcmp(name,"nbonds") == 0) return (void *) &lmp->atom->nbonds;
  if (strcmp(name,"nangles") == 0) return (void *) &lmp->atom->nangles;
  if (strcmp(name,"ndihedrals") == 0) return (void *) &lmp->atom->ndihedrals;
  if (strcmp(name,"nimpropers") == 0) return (void *) &lmp->atom->nimpropers;
  if (strcmp(name,"nlocal") == 0) return (void *) &lmp->atom->nlocal;
  if (strcmp(name,"nghost") == 0) return (void *) &lmp->atom->nghost;
  if (strcmp(name,"nmax") == 0) return (void *) &lmp->atom->nmax;
  if (strcmp(name,"ntypes") == 0) return (void *) &lmp->atom->ntypes;

  if (strcmp(name,"q_flag") == 0) return (void *) &lmp->atom->q_flag;

  // update->atime can be referenced as a pointer
  // thermo "timer" data cannot be, since it is computed on request
  // lammps_get_thermo() can access all thermo keywords by value

  if (strcmp(name,"atime") == 0) return (void *) &lmp->update->atime;
  if (strcmp(name,"atimestep") == 0) return (void *) &lmp->update->atimestep;

  // global constants defined by units

  if (strcmp(name,"boltz") == 0) return (void *) &lmp->force->boltz;
  if (strcmp(name,"hplanck") == 0) return (void *) &lmp->force->hplanck;
  if (strcmp(name,"mvv2e") == 0) return (void *) &lmp->force->mvv2e;
  if (strcmp(name,"ftm2v") == 0) return (void *) &lmp->force->ftm2v;
  if (strcmp(name,"mv2d") == 0) return (void *) &lmp->force->mv2d;
  if (strcmp(name,"nktv2p") == 0) return (void *) &lmp->force->nktv2p;
  if (strcmp(name,"qqr2e") == 0) return (void *) &lmp->force->qqr2e;
  if (strcmp(name,"qe2f") == 0) return (void *) &lmp->force->qe2f;
  if (strcmp(name,"vxmu2f") == 0) return (void *) &lmp->force->vxmu2f;
  if (strcmp(name,"xxt2kmu") == 0) return (void *) &lmp->force->xxt2kmu;
  if (strcmp(name,"dielectric") == 0) return (void *) &lmp->force->dielectric;
  if (strcmp(name,"qqrd2e") == 0) return (void *) &lmp->force->qqrd2e;
  if (strcmp(name,"e_mass") == 0) return (void *) &lmp->force->e_mass;
  if (strcmp(name,"hhmrr2e") == 0) return (void *) &lmp->force->hhmrr2e;
  if (strcmp(name,"mvh2r") == 0) return (void *) &lmp->force->mvh2r;

  if (strcmp(name,"angstrom") == 0) return (void *) &lmp->force->angstrom;
  if (strcmp(name,"femtosecond") == 0) return (void *) &lmp->force->femtosecond;
  if (strcmp(name,"qelectron") == 0) return (void *) &lmp->force->qelectron;

  return NULL;
}

/** \brief Extract simulation box parameters

\verbatim embed:rst
This function will (re-)initialize the simulation box and boundary information
and then assign the designated data to the locations in the pointers passed as
arguments.
\endverbatim

 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``
 * \param boxlo pointer to 3 doubles where the lower box boundary is stored
 * \param boxhi pointer to 3 doubles where the upper box boundary is stored
 * \param xy pointer to a double where the xy tilt factor is stored
 * \param yz pointer to a double where the yz tilt factor is stored
 * \param xz pointer to a double where the xz tilt factor is stored
 * \param periodicity pointer to 3 ints, set to 1 for periodic boundaries and 0 for non-periodic
 * \param box_change pointer to an int, which is set to 1 if the box will be
 *        changed during a simulation by a fix and 0 if not.
 */
void lammps_extract_box(void *ptr, double *boxlo, double *boxhi,
                        double *xy, double *yz, double *xz,
                        int *periodicity, int *box_change)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  Domain *domain = lmp->domain;

  // domain->init() is needed to set box_change
  domain->init();

  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
  boxhi[0] = domain->boxhi[0];
  boxhi[1] = domain->boxhi[1];
  boxhi[2] = domain->boxhi[2];

  *xy = domain->xy;
  *yz = domain->yz;
  *xz = domain->xz;

  periodicity[0] = domain->periodicity[0];
  periodicity[1] = domain->periodicity[1];
  periodicity[2] = domain->periodicity[2];

  *box_change = domain->box_change;
}


/** \brief Get pointer to a LAMMPS per-atom property.
 *
\verbatim embed:rst
This function returns a pointer to the location of per-atom properties
(and per-atomtype properties as in the case of the 'mass' keyword).
This is data that is distributed across sub-domains and thus MPI ranks.
The returned pointer is cast to ``void *`` and needs to be cast to a
pointer of type that the entity represents.

.. note::

   The pointers returned by this function are generally not persistent
   since per-atom data may be re-distributed, re-allocated, and
   re-ordered at every re-neighboring operation.

This table lists a large part of the supported names, their data types,
length of the data area, and a short description.  You can look up
additional supported keywords and their data types in
:ref:`Atom::extract() <cpp_atom_extract>` and the ``src/atom.h`` header file.
The ``bigint``, ``tagint`, or ``imageint`` types may be defined to be
either an ``int`` or an ``int64_t``.  This is selected at
:ref:`compile time <size>`.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Items per atom
     - Description
   * - mass
     - double
     - 1
     - per-type mass. This array is **NOT** a per-atom array but of length ``ntypes+1``, element 0 is ignored.
   * - id
     - tagint
     - 1
     - atom ID of the particles
   * - type
     - int
     - 1
     - atom type of the particles
   * - mask
     - int
     - 1
     - bitmask for mapping to groups. Individual bits are set to 0 or 1 for each group.
   * - image
     - imageint
     - 1
     - 3 image flags encoded into a single integer using either 10 or 21 bits per direction.
   * - x
     - double
     - 3
     - x-, y-, and z-coordinate of the particles
   * - v
     - double
     - 3
     - x-, y-, and z-component of the velocity of the particles
   * - f
     - double
     - 3
     - x-, y-, and z-component of the force on the particles
   * - molecule
     - int
     - 1
     - molecule ID of the particles
   * - q
     - double
     - 1
     - charge of the particles
   * - mu
     - double
     - 3
     - dipole moment of the particles
   * - omega
     - double
     - 3
     - x-, y-, and z-component of rotational velocity of the particles
   * - angmom
     - double
     - 3
     - x-, y-, and z-component of angular momentum of the particles
   * - torque
     - double
     - 3
     - x-, y-, and z-component of the torque on the particles
   * - radius
     - double
     - 1
     - radius of the (extended) particles
   * - rmass
     - double
     - 1
     - per-atom mass of the particles. ``NULL`` if per-type masses are used. See :cpp:func:`'rmass_flag' setting <lammps_extract_setting>`.
   * - ellipsoid
     - int
     - 1
     - 1 if the particle is an ellipsoidal particle, 0 if not
   * - line
     - int
     - 1
     - 1 if the particle is a line particle, 0 if not
   * - tri
     - int
     - 1
     - 1 if the particle is a triangulated particle, 0 if not
   * - body
     - int
     - 1
     - 1 if the particle is a body particle, 0 if not

\endverbatim
 *
 * \param ptr pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param name string with name of the entity
 * \return pointer cast to ``void *`` to the location of the requested data or NULL if not found
 */
void *lammps_extract_atom(void *ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->extract(name);
}
/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS atom-based entity
   name = desired quantity, e.g. x or mass
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if Atom::extract() does not recognize the name
   customize by adding names to Atom::extract()
------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS compute-based entity
   the compute is invoked if its value(s) is not current
   id = compute ID
   style = 0 for global data, 1 for per-atom data, 2 for local data
   type = 0 for scalar, 1 for vector, 2 for array
   for global data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a scalar or vector
     caller should cast it to (double **) for an array
   for per-atom or local vector/array data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a vector
     caller should cast it to (double **) for an array
   for local data, accessing scalar data for the compute (type = 0),
   returns a pointer that should be cast to (int *) which points to
   an int with the number of local rows, i.e. the length of the local array.
   returns a void pointer to the compute's internal data structure
     for the entity which the caller can cast to the proper data type
   returns a NULL if id is not recognized or style/type not supported
   the returned pointer is not a permanent valid reference to the
     compute data, this function should be re-invoked
   IMPORTANT: if the compute is not current it will be invoked
     LAMMPS cannot easily check here if it is valid to invoke the compute,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_compute(void *ptr, char *id, int style, int type)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int icompute = lmp->modify->find_compute(id);
    if (icompute < 0) return NULL;
    Compute *compute = lmp->modify->compute[icompute];

    if (style == 0) {
      if (type == 0) {
        if (!compute->scalar_flag) return NULL;
        if (compute->invoked_scalar != lmp->update->ntimestep)
          compute->compute_scalar();
        return (void *) &compute->scalar;
      }
      if (type == 1) {
        if (!compute->vector_flag) return NULL;
        if (compute->invoked_vector != lmp->update->ntimestep)
          compute->compute_vector();
        return (void *) compute->vector;
      }
      if (type == 2) {
        if (!compute->array_flag) return NULL;
        if (compute->invoked_array != lmp->update->ntimestep)
          compute->compute_array();
        return (void *) compute->array;
      }
    }

    if (style == 1) {
      if (!compute->peratom_flag) return NULL;
      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();
      if (type == 1) return (void *) compute->vector_atom;
      if (type == 2) return (void *) compute->array_atom;
    }

    if (style == 2) {
      if (!compute->local_flag) return NULL;
      if (compute->invoked_local != lmp->update->ntimestep)
        compute->compute_local();
      if (type == 0) return (void *) &compute->size_local_rows;
      if (type == 1) return (void *) compute->vector_local;
      if (type == 2) return (void *) compute->array_local;
    }
  }
  END_CAPTURE

  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS fix-based entity
   id = fix ID
   style = 0 for global data, 1 for per-atom data, 2 for local data
   type = 0 for scalar, 1 for vector, 2 for array
   i,j = indices needed only to specify which global vector or array value
   for global data, returns a pointer to a memory location
     which is allocated by this function
     which the caller can cast to a (double *) which points to the value
   for per-atom or local data, returns a pointer to the
     fix's internal data structure for the entity
     caller should cast it to (double *) for a vector
     caller should cast it to (double **) for an array
   returns a NULL if id is not recognized or style/type not supported
   IMPORTANT: for global data,
     this function allocates a double to store the value in,
     so the caller must free this memory to avoid a leak, e.g.
       double *dptr = (double *) lammps_extract_fix();
       double value = *dptr;
       lammps_free(dptr);
   IMPORTANT: LAMMPS cannot easily check here when info extracted from
     the fix is valid, so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_fix(void *ptr, char *id, int style, int type,
                         int i, int j)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int ifix = lmp->modify->find_fix(id);
    if (ifix < 0) return NULL;
    Fix *fix = lmp->modify->fix[ifix];

    if (style == 0) {
      if (type == 0) {
        if (!fix->scalar_flag) return NULL;
        double *dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_scalar();
        return (void *) dptr;
      }
      if (type == 1) {
        if (!fix->vector_flag) return NULL;
        double *dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_vector(i);
        return (void *) dptr;
      }
      if (type == 2) {
        if (!fix->array_flag) return NULL;
        double *dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_array(i,j);
        return (void *) dptr;
      }
    }

    if (style == 1) {
      if (!fix->peratom_flag) return NULL;
      if (type == 1) return (void *) fix->vector_atom;
      if (type == 2) return (void *) fix->array_atom;
    }

    if (style == 2) {
      if (!fix->local_flag) return NULL;
      if (type == 0) return (void *) &fix->size_local_rows;
      if (type == 1) return (void *) fix->vector_local;
      if (type == 2) return (void *) fix->array_local;
    }
  }
  END_CAPTURE

  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS evaluated variable
   name = variable name, must be equal-style or atom-style variable
   group = group ID for evaluating an atom-style variable, else NULL
   for equal-style variable, returns a pointer to a memory location
     which is allocated by this function
     which the caller can cast to a (double *) which points to the value
   for atom-style variable, returns a pointer to the
     vector of per-atom values on each processor,
     which the caller can cast to a (double *) which points to the values
   returns a NULL if name is not recognized or not equal-style or atom-style
   IMPORTANT: for both equal-style and atom-style variables,
     this function allocates memory to store the variable data in
     so the caller must free this memory to avoid a leak
     e.g. for equal-style variables
       double *dptr = (double *) lammps_extract_variable();
       double value = *dptr;
       lammps_free(dptr);
     e.g. for atom-style variables
       double *vector = (double *) lammps_extract_variable();
       use the vector values
       lammps_free(vector);
   IMPORTANT: LAMMPS cannot easily check here when it is valid to evaluate
     the variable or any fixes or computes or thermodynamic info it references,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_variable(void *ptr, char *name, char *group)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int ivar = lmp->input->variable->find(name);
    if (ivar < 0) return NULL;

    if (lmp->input->variable->equalstyle(ivar)) {
      double *dptr = (double *) malloc(sizeof(double));
      *dptr = lmp->input->variable->compute_equal(ivar);
      return (void *) dptr;
    }

    if (lmp->input->variable->atomstyle(ivar)) {
      int igroup = lmp->group->find(group);
      if (igroup < 0) return NULL;
      int nlocal = lmp->atom->nlocal;
      double *vector = (double *) malloc(nlocal*sizeof(double));
      lmp->input->variable->compute_atom(ivar,igroup,vector,1,0);
      return (void *) vector;
    }
  }
  END_CAPTURE

  return NULL;
}

/* ----------------------------------------------------------------------
   return the current value of a thermo keyword as a double
   unlike lammps_extract_global() this does not give access to the
     storage of the data in question
   instead it triggers the Thermo class to compute the current value
     and returns it
------------------------------------------------------------------------- */

double lammps_get_thermo(void *ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  double dval = 0.0;

  BEGIN_CAPTURE
  {
    lmp->output->thermo->evaluate_keyword(name,&dval);
  }
  END_CAPTURE

  return dval;
}

/** \brief Return the total number of atoms in the system

\verbatim embed:rst
This is particularly useful before making calls to
:cpp:func:`lammps_get_natoms` so one can
pre-allocate the correct amount of storage for the result vector.

.. note::

   This function returns a 32-bit signed integer and thus will
   not work for systems with more than about 2 billion atoms.

\endverbatim

 * \param ptr pointer to a previously created LAMMPS instance
 * \return total number of atoms in the system or 0 if value too large.
 */
int lammps_get_natoms(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  if (lmp->atom->natoms > MAXSMALLINT) return 0;
  int natoms = static_cast<int> (lmp->atom->natoms);
  return natoms;
}

/* ----------------------------------------------------------------------
   set the value of a STRING variable to str
   return -1 if variable doesn't exist or not a STRING variable
   return 0 for success
------------------------------------------------------------------------- */

int lammps_set_variable(void *ptr, char *name, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  int err = -1;

  BEGIN_CAPTURE
  {
    err = lmp->input->variable->set_string(name,str);
  }
  END_CAPTURE

  return err;
}

/* ----------------------------------------------------------------------
   reset simulation box parameters
   see domain.h for definition of these arguments
   assumes domain->set_initial_box() has been invoked previously
------------------------------------------------------------------------- */

void lammps_reset_box(void *ptr, double *boxlo, double *boxhi,
                      double xy, double yz, double xz)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  Domain *domain = lmp->domain;

  domain->boxlo[0] = boxlo[0];
  domain->boxlo[1] = boxlo[1];
  domain->boxlo[2] = boxlo[2];
  domain->boxhi[0] = boxhi[0];
  domain->boxhi[1] = boxhi[1];
  domain->boxhi[2] = boxhi[2];

  domain->xy = xy;
  domain->yz = yz;
  domain->xz = xz;

  domain->set_global_box();
  lmp->comm->set_proc_grid();
  domain->set_local_box();
}

/* ----------------------------------------------------------------------
   gather the named atom-based entity for all atoms
     return it in user-allocated data
   data will be ordered by atom ID
     requirement for consecutive atom IDs (1 to N)
   see gather_atoms_concat() to return data for all atoms, unordered
   see gather_atoms_subset() to return data for only a subset of atoms
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   return atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Natoms, as queried by get_natoms()
   method:
     alloc and zero count*Natom length vector
     loop over Nlocal to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms(void *ptr, char * /*name */,
                         int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms(void *ptr, char *name,
                         int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,offset;

    // error if tags are not defined or not consecutive
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lmp->atom->extract(name);
    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0;

      tagint *tag = lmp->atom->tag;
      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < nlocal; i++)
          copy[tag[i]-1] = vector[i];

      } else if (imgunpack) {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          const int image = vector[i];
          copy[offset++] = (image & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
        }

      } else {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          for (j = 0; j < count; j++)
            copy[offset++] = array[i][j];
        }
      }

      MPI_Allreduce(copy,data,count*natoms,MPI_INT,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      double *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0.0;

      tagint *tag = lmp->atom->tag;
      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < nlocal; i++)
          copy[tag[i]-1] = vector[i];

      } else {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          for (j = 0; j < count; j++)
            copy[offset++] = array[i][j];
        }
      }

      MPI_Allreduce(copy,data,count*natoms,MPI_DOUBLE,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   gather the named atom-based entity for all atoms
     return it in user-allocated data
   data will be a concatenation of chunks of each proc's atoms,
     in whatever order the atoms are on each proc
     no requirement for consecutive atom IDs (1 to N)
     can do a gather_atoms_concat for "id" if need to know atom IDs
   see gather_atoms() to return data ordered by consecutive atom IDs
   see gather_atoms_subset() to return data for only a subset of atoms
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   return atom-based values in 1d data, ordered by count, then by atom
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Natoms, as queried by get_natoms()
   method:
     Allgather Nlocal atoms from each proc into data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms_concat(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_concat() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms_concat(void *ptr, char *name,
                                int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,offset;

    // error if tags are not defined
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lmp->atom->extract(name);
    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // perform MPI_Allgatherv on each proc's chunk of Nlocal atoms

    int nprocs = lmp->comm->nprocs;

    int *recvcounts,*displs;
    lmp->memory->create(recvcounts,nprocs,"lib/gather:recvcounts");
    lmp->memory->create(displs,nprocs,"lib/gather:displs");

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(vector,nlocal,MPI_INT,data,recvcounts,displs,
                       MPI_INT,lmp->world);

      } else if (imgunpack) {
        int *copy;
        lmp->memory->create(copy,count*nlocal,"lib/gather:copy");
        offset = 0;
        for (i = 0; i < nlocal; i++) {
          const int image = vector[i];
          copy[offset++] = (image & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
        }
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(copy,count*nlocal,MPI_INT,
                       data,recvcounts,displs,MPI_INT,lmp->world);
        lmp->memory->destroy(copy);

      } else {
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(&array[0][0],count*nlocal,MPI_INT,
                       data,recvcounts,displs,MPI_INT,lmp->world);
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(vector,nlocal,MPI_DOUBLE,data,recvcounts,displs,
                       MPI_DOUBLE,lmp->world);

      } else {
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(&array[0][0],count*nlocal,MPI_DOUBLE,
                       data,recvcounts,displs,MPI_DOUBLE,lmp->world);
      }
    }

    lmp->memory->destroy(recvcounts);
    lmp->memory->destroy(displs);
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   gather the named atom-based entity for a subset of atoms
     return it in user-allocated data
   data will be ordered by requested atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see gather_atoms() to return data for all atoms, ordered by consecutive IDs
   see gather_atoms_concat() to return data for all atoms, unordered
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   ndata = # of atoms to return data for (could be all atoms)
   ids = list of Ndata atom IDs to return data for
   return atom-based values in 1d data, ordered by count, then by atom
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Ndata
   method:
     alloc and zero count*Ndata length vector
     loop over Ndata to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms_subset(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/,
                                int /*ndata*/, int * /*ids*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms_subset(void *ptr, char *name,
                                int type, int count,
                                int ndata, int *ids, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms_subset");
      return;
    }

    void *vptr = lmp->atom->extract(name);
    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms_subset: "
                          "unknown property name");
      return;
    }

    // copy = Ndata length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*ndata,"lib/gather:copy");
      for (i = 0; i < count*ndata; i++) copy[i] = 0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal)
            copy[i] = vector[m];
        }

      } else if (imgunpack) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            const int image = vector[m];
            copy[offset++] = (image & IMGMASK) - IMGMAX;
            copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
            copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
          }
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            for (j = 0; j < count; j++)
              copy[offset++] = array[m][j];
          }
        }
      }

      MPI_Allreduce(copy,data,count*ndata,MPI_INT,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      double *copy;
      lmp->memory->create(copy,count*ndata,"lib/gather:copy");
      for (i = 0; i < count*ndata; i++) copy[i] = 0.0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal)
            copy[i] = vector[m];
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            for (j = 0; j < count; j++)
              copy[offset++] = array[m][j];
          }
        }
      }

      MPI_Allreduce(copy,data,count*ndata,MPI_DOUBLE,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   scatter the named atom-based entity in data to all atoms
   data is ordered by atom ID
     requirement for consecutive atom IDs (1 to N)
   see scatter_atoms_subset() to scatter data for some (or all) atoms, unordered
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" for xyz to be packed into single image flag
   data = atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be correct length = count*Natoms, as queried by get_natoms()
   method:
     loop over Natoms, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_scatter_atoms(void *ptr, char * /*name */,
                          int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_scatter_atoms(void *ptr, char *name,
                          int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;

    // error if tags are not defined or not consecutive or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == 0) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lmp->atom->extract(name);
    if(vptr == NULL) {
        lmp->error->warning(FLERR,
                            "lammps_scatter_atoms: unknown property name");
        return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgpack) vector = (int *) vptr;
      else array = (int **) vptr;
      int *dptr = (int *) data;

      if (count == 1) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0)
            vector[m] = dptr[i];

      } else if (imgpack) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            int image = dptr[offset++] + IMGMAX;
            image += (dptr[offset++] + IMGMAX) << IMGBITS;
            image += (dptr[offset++] + IMGMAX) << IMG2BITS;
            vector[m] = image;
          }

      } else {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      double *dptr = (double *) data;

      if (count == 1) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0)
            vector[m] = dptr[i];

      } else {
        for (i = 0; i < natoms; i++) {
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   scatter the named atom-based entity in data to a subset of atoms
   data is ordered by provided atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see scatter_atoms() to scatter data for all atoms, ordered by consecutive IDs
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" for xyz to be packed into single image flag
   ndata = # of atoms in ids and data (could be all atoms)
   ids = list of Ndata atom IDs to scatter data to
   data = atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be correct length = count*Ndata
   method:
     loop over Ndata, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_scatter_atoms_subset(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/,
                                int /*ndata*/, int * /*ids*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_scatter_atoms_subset(void *ptr, char *name,
                                 int type, int count,
                                 int ndata, int *ids, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == 0) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms_subset");
      return;
    }

    void *vptr = lmp->atom->extract(name);
    if(vptr == NULL) {
        lmp->error->warning(FLERR,
                            "lammps_scatter_atoms_subset: unknown property name");
        return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgpack) vector = (int *) vptr;
      else array = (int **) vptr;
      int *dptr = (int *) data;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0)
            vector[m] = dptr[i];
        }

      } else if (imgpack) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            int image = dptr[offset++] + IMGMAX;
            image += (dptr[offset++] + IMGMAX) << IMGBITS;
            image += (dptr[offset++] + IMGMAX) << IMG2BITS;
            vector[m] = image;
          }
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      double *dptr = (double *) data;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0)
            vector[m] = dptr[i];
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   create N atoms and assign them to procs based on coords
   id = atom IDs (optional, NULL will generate 1 to N)
   type = N-length vector of atom types (required)
   x = 3N-length 1d vector of atom coords (required)
   v = 3N-length 1d vector of atom velocities (optional, NULL if just 0.0)
   image flags can be treated in two ways:
     (a) image = vector of current image flags
         each atom will be remapped into periodic box by domain->ownatom()
         image flag will be incremented accordingly and stored with atom
     (b) image = NULL
         each atom will be remapped into periodic box by domain->ownatom()
         image flag will be set to 0 by atom->avec->create_atom()
   shrinkexceed = 1 allows atoms to be outside a shrinkwrapped boundary
     passed to ownatom() which will assign them to boundary proc
     important if atoms may be (slightly) outside non-periodic dim
     e.g. due to restoring a snapshot from a previous run and previous box
   id and image must be 32-bit integers
   x,v = ordered by xyz, then by atom
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
------------------------------------------------------------------------- */

void lammps_create_atoms(void *ptr, int n, tagint *id, int *type,
                         double *x, double *v, imageint *image,
                         int shrinkexceed)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    // error if box does not exist or tags not defined

    int flag = 0;
    if (lmp->domain->box_exist == 0) flag = 1;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_create_atoms");
      return;
    }

    // loop over N atoms of entire system
    // if this proc owns it based on coords, invoke create_atom()
    // optionally set atom tags and velocities

    Atom *atom = lmp->atom;
    Domain *domain = lmp->domain;
    int nlocal = atom->nlocal;

    bigint natoms_prev = atom->natoms;
    int nlocal_prev = nlocal;
    double xdata[3];

    for (int i = 0; i < n; i++) {
      xdata[0] = x[3*i];
      xdata[1] = x[3*i+1];
      xdata[2] = x[3*i+2];
      imageint * img = image ? &image[i] : NULL;
      tagint     tag = id    ? id[i]     : -1;
      if (!domain->ownatom(tag, xdata, img, shrinkexceed)) continue;

      atom->avec->create_atom(type[i],xdata);
      if (id) atom->tag[nlocal] = id[i];
      else atom->tag[nlocal] = i+1;
      if (v) {
        atom->v[nlocal][0] = v[3*i];
        atom->v[nlocal][1] = v[3*i+1];
        atom->v[nlocal][2] = v[3*i+2];
      }
      if (image) atom->image[nlocal] = image[i];
      nlocal++;
    }

    // need to reset atom->natoms inside LAMMPS

    bigint ncurrent = nlocal;
    MPI_Allreduce(&ncurrent,&lmp->atom->natoms,1,MPI_LMP_BIGINT,
                  MPI_SUM,lmp->world);

    // init per-atom fix/compute/variable values for created atoms

    atom->data_fix_compute_variable(nlocal_prev,nlocal);

    // if global map exists, reset it
    // invoke map_init() b/c atom count has grown

    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // warn if new natoms is not correct

    if (lmp->atom->natoms != natoms_prev + n) {
      char str[128];
      snprintf(str, 128, "Library warning in lammps_create_atoms, "
              "invalid total atoms " BIGINT_FORMAT " " BIGINT_FORMAT,
              lmp->atom->natoms,natoms_prev+n);
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,str);
    }
  }
  END_CAPTURE
}

/* ----------------------------------------------------------------------
   find fix external with given ID and set the callback function
   and caller pointer
------------------------------------------------------------------------- */

void lammps_set_fix_external_callback(void *ptr, char *id, FixExternalFnPtr callback_ptr, void * caller)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  FixExternal::FnPtr callback = (FixExternal::FnPtr) callback_ptr;

  BEGIN_CAPTURE
  {
    int ifix = lmp->modify->find_fix(id);
    if (ifix < 0) {
      char str[128];
      snprintf(str, 128, "Can not find fix with ID '%s'!", id);
      lmp->error->all(FLERR,str);
    }

    Fix *fix = lmp->modify->fix[ifix];

    if (strcmp("external",fix->style) != 0){
      char str[128];
      snprintf(str, 128, "Fix '%s' is not of style external!", id);
      lmp->error->all(FLERR,str);
    }

    FixExternal * fext = (FixExternal*) fix;
    fext->set_callback(callback, caller);
  }
  END_CAPTURE
}


// ----------------------------------------------------------------------
// library API functions for accessing LAMMPS configuration
// ----------------------------------------------------------------------

int lammps_config_has_package(char * package_name) {
  return Info::has_package(package_name);
}

int lammps_config_package_count() {
  int i = 0;
  while(LAMMPS::installed_packages[i] != NULL) {
    ++i;
  }
  return i;
}

int lammps_config_package_name(int index, char * buffer, int max_size) {
  int i = 0;
  while(LAMMPS::installed_packages[i] != NULL && i < index) {
    ++i;
  }

  if(LAMMPS::installed_packages[i] != NULL) {
    strncpy(buffer, LAMMPS::installed_packages[i], max_size);
    return true;
  }

  return false;
}

int lammps_has_style(void * ptr, char * category, char * name) {
  LAMMPS *lmp = (LAMMPS *) ptr;
  Info info(lmp);
  return info.has_style(category, name);
}

int lammps_style_count(void * ptr, char * category) {
  LAMMPS *lmp = (LAMMPS *) ptr;
  Info info(lmp);
  return info.get_available_styles(category).size();
}

int lammps_style_name(void* ptr, char * category, int index, char * buffer, int max_size) {
  LAMMPS *lmp = (LAMMPS *) ptr;
  Info info(lmp);
  auto styles = info.get_available_styles(category);

  if (index < styles.size()) {
    strncpy(buffer, styles[index].c_str(), max_size);
    return true;
  }

  return false;
}

int lammps_config_has_gzip_support() {
  return Info::has_gzip_support();
}

int lammps_config_has_png_support() {
  return Info::has_png_support();
}

int lammps_config_has_jpeg_support() {
  return Info::has_jpeg_support();
}

int lammps_config_has_ffmpeg_support() {
  return Info::has_ffmpeg_support();
}

int lammps_config_has_exceptions() {
  return Info::has_exceptions();
}

// ----------------------------------------------------------------------
// library API functions for error handling
// ----------------------------------------------------------------------

#ifdef LAMMPS_EXCEPTIONS

/* ----------------------------------------------------------------------
   check if a new error message
------------------------------------------------------------------------- */

int lammps_has_error(void *ptr) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Error * error = lmp->error;
  return error->get_last_error() ? 1 : 0;
}

/* ----------------------------------------------------------------------
   copy the last error message of LAMMPS into a character buffer
   return value encodes which type of error:
   1 = normal error (recoverable)
   2 = abort error (non-recoverable)
------------------------------------------------------------------------- */

int lammps_get_last_error_message(void *ptr, char * buffer, int buffer_size) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Error * error = lmp->error;

  if(error->get_last_error()) {
    int error_type = error->get_last_error_type();
    strncpy(buffer, error->get_last_error(), buffer_size-1);
    error->set_last_error(NULL, ERROR_NONE);
    return error_type;
  }
  return 0;
}

#endif

/*******************************************************************************
 * Find neighbor list index of pair style neighbor list
 *
 * Try finding pair instance that matches style. If exact is set, the pair must
 * match style exactly. If exact is 0, style must only be contained. If pair is
 * of style pair/hybrid, style is instead matched the nsub-th hybrid sub-style.
 *
 * Once the pair instance has been identified, multiple neighbor list requests
 * may be found. Every neighbor list is uniquely identified by its request
 * index. Thus, providing this request index ensures that the correct neighbor
 * list index is returned.
 *
 * @param ptr      Pointer to LAMMPS instance
 * @param style    String used to search for pair style instance
 * @param exact    Flag to control whether style should match exactly or only
 *                 must be contained in pair style name
 * @param nsub     match nsub-th hybrid sub-style
 * @param request  request index that specifies which neighbor list should be
 *                 returned, in case there are multiple neighbor lists requests
 *                 for the found pair style
 * @return         return neighbor list index if found, otherwise -1
 ******************************************************************************/
int lammps_find_pair_neighlist(void* ptr, char * style, int exact, int nsub, int request) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Pair* pair = lmp->force->pair_match(style, exact, nsub);

  if (pair != NULL) {
    // find neigh list
    for (int i = 0; i < lmp->neighbor->nlist; i++) {
      NeighList * list = lmp->neighbor->lists[i];
      if (list->requestor_type != NeighList::PAIR || pair != list->requestor) continue;

      if (list->index == request) {
          return i;
      }
    }
  }
  return -1;
}

/*******************************************************************************
 * Find neighbor list index of fix neighbor list
 *
 * @param ptr      Pointer to LAMMPS instance
 * @param id       Identifier of fix instance
 * @param request  request index that specifies which request should be returned,
 *                 in case there are multiple neighbor lists for this fix
 * @return         return neighbor list index if found, otherwise -1
 ******************************************************************************/
int lammps_find_fix_neighlist(void* ptr, char * id, int request) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Fix* fix = NULL;
  const int nfix = lmp->modify->nfix;

  // find fix with name
  for (int ifix = 0; ifix < nfix; ifix++) {
    if (strcmp(lmp->modify->fix[ifix]->id, id) == 0) {
        fix = lmp->modify->fix[ifix];
        break;
    }
  }

  if (fix != NULL) {
    // find neigh list
    for (int i = 0; i < lmp->neighbor->nlist; i++) {
      NeighList * list = lmp->neighbor->lists[i];
      if (list->requestor_type != NeighList::FIX || fix != list->requestor) continue;

      if (list->index == request) {
          return i;
      }
    }
  }
  return -1;
}

/*******************************************************************************
 * Find neighbor list index of compute neighbor list
 *
 * @param ptr      Pointer to LAMMPS instance
 * @param id       Identifier of fix instance
 * @param request  request index that specifies which request should be returned,
 *                 in case there are multiple neighbor lists for this fix
 * @return         return neighbor list index if found, otherwise -1
 ******************************************************************************/
int lammps_find_compute_neighlist(void* ptr, char * id, int request) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Compute* compute = NULL;
  const int ncompute = lmp->modify->ncompute;

  // find compute with name
  for (int icompute = 0; icompute < ncompute; icompute++) {
    if (strcmp(lmp->modify->compute[icompute]->id, id) == 0) {
        compute = lmp->modify->compute[icompute];
        break;
    }
  }

  if (compute == NULL) {
    // find neigh list
    for (int i = 0; i < lmp->neighbor->nlist; i++) {
      NeighList * list = lmp->neighbor->lists[i];
      if (list->requestor_type != NeighList::COMPUTE || compute != list->requestor) continue;

      if (list->index == request) {
          return i;
      }
    }
  }
  return -1;
}

/*******************************************************************************
 * Return the number of entries in the neighbor list with given index
 *
 * @param ptr      Pointer to LAMMPS instance
 * @param idx      neighbor list index
 * @return         return number of entries in neighbor list, -1 if idx is
 *                 not a valid index
 ******************************************************************************/
int lammps_neighlist_num_elements(void * ptr, int idx) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Neighbor * neighbor = lmp->neighbor;

  if(idx < 0 || idx >= neighbor->nlist) {
    return -1;
  }

  NeighList * list = neighbor->lists[idx];
  return list->inum;
}

/*******************************************************************************
 * Return atom local index, number of neighbors, and array of neighbor local
 * atom indices of neighbor list entry
 *
 * @param ptr             Pointer to LAMMPS instance
 * @param idx             neighbor list index
 * @param element         neighbor list element index
 * @param[out] iatom      atom local index in range [0, nlocal + nghost), -1 if
                          invalid idx or element index
 * @param[out] numneigh   number of neighbors of atom i or 0
 * @param[out] neighbors  pointer to array of neighbor atom local indices or
 *                        NULL
 ******************************************************************************/
void lammps_neighlist_element_neighbors(void * ptr, int idx, int element, int * iatom, int * numneigh, int ** neighbors) {
  LAMMPS *  lmp = (LAMMPS *) ptr;
  Neighbor * neighbor = lmp->neighbor;
  *iatom = -1;
  *numneigh = 0;
  *neighbors = NULL;

  if(idx < 0 || idx >= neighbor->nlist) {
    return;
  }

  NeighList * list = neighbor->lists[idx];

  if(element < 0 || element >= list->inum) {
    return;
  }

  int i = list->ilist[element];
  *iatom     = i;
  *numneigh  = list->numneigh[i];
  *neighbors = list->firstneigh[i];
}
