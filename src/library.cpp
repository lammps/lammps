// clang-format off
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

// C style library interface to LAMMPS.
// See the manual for detailed documentation.

#define LAMMPS_LIB_MPI 1
#include "library.h"
#include <mpi.h>

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "exceptions.h"
#include "fix.h"
#include "fix_external.h"
#include "force.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "lmppython.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#if defined(LMP_PLUGIN)
#include "plugin.h"
#endif
#include "region.h"
#include "respa.h"
#include "thermo.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cstring>

#if defined(LMP_PYTHON)
#include <Python.h>
#endif

/// string buffer for error messages of global errors
static std::string lammps_last_global_errormessage;

using namespace LAMMPS_NS;

// for printing the non-null pointer argument warning only once

static int ptr_argument_flag = 1;
static void ptr_argument_warning()
{
  if (!ptr_argument_flag) return;
  fprintf(stderr,"Using a 'void **' argument to return the LAMMPS handle "
          "is deprecated.  Please use the return value instead.\n");
  ptr_argument_flag = 0;
}

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

#define BEGIN_CAPTURE \
  Error *error = lmp->error; \
  try

#define END_CAPTURE \
  catch(LAMMPSAbortException &ae) { \
    int nprocs = 0; \
    MPI_Comm_size(ae.get_universe(), &nprocs ); \
    \
    if (nprocs > 1) { \
      error->set_last_error(ae.what(), ERROR_ABORT); \
    } else { \
      error->set_last_error(ae.what(), ERROR_NORMAL); \
    } \
  } catch(LAMMPSException &e) { \
    error->set_last_error(e.what(), ERROR_NORMAL); \
  }

// ----------------------------------------------------------------------
// Library functions to create/destroy an instance of LAMMPS
// ----------------------------------------------------------------------

/** Create instance of the LAMMPS class and return pointer to it.
 *
\verbatim embed:rst

The :cpp:func:`lammps_open` function creates a new :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` class instance while passing in a list of strings
as if they were :doc:`command-line arguments <Run_options>` for the
LAMMPS executable, and an MPI communicator for LAMMPS to run under.
Since the list of arguments is **exactly** as when called from the
command line, the first argument would be the name of the executable and
thus is otherwise ignored.  However ``argc`` may be set to 0 and then
``argv`` may be ``NULL``.  If MPI is not yet initialized, ``MPI_Init()``
will be called during creation of the LAMMPS class instance.

If for some reason the creation or initialization of the LAMMPS instance
fails a null pointer is returned.

.. versionchanged:: 18Sep2020

   This function now has the pointer to the created LAMMPS class
   instance as return value.  For backward compatibility it is still
   possible to provide the address of a pointer variable as final
   argument *ptr*\ .

.. deprecated:: 18Sep2020

   The *ptr* argument will be removed in a future release of LAMMPS.
   It should be set to ``NULL`` instead.

.. note::

   This function is **only** declared when the code using the LAMMPS
   ``library.h`` include file is compiled with ``-DLAMMPS_LIB_MPI``,
   or contains a ``#define LAMMPS_LIB_MPI 1`` statement before
   ``#include "library.h"``.  Otherwise you can only use the
   :cpp:func:`lammps_open_no_mpi` or :cpp:func:`lammps_open_fortran`
   functions.

*See also*
   :cpp:func:`lammps_open_no_mpi`, :cpp:func:`lammps_open_fortran`

\endverbatim
 *
 * \param  argc  number of command line arguments
 * \param  argv  list of command line argument strings
 * \param  comm  MPI communicator for this LAMMPS instance
 * \param  ptr   pointer to a void pointer variable which serves
 *               as a handle; may be ``NULL``
 * \return       pointer to new LAMMPS instance cast to ``void *`` */

void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr)
{
  LAMMPS *lmp = nullptr;
  lammps_mpi_init();
  if (ptr) ptr_argument_warning();

  try {
    lammps_last_global_errormessage.clear();
    lmp = new LAMMPS(argc, argv, comm);
    if (ptr) *ptr = (void *) lmp;
  } catch (fmt::format_error &fe) {
    lammps_last_global_errormessage = fe.what();
    fprintf(stderr, "fmt::format_error: %s\n", fe.what());
    if (ptr) *ptr = nullptr;
  } catch(LAMMPSException &e) {
    lammps_last_global_errormessage = e.what();

    fmt::print(stderr, "LAMMPS Exception: {}", e.what());
    if (ptr) *ptr = nullptr;
  }
  return (void *) lmp;
}

/* ---------------------------------------------------------------------- */

/** Variant of ``lammps_open()`` that implicitly uses ``MPI_COMM_WORLD``.
 *
\verbatim embed:rst

This function is a version of :cpp:func:`lammps_open`, that is missing
the MPI communicator argument.  It will use ``MPI_COMM_WORLD`` instead.
The type and purpose of arguments and return value are otherwise the
same.

Outside of the convenience, this function is useful, when the LAMMPS
library was compiled in serial mode, but the calling code runs in
parallel and the ``MPI_Comm`` data type of the STUBS library would not
be compatible with that of the calling code.

If for some reason the creation or initialization of the LAMMPS instance
fails a null pointer is returned.

.. versionchanged:: 18Sep2020

   This function now has the pointer to the created LAMMPS class
   instance as return value.  For backward compatibility it is still
   possible to provide the address of a pointer variable as final
   argument *ptr*\ .

.. deprecated:: 18Sep2020

   The *ptr* argument will be removed in a future release of LAMMPS.
   It should be set to ``NULL`` instead.


*See also*
   :cpp:func:`lammps_open`, :cpp:func:`lammps_open_fortran`

\endverbatim
 *
 * \param  argc  number of command line arguments
 * \param  argv  list of command line argument strings
 * \param  ptr   pointer to a void pointer variable
 *               which serves as a handle; may be ``NULL``
 * \return       pointer to new LAMMPS instance cast to ``void *`` */

void *lammps_open_no_mpi(int argc, char **argv, void **ptr)
{
  return lammps_open(argc, argv, MPI_COMM_WORLD, ptr);
}

/* ---------------------------------------------------------------------- */

/** Variant of ``lammps_open()`` using a Fortran MPI communicator.
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

This function is a version of :cpp:func:`lammps_open`, that uses an
integer for the MPI communicator as the MPI Fortran interface does.  It
is used in the :f:func:`lammps` constructor of the LAMMPS Fortran
module.  Internally it converts the *f_comm* argument into a C-style MPI
communicator with ``MPI_Comm_f2c()`` and then calls
:cpp:func:`lammps_open`.

If for some reason the creation or initialization of the LAMMPS instance
fails a null pointer is returned.

*See also*
   :cpp:func:`lammps_open_fortran`, :cpp:func:`lammps_open_no_mpi`

\endverbatim
 *
 * \param  argc   number of command line arguments
 * \param  argv   list of command line argument strings
 * \param  f_comm Fortran style MPI communicator for this LAMMPS instance
 * \return        pointer to new LAMMPS instance cast to ``void *`` */

void *lammps_open_fortran(int argc, char **argv, int f_comm)
{
  lammps_mpi_init();
  MPI_Comm c_comm = MPI_Comm_f2c((MPI_Fint)f_comm);
  return lammps_open(argc, argv, c_comm, nullptr);
}

/* ---------------------------------------------------------------------- */

/** Delete a LAMMPS instance created by lammps_open() or its variants.
 *
\verbatim embed:rst

This function deletes the LAMMPS class instance pointed to by ``handle``
that was created by one of the :cpp:func:`lammps_open` variants.  It
does **not** call ``MPI_Finalize()`` to allow creating and deleting
multiple LAMMPS instances concurrently or sequentially.  See
:cpp:func:`lammps_mpi_finalize` for a function performing this operation.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance */

void lammps_close(void *handle)
{
  auto lmp = (LAMMPS *) handle;
  delete lmp;
}

/* ---------------------------------------------------------------------- */

/** Ensure the MPI environment is initialized.
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

The MPI standard requires that any MPI application must call
``MPI_Init()`` exactly once before performing any other MPI function
calls.  This function checks, whether MPI is already initialized and
calls ``MPI_Init()`` in case it is not.

\endverbatim */

void lammps_mpi_init()
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    // provide a dummy argc and argv for MPI_Init().
    int argc = 1;
    char *args[] = { (char *)"liblammps" , nullptr  };
    char **argv = args;
    MPI_Init(&argc,&argv);
  }
}

/* ---------------------------------------------------------------------- */

/** Shut down the MPI infrastructure.
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

The MPI standard requires that any MPI application calls
``MPI_Finalize()`` before exiting.  Even if a calling program does not
do any MPI calls, MPI is still initialized internally to avoid errors
accessing any MPI functions.  This function should then be called right
before exiting the program to wait until all (parallel) tasks are
completed and then MPI is cleanly shut down.  After calling this
function no more MPI calls may be made.

*See also*
   :cpp:func:`lammps_kokkos_finalize`, :cpp:func:`lammps_python_finalize`
\endverbatim */

void lammps_mpi_finalize()
{
  int flag;
  MPI_Initialized(&flag);
  if (flag) {
    MPI_Finalized(&flag);
    if (!flag) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
    }
  }
}

/* ---------------------------------------------------------------------- */

/** Shut down the Kokkos library environment.
 *
\verbatim embed:rst

.. versionadded:: 2Jul2021

The Kokkos library may only be initialized once during the execution of
a process.  This is done automatically the first time Kokkos
functionality is used.  This requires that the Kokkos environment
must be explicitly shut down after any LAMMPS instance using it is
closed (to release associated resources).
After calling this function no Kokkos functionality may be used.

*See also*
   :cpp:func:`lammps_mpi_finalize`, :cpp:func:`lammps_python_finalize`
\endverbatim */

void lammps_kokkos_finalize()
{
  KokkosLMP::finalize();
}

/* ---------------------------------------------------------------------- */

/** Clear the embedded Python environment
 *
\verbatim embed:rst

.. versionadded:: 20Sep2021

This function resets and clears an embedded Python environment
by calling the `Py_Finalize() function
<https://docs.python.org/3/c-api/init.html#c.Py_FinalizeEx>`_
of the embedded Python library, if enabled.
This call would free up all allocated resources and release
loaded shared objects.

However, this is **not** done when a LAMMPS instance is deleted because
a) LAMMPS may have been used through the Python module and thus
the Python interpreter is external and not embedded into LAMMPS
and therefore may not be reset by LAMMPS b) some Python modules
and extensions, most notably NumPy, are not compatible with being
initialized multiple times, which would happen if additional
LAMMPS instances using Python would be created *after*
after calling Py_Finalize().

This function can be called to explicitly clear the Python
environment in case it is safe to do so.

*See also*
   :cpp:func:`lammps_mpi_finalize`, :cpp:func:`lammps_kokkos_finalize`
\endverbatim */

void lammps_python_finalize()
{
  Python::finalize();
}


/* ---------------------------------------------------------------------- */

/** Call a LAMMPS Error class function
 *
\verbatim embed:rst

.. versionadded:: 3Nov2022

This function is a wrapper around functions in the ``Error`` to print an
error message and then stop LAMMPS.

The *error_type* parameter selects which function to call.  It is a sum
of constants from :cpp:enum:`_LMP_ERROR_CONST`.  If the value does not
match any valid combination of constants a warning is printed and the
function returns.

\endverbatim
 *
 * \param  handle       pointer to a previously created LAMMPS instance
 * \param  error_type   parameter to select function in the Error class
 * \param  error_text   error message */

void lammps_error(void *handle, int error_type, const char *error_text)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    switch (error_type) {
    case LMP_ERROR_WARNING:
      lmp->error->warning("(library)", 0, error_text);
      break;
    case LMP_ERROR_ONE:
      lmp->error->one("(library)", 0, error_text);
      break;
    case LMP_ERROR_ALL:
      lmp->error->all("(library)", 0, error_text);
      break;
    case LMP_ERROR_WARNING|LMP_ERROR_WORLD:
      lmp->error->warning("(library)", 0, error_text);
      break;
    case LMP_ERROR_ONE|LMP_ERROR_WORLD:
      lmp->error->one("(library)", 0, error_text);
      break;
    case LMP_ERROR_ALL|LMP_ERROR_WORLD:
      lmp->error->all("(library)", 0, error_text);
      break;
    case LMP_ERROR_WARNING|LMP_ERROR_UNIVERSE:
      lmp->error->universe_warn("(library)", 0, error_text);
      break;
    case LMP_ERROR_ONE|LMP_ERROR_UNIVERSE:
      lmp->error->universe_one("(library)", 0, error_text);
      break;
    case LMP_ERROR_ALL|LMP_ERROR_UNIVERSE:
      lmp->error->universe_all("(library)", 0, error_text);
      break;
    default:
      auto mesg = fmt::format("Unknown error type {} for message: {}", error_type, error_text);
      lmp->error->warning("(library)", 0, mesg);
    }
  }
  END_CAPTURE

    // with enabled exceptions the above code will simply throw an
    // exception and record the error message. So we have to explicitly
    // stop here like we do in main.cpp
  if (lammps_has_error(handle)) {
    if (error_type & 1) {
      lammps_kokkos_finalize();
      lammps_python_finalize();
      MPI_Abort(lmp->universe->uworld, 1);
    } else if (error_type & 2) {
      lammps_kokkos_finalize();
      lammps_python_finalize();
      lammps_mpi_finalize();
      exit(1);
    }
  }
}

// ----------------------------------------------------------------------
// Library functions to process commands
// ----------------------------------------------------------------------

/** Process LAMMPS input from a file.
 *
\verbatim embed:rst

This function processes commands in the file pointed to by *filename*
line by line and thus functions very similar to the :doc:`include
<include>` command. The function returns when the end of the file is
reached and the commands have completed.

The actual work is done by the functions
:cpp:func:`Input::file(const char *)<void LAMMPS_NS::Input::file(const char *)>`
and :cpp:func:`Input::file()<void LAMMPS_NS::Input::file()>`.

\endverbatim
 *
 * \param  handle    pointer to a previously created LAMMPS instance
 * \param  filename  name of a file with LAMMPS input */

void lammps_file(void *handle, const char *filename)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0)
      lmp->error->all(FLERR, "Library error: issuing LAMMPS commands during a run is not allowed");
    else
      lmp->input->file(filename);
  }
  END_CAPTURE
}

/* ---------------------------------------------------------------------- */

/** Process a single LAMMPS input command from a string.
 *
\verbatim embed:rst

This function tells LAMMPS to execute the single command in the string
*cmd*.  The entire string is considered as command and need not have a
(final) newline character.  Newline characters in the body of the
string, however, will be treated as part of the command and will **not**
start a second command.  The function :cpp:func:`lammps_commands_string`
processes a string with multiple command lines.

The function returns the name of the command on success or ``NULL`` when
passing a string without a command.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  cmd     string with a single LAMMPS command
 * \return         string with parsed command name or ``NULL`` */

char *lammps_command(void *handle, const char *cmd)
{
  auto lmp = (LAMMPS *) handle;
  char *result = nullptr;

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0)
      lmp->error->all(FLERR,"Library error: issuing LAMMPS commands "
                      "during a run is not allowed.");
    else
      result = lmp->input->one(cmd);
  }
  END_CAPTURE

  return result;
}

/* ---------------------------------------------------------------------- */

/** Process multiple LAMMPS input commands from list of strings.
 *
\verbatim embed:rst

This function processes multiple commands from a list of strings by
first concatenating the individual strings in *cmds* into a single
string, inserting newline characters as needed.  The combined string
is passed to :cpp:func:`lammps_commands_string` for processing.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  ncmd    number of lines in *cmds*
 * \param  cmds    list of strings with LAMMPS commands */

void lammps_commands_list(void *handle, int ncmd, const char **cmds)
{
  std::string allcmds;

  for (int i = 0; i < ncmd; i++) {
    allcmds.append(cmds[i]);
    if (allcmds.empty() || (allcmds.back() != '\n')) allcmds.append(1,'\n');
  }

  lammps_commands_string(handle,allcmds.c_str());
}

/* ---------------------------------------------------------------------- */

/** Process a block of LAMMPS input commands from a single string.
 *
\verbatim embed:rst

This function processes a multi-line string similar to a block of
commands from a file.  The string may have multiple lines (separated by
newline characters) and also single commands may be distributed over
multiple lines with continuation characters ('&').  Those lines are
combined by removing the '&' and the following newline character.  After
this processing the string is handed to LAMMPS for parsing and
executing.

.. note::

   Multi-line commands enabled by triple quotes will NOT work with
   this function.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  str     string with block of LAMMPS input commands */

void lammps_commands_string(void *handle, const char *str)
{
  auto lmp = (LAMMPS *) handle;

  // copy str and convert from CR-LF (DOS-style) to LF (Unix style) line
  int n = strlen(str);
  char *ptr, *copy = new char[n+1];

  for (ptr = copy; *str != '\0'; ++str) {
    if ((str[0] == '\r') && (str[1] == '\n')) continue;
    *ptr++ = *str;
  }
  *ptr = '\0';

  BEGIN_CAPTURE
  {
    if (lmp->update->whichflag != 0) {
      lmp->error->all(FLERR,"Library error: issuing LAMMPS command during run");
    }

    n = strlen(copy);
    ptr = copy;
    for (int i=0; i < n; ++i) {

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

  delete[] copy;
}

// -----------------------------------------------------------------------
// Library functions to extract info from LAMMPS or set data in LAMMPS
// -----------------------------------------------------------------------

/** Return the total number of atoms in the system.
 *
\verbatim embed:rst

This number may be very large when running large simulations across
multiple processes.  Depending on compile time choices, LAMMPS may be
using either 32-bit or a 64-bit integer to store this number. For
portability this function returns thus a double precision
floating point number, which can represent up to a 53-bit signed
integer exactly (:math:`\approx 10^{16}`).

As an alternative, you can use :cpp:func:`lammps_extract_global`
and cast the resulting pointer to an integer pointer of the correct
size and dereference it.  The size of that integer (in bytes) can be
queried by calling :cpp:func:`lammps_extract_setting` to return
the size of a ``bigint`` integer.

.. versionchanged:: 18Sep2020

   The type of the return value was changed from ``int`` to ``double``
   to accommodate reporting atom counts for larger systems that would
   overflow a 32-bit int without having to depend on a 64-bit bit
   integer type definition.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \return         total number of atoms or 0 if value is too large */

double lammps_get_natoms(void *handle)
{
  auto lmp = (LAMMPS *) handle;

  auto  natoms = static_cast<double>(lmp->atom->natoms);
  if (natoms > 9.0e15) return 0; // TODO:XXX why not -1?
  return natoms;
}

/* ---------------------------------------------------------------------- */

/** Evaluate a thermo keyword.
 *
\verbatim embed:rst

This function returns the current value of a :doc:`thermo keyword <thermo_style>`.
Unlike :cpp:func:`lammps_extract_global` it does not give access to the
storage of the desired data but returns its value as a ``double``, so it
can also return information that is computed on-the-fly.
Use :cpp:func:`lammps_last_thermo` to get access to the cached data from
the last thermo output.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  keyword  string with the name of the thermo keyword
 * \return          value of the requested thermo property or 0.0 */

double lammps_get_thermo(void *handle, const char *keyword)
{
  auto lmp = (LAMMPS *) handle;
  double dval = 0.0;

  BEGIN_CAPTURE
  {
    lmp->output->thermo->evaluate_keyword(keyword,&dval);
  }
  END_CAPTURE

  return dval;
}

/* ---------------------------------------------------------------------- */

/** Access cached data from last thermo output
 *
\verbatim embed:rst

.. versionadded:: 15Jun2023

This function provides access to cached data from the last thermo output.
This differs from :cpp:func:`lammps_get_thermo` in that it does not trigger
an evaluation.  Instead it provides direct access to a read-only location
of the last thermo output data and the corresponding keyword strings.
The how to handle the return value depends on the value of the *what*
argument string.

.. note::

   The *type* property points to a static location that is reassigned
   with every call, so the returned pointer should be recast,
   dereferenced, and assigned immediately. Otherwise, its value may be
   changed with the next invocation of the function.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Value of *what*
     - Description of return value
     - Data type
     - Uses index
   * - step
     - timestep when the last thermo output was generated or -1
     - pointer to bigint
     - no
   * - num
     - number of fields in thermo output
     - pointer to int
     - no
   * - keyword
     - column keyword for thermo output
     - pointer to 0-terminated const char array
     - yes
   * - type
     - data type of thermo output column; see :cpp:enum:`_LMP_DATATYPE_CONST`
     - pointer to int
     - yes
   * - data
     - actual field data for column
     - pointer to int, int64_t or double
     - yes

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  what     string with the kind of data requested
 * \param  index    integer with index into data arrays, ignored for scalar data
 * \return          pointer to location of requested data cast to void or NULL */

void *lammps_last_thermo(void *handle, const char *what, int index)
{
  auto lmp = (LAMMPS *) handle;
  void *val = nullptr;
  Thermo *th = lmp->output->thermo;
  if (!th) return nullptr;
  const int nfield = *th->get_nfield();

  BEGIN_CAPTURE
  {
    if (strcmp(what, "step") == 0) {
      val = (void *) th->get_timestep();

    } else if (strcmp(what, "num") == 0) {
      val = (void *) th->get_nfield();

    } else if (strcmp(what, "keyword") == 0) {
      if ((index < 0) || (index >= nfield)) return nullptr;
      const auto &keywords = th->get_keywords();
      val = (void *) keywords[index].c_str();

    } else if (strcmp(what, "type") == 0) {
      if ((index < 0) || (index >= nfield)) return nullptr;
      const auto &field = th->get_fields()[index];
      val = (void *) &field.type;
    } else if (strcmp(what, "data") == 0) {
      if ((index < 0) || (index >= nfield)) return nullptr;
      const auto &field = th->get_fields()[index];
      if (field.type == multitype::LAMMPS_INT) {
        val = (void *) &field.data.i;
      } else if (field.type == multitype::LAMMPS_INT64) {
        val = (void *) &field.data.b;
      } else if (field.type == multitype::LAMMPS_DOUBLE) {
        val = (void *) &field.data.d;
      }

    } else val = nullptr;
  }
  END_CAPTURE
  return val;
}

/* ---------------------------------------------------------------------- */

/** Extract simulation box parameters.
 *
\verbatim embed:rst

This function (re-)initializes the simulation box and boundary
information and then assign the designated data to the locations in the
pointers passed as arguments. Any argument (except the first) may be
a NULL pointer and then will not be assigned.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  boxlo    pointer to 3 doubles where the lower box boundary is stored
 * \param  boxhi    pointer to 3 doubles where the upper box boundary is stored
 * \param  xy       pointer to a double where the xy tilt factor is stored
 * \param  yz       pointer to a double where the yz tilt factor is stored
 * \param  xz       pointer to a double where the xz tilt factor is stored
 * \param  pflags   pointer to 3 ints, set to 1 for periodic boundaries
                    and 0 for non-periodic
 * \param  boxflag  pointer to an int, which is set to 1 if the box will be
 *                  changed during a simulation by a fix and 0 if not. */

void lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                        double *xy, double *yz, double *xz,
                        int *pflags, int *boxflag)
{
  auto lmp = (LAMMPS *) handle;
  Domain *domain = lmp->domain;

  BEGIN_CAPTURE
  {
    // do nothing if box does not yet exist
    if ((lmp->domain->box_exist == 0)
        && (lmp->comm->me == 0)) {
      lmp->error->warning(FLERR,"Calling lammps_extract_box without a box");
      return;
    }

    // domain->init() is needed to update domain->box_change
    domain->init();

    if (boxlo) {
      boxlo[0] = domain->boxlo[0];
      boxlo[1] = domain->boxlo[1];
      boxlo[2] = domain->boxlo[2];
    }
    if (boxhi) {
      boxhi[0] = domain->boxhi[0];
      boxhi[1] = domain->boxhi[1];
      boxhi[2] = domain->boxhi[2];
    }
    if (xy) *xy = domain->xy;
    if (yz) *yz = domain->yz;
    if (xz) *xz = domain->xz;

    if (pflags) {
      pflags[0] = domain->periodicity[0];
      pflags[1] = domain->periodicity[1];
      pflags[2] = domain->periodicity[2];
    }
    if (boxflag) *boxflag = domain->box_change;
  }
  END_CAPTURE
}

/* ---------------------------------------------------------------------- */

/** Reset simulation box parameters.
 *
\verbatim embed:rst

This function sets the simulation box dimensions (upper and lower bounds
and tilt factors) from the provided data and then re-initializes the box
information and all derived settings. It may only be called before atoms
are created.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  boxlo    pointer to 3 doubles containing the lower box boundary
 * \param  boxhi    pointer to 3 doubles containing the upper box boundary
 * \param  xy       xy tilt factor
 * \param  yz       yz tilt factor
 * \param  xz       xz tilt factor */

void lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                      double xy, double yz, double xz)
{
  auto lmp = (LAMMPS *) handle;
  Domain *domain = lmp->domain;

  BEGIN_CAPTURE
  {
    if (lmp->atom->natoms > 0)
      lmp->error->all(FLERR,"Calling lammps_reset_box not supported when atoms exist");

    // warn and do nothing if no box exists
    if (lmp->domain->box_exist == 0) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Ignoring call to lammps_reset_box without a box");
      return;
    }

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
  END_CAPTURE
}

/* ---------------------------------------------------------------------- */

/** Get memory usage information
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

This function will retrieve memory usage information for the current
LAMMPS instance or process.  The *meminfo* buffer will be filled with
3 different numbers (if supported by the operating system).  The first
is the tally (in MBytes) of all large memory allocations made by LAMMPS.
This is a lower boundary of how much memory is requested and does not
account for memory allocated on the stack or allocations via ``new``.
The second number is the current memory allocation of the current process
as returned by a memory allocation reporting in the system library.  The
third number is the maximum amount of RAM (not swap) used by the process
so far. If any of the two latter parameters is not supported by the operating
system it will be set to zero.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  meminfo  buffer with space for at least 3 double to store
 * data in. */

void lammps_memory_usage(void *handle, double *meminfo)
{
  auto lmp = (LAMMPS *) handle;
  Info info(lmp);
  info.get_memory_info(meminfo);
}

/* ---------------------------------------------------------------------- */

/** Return current LAMMPS world communicator as integer
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

This will take the LAMMPS "world" communicator and convert it to an
integer using ``MPI_Comm_c2f()``, so it is equivalent to the
corresponding MPI communicator in Fortran. This way it can be safely
passed around between different programming languages.  To convert it
to the C language representation use ``MPI_Comm_f2c()``.

If LAMMPS was compiled with MPI_STUBS, this function returns -1.

*See also*
   :cpp:func:`lammps_open_fortran`

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \return         Fortran representation of the LAMMPS world communicator */

int lammps_get_mpi_comm(void *handle)
{
#ifdef MPI_STUBS
  return -1;
#else
  LAMMPS *lmp = (LAMMPS *) handle;
  MPI_Fint f_comm = MPI_Comm_c2f(lmp->world);
  return f_comm;
#endif
}

/* ---------------------------------------------------------------------- */

/** Query LAMMPS about global settings.
 *
\verbatim embed:rst

This function will retrieve or compute global properties. In contrast to
:cpp:func:`lammps_get_thermo` this function returns an ``int``.  The
following tables list the currently supported keyword.  If a keyword is
not recognized, the function returns -1.  The integer sizes functions may
be called without a valid LAMMPS object handle (it is ignored).

* :ref:`Integer sizes <extract_integer_sizes>`
* :ref:`Image masks <extract_image_masks>`
* :ref:`System status <extract_system_status>`
* :ref:`System sizes <extract_system_sizes>`
* :ref:`Atom style flags <extract_atom_flags>`

.. _extract_integer_sizes:

**Integer sizes**

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

.. _extract_image_masks:

**Image masks**

These settings are related to how LAMMPS stores and interprets periodic images. The values are used
internally by the :doc:`Fortran interface <Fortran>` and are not likely to be useful to users.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Description / Return value
   * - IMGMASK
     - Bit-mask used to convert image flags to a single integer
   * - IMGMAX
     - Maximum allowed image number for a particular atom
   * - IMGBITS
     - Bits used in image counts
   * - IMG2BITS
     - Second bitmask used in image counts

.. _extract_system_status:

**System status**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Description / Return value
   * - dimension
     - Number of dimensions: 2 or 3. See :doc:`dimension`.
   * - box_exist
     - 1 if the simulation box is defined, 0 if not.
       See :doc:`create_box`.
   * - kokkos_active
     - 1 if the KOKKOS package is compiled in **and** activated, 0 if not.
       See :doc:`Speed_kokkos`.
   * - kokkos_nthreads
     - Number of Kokkos threads per MPI process, 0 if Kokkos is not active.
       See :doc:`Speed_kokkos`.
   * - kokkos_ngpus
     - Number of Kokkos gpus per physical node, 0 if Kokkos is not active or no GPU support.
       See :doc:`Speed_kokkos`.
   * - nthreads
     - Number of requested OpenMP threads per MPI process for LAMMPS' execution
   * - newton_bond
     - 1 if Newton's 3rd law is applied to bonded interactions, 0 if not.
   * - newton_pair
     - 1 if Newton's 3rd law is applied to non-bonded interactions, 0 if not.
   * - triclinic
     - 1 if the the simulation box is triclinic, 0 if orthogonal.
       See :doc:`change_box`.
   * - universe_rank
     - MPI rank on LAMMPS' universe communicator (0 <= universe_rank < universe_size)
   * - universe_size
     - Number of ranks on LAMMPS' universe communicator (world_size <= universe_size)
   * - world_rank
     - MPI rank on LAMMPS' world communicator (0 <= world_rank < world_size)
   * - world_size
     - Number of ranks on LAMMPS' world communicator

.. _extract_system_sizes:

**System sizes**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Description / Return value
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
   * - nbondtypes
     - number of bond types
   * - nangletypes
     - number of angle types
   * - ndihedraltypes
     - number of dihedral types
   * - nimpropertypes
     - number of improper types
   * - nellipsoids
     - number of atoms that have ellipsoid data
   * - nlines
     - number of atoms that have line data (see :doc:`pair style line/lj <pair_line_lj>`)
   * - ntris
     - number of atoms that have triangle data (see :doc:`pair style tri/lj <pair_tri_lj>`)
   * - nbodies
     - number of atoms that have body data (see :doc:`the Body particle HowTo <Howto_body>`)

.. _extract_atom_flags:

**Atom style flags**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Keyword
     - Description / Return value
   * - molecule_flag
     - 1 if the atom style includes molecular topology data. See :doc:`atom_style`.
   * - q_flag
     - 1 if the atom style includes point charges. See :doc:`atom_style`.
   * - mu_flag
     - 1 if the atom style includes point dipoles. See :doc:`atom_style`.
   * - rmass_flag
     - 1 if the atom style includes per-atom masses, 0 if there are per-type masses. See :doc:`atom_style`.
   * - radius_flag
     - 1 if the atom style includes a per-atom radius. See :doc:`atom_style`.
   * - sphere_flag
     - 1 if the atom style describes extended particles that can rotate. See :doc:`atom_style`.
   * - ellipsoid_flag
     - 1 if the atom style describes extended particles that may be ellipsoidal. See :doc:`atom_style`.
   * - omega_flag
     - 1 if the atom style can store per-atom rotational velocities. See :doc:`atom_style`.
   * - torque_flag
     - 1 if the atom style can store per-atom torques. See :doc:`atom_style`.
   * - angmom_flag
     - 1 if the atom style can store per-atom angular momentum. See :doc:`atom_style`.

*See also*
   :cpp:func:`lammps_extract_global`

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  keyword  string with the name of the thermo keyword
 * \return          value of the queried setting or -1 if unknown */

int lammps_extract_setting(void *handle, const char *keyword)
{
  auto lmp = (LAMMPS *) handle;

// This can be customized by adding keywords and documenting them in the section above.
  if (strcmp(keyword,"bigint") == 0) return sizeof(bigint);
  if (strcmp(keyword,"tagint") == 0) return sizeof(tagint);
  if (strcmp(keyword,"imageint") == 0) return sizeof(imageint);

  if (strcmp(keyword,"IMGMASK") == 0) return IMGMASK;
  if (strcmp(keyword,"IMGBITS") == 0) return IMGBITS;
  if (strcmp(keyword,"IMG2BITS") == 0) return IMG2BITS;
  if (strcmp(keyword,"IMGMAX") == 0) return IMGMAX;

  if (strcmp(keyword,"dimension") == 0) return lmp->domain->dimension;
  if (strcmp(keyword,"box_exist") == 0) return lmp->domain->box_exist;
  if (strcmp(keyword,"kokkos_active") == 0) return (lmp->kokkos) ? 1 : 0;
  if (strcmp(keyword,"kokkos_nthreads") == 0) return (lmp->kokkos) ? lmp->kokkos->nthreads : 0;
  if (strcmp(keyword,"kokkos_ngpus") == 0) return (lmp->kokkos) ? lmp->kokkos->ngpus : 0;
  if (strcmp(keyword,"newton_bond") == 0) return lmp->force->newton_bond;
  if (strcmp(keyword,"newton_pair") == 0) return lmp->force->newton_pair;
  if (strcmp(keyword,"triclinic") == 0) return lmp->domain->triclinic;

  if (strcmp(keyword,"universe_rank") == 0) return lmp->universe->me;
  if (strcmp(keyword,"universe_size") == 0) return lmp->universe->nprocs;
  if (strcmp(keyword,"world_rank") == 0) return lmp->comm->me;
  if (strcmp(keyword,"world_size") == 0) return lmp->comm->nprocs;
  if (strcmp(keyword,"nthreads") == 0) return lmp->comm->nthreads;

  if (strcmp(keyword,"nlocal") == 0) return lmp->atom->nlocal;
  if (strcmp(keyword,"nghost") == 0) return lmp->atom->nghost;
  if (strcmp(keyword,"nall") == 0) return lmp->atom->nlocal+lmp->atom->nghost;
  if (strcmp(keyword,"nmax") == 0) return lmp->atom->nmax;
  if (strcmp(keyword,"ntypes") == 0) return lmp->atom->ntypes;
  if (strcmp(keyword,"nbondtypes") == 0) return lmp->atom->nbondtypes;
  if (strcmp(keyword,"nangletypes") == 0) return lmp->atom->nangletypes;
  if (strcmp(keyword,"ndihedraltypes") == 0) return lmp->atom->ndihedraltypes;
  if (strcmp(keyword,"nimpropertypes") == 0) return lmp->atom->nimpropertypes;
  if (strcmp(keyword,"nellipsoids") == 0) return lmp->atom->nellipsoids;
  if (strcmp(keyword,"nlines") == 0) return lmp->atom->nlines;
  if (strcmp(keyword,"ntris") == 0) return lmp->atom->ntris;
  if (strcmp(keyword,"nbodies") == 0) return lmp->atom->nbodies;

  if (strcmp(keyword,"molecule_flag") == 0) return lmp->atom->molecule_flag;
  if (strcmp(keyword,"q_flag") == 0) return lmp->atom->q_flag;
  if (strcmp(keyword,"mu_flag") == 0) return lmp->atom->mu_flag;
  if (strcmp(keyword,"rmass_flag") == 0) return lmp->atom->rmass_flag;
  if (strcmp(keyword,"radius_flag") == 0) return lmp->atom->radius_flag;
  if (strcmp(keyword,"sphere_flag") == 0) return lmp->atom->sphere_flag;
  if (strcmp(keyword,"ellipsoid_flag") == 0) return lmp->atom->ellipsoid_flag;
  if (strcmp(keyword,"omega_flag") == 0) return lmp->atom->omega_flag;
  if (strcmp(keyword,"torque_flag") == 0) return lmp->atom->torque_flag;
  if (strcmp(keyword,"angmom_flag") == 0) return lmp->atom->angmom_flag;
  if (strcmp(keyword,"peri_flag") == 0) return lmp->atom->peri_flag;

  return -1;
}

/* ---------------------------------------------------------------------- */

/** Get data type of internal global LAMMPS variables or arrays.
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

This function returns an integer that encodes the data type of the global
property with the specified name. See :cpp:enum:`_LMP_DATATYPE_CONST` for valid
values. Callers of :cpp:func:`lammps_extract_global` can use this information
to then decide how to cast the ``void *`` pointer and access the data.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance (unused)
 * \param  name     string with the name of the extracted property
 * \return          integer constant encoding the data type of the property
 *                  or -1 if not found. */

int lammps_extract_global_datatype(void * /*handle*/, const char *name)
{
  if (strcmp(name,"dt") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ntimestep") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"atime") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"atimestep") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"respa_levels") == 0) return LAMMPS_INT;
  if (strcmp(name,"respa_dt") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"boxlo") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxhi") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"sublo") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"subhi") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"sublo_lambda") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"subhi_lambda") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxxlo") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxxhi") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxylo") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxyhi") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxzlo") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"boxzhi") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"periodicity") == 0) return LAMMPS_INT;
  if (strcmp(name,"triclinic") == 0) return LAMMPS_INT;
  if (strcmp(name,"xy") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"xz") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"yz") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"natoms") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"nbonds") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"nangles") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"ndihedrals") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"nimpropers") == 0) return LAMMPS_BIGINT;
  if (strcmp(name,"nlocal") == 0) return LAMMPS_INT;
  if (strcmp(name,"nghost") == 0) return LAMMPS_INT;
  if (strcmp(name,"nmax") == 0) return LAMMPS_INT;
  if (strcmp(name,"ntypes") == 0) return LAMMPS_INT;
  if (strcmp(name,"special_lj") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"special_coul") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"q_flag") == 0) return LAMMPS_INT;

  if (strcmp(name,"units") == 0) return LAMMPS_STRING;
  if (strcmp(name,"atom_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"pair_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"bond_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"angle_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"dihedral_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"improper_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"kspace_style") == 0) return LAMMPS_STRING;
  if (strcmp(name,"boltz") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"hplanck") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"mvv2e") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ftm2v") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"mv2d") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"nktv2p") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"qqr2e") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"qe2f") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"vxmu2f") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"xxt2kmu") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"dielectric") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"qqrd2e") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"e_mass") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"hhmrr2e") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"mvh2r") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"angstrom") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"femtosecond") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"qelectron") == 0) return LAMMPS_DOUBLE;

  return -1;
}

/* ---------------------------------------------------------------------- */

/** Get pointer to internal global LAMMPS variables or arrays.
 *
\verbatim embed:rst

This function returns a pointer to the location of some global property
stored in one of the constituent classes of a LAMMPS instance.  The
returned pointer is cast to ``void *`` and needs to be cast to a pointer
of the type that the entity represents. The pointers returned by this
function are generally persistent; therefore it is not necessary to call
the function again, unless a :doc:`clear` command is issued which wipes
out and recreates the contents of the :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` class.

Please also see :cpp:func:`lammps_extract_setting`,
:cpp:func:`lammps_get_thermo`, and :cpp:func:`lammps_extract_box`.

.. warning::

   Modifying the data in the location pointed to by the returned pointer
   may lead to inconsistent internal data and thus may cause failures or
   crashes or bogus simulations.  In general it is thus usually better
   to use a LAMMPS input command that sets or changes these parameters.
   Those will take care of all side effects and necessary updates of
   settings derived from such settings.  Where possible, a reference to
   such a command or a relevant section of the manual is given below.

The following tables list the supported names, their data types, length
of the data area, and a short description.  The data type can also be
queried through calling :cpp:func:`lammps_extract_global_datatype`.
The ``bigint`` type may be defined to be either an ``int`` or an
``int64_t``.  This is set at :ref:`compile time <size>` of the LAMMPS
library and can be queried through calling
:cpp:func:`lammps_extract_setting`.
The function :cpp:func:`lammps_extract_global_datatype` will directly
report the "native" data type.  The following tables are provided:

* :ref:`Timestep settings <extract_timestep_settings>`
* :ref:`Simulation box settings <extract_box_settings>`
* :ref:`System property settings <extract_system_settings>`
* :ref:`Unit settings <extract_unit_settings>`

.. _extract_timestep_settings:

**Timestep settings**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Length
     - Description
   * - dt
     - double
     - 1
     - length of the time step. See :doc:`timestep`.
   * - ntimestep
     - bigint
     - 1
     - current time step number. See :doc:`reset_timestep`.
   * - atime
     - double
     - 1
     - accumulated simulation time in time units.
   * - atimestep
     - bigint
     - 1
     - the number of the timestep when "atime" was last updated.
   * - respa_levels
     - int
     - 1
     - number of r-RESPA levels. See :doc:`run_style`.
   * - respa_dt
     - double
     - number of r-RESPA levels
     - length of the time steps with r-RESPA. See :doc:`run_style`.

.. _extract_box_settings:

**Simulation box settings**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Length
     - Description
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
   * - sublo
     - double
     - 3
     - subbox lower boundaries
   * - subhi
     - double
     - 3
     - subbox upper boundaries
   * - sublo_lambda
     - double
     - 3
     - subbox lower boundaries in fractional coordinates (for triclinic cells)
   * - subhi_lambda
     - double
     - 3
     - subbox upper boundaries in fractional coordinates (for triclinic cells)
   * - periodicity
     - int
     - 3
     - 0 if non-periodic, 1 if periodic for x, y, and z;
       See :doc:`boundary`.
   * - triclinic
     - int
     - 1
     - 1 if the the simulation box is triclinic, 0 if orthogonal;
       See :doc:`change_box`.
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

.. _extract_system_settings:

**System property settings**

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Length
     - Description
   * - ntypes
     - int
     - 1
     - number of atom types
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
   * - natoms
     - bigint
     - 1
     - total number of atoms in the simulation.
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
   * - special_lj
     - double
     - 4
     - special :doc:`pair weighting factors <special_bonds>` for LJ interactions (first element is always 1.0)
   * - special_coul
     - double
     - 4
     - special :doc:`pair weighting factors <special_bonds>` for Coulomb interactions (first element is always 1.0)
   * - q_flag
     - int
     - 1
     - **deprecated**. Use :cpp:func:`lammps_extract_setting` instead.
   * - atom_style
     - char \*
     - 1
     - string with the current atom style.
   * - pair_style
     - char \*
     - 1
     - string with the current pair style.
   * - bond_style
     - char \*
     - 1
     - string with the current bond style.
   * - angle_style
     - char \*
     - 1
     - string with the current angle style.
   * - dihedral_style
     - char \*
     - 1
     - string with the current dihedral style.
   * - improper_style
     - char \*
     - 1
     - string with the current improper style.
   * - kspace_style
     - char \*
     - 1
     - string with the current KSpace style.

.. _extract_unit_settings:

**Unit settings**

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
     - factor to convert :math:`\frac{1}{2}mv^2` for a particle to
       the current energy unit; See :doc:`units`.
   * - ftm2v
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - mv2d
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - nktv2p
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - qqr2e
     - double
     - 1
     - factor to convert :math:`\frac{q_i q_j}{r}` to energy units;
       See :doc:`units`.
   * - qe2f
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - vxmu2f
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - xxt2kmu
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - dielectric
     - double
     - 1
     - value of the dielectric constant. See :doc:`dielectric`.
   * - qqrd2e
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - e_mass
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - hhmrr2e
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - mvh2r
     - double
     - 1
     - (description missing) See :doc:`units`.
   * - angstrom
     - double
     - 1
     - constant to convert current length unit to angstroms;
       1.0 for reduced (aka "lj") units. See :doc:`units`.
   * - femtosecond
     - double
     - 1
     - constant to convert current time unit to femtoseconds;
       1.0 for reduced (aka "lj") units
   * - qelectron
     - double
     - 1
     - (description missing) See :doc:`units`.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  name     string with the name of the extracted property
 * \return          pointer (cast to ``void *``) to the location of the
                    requested property. NULL if name is not known. */

void *lammps_extract_global(void *handle, const char *name)
{
  auto lmp = (LAMMPS *) handle;

  if (strcmp(name,"units") == 0) return (void *) lmp->update->unit_style;
  if (strcmp(name,"atom_style") == 0) return (void *) lmp->atom->atom_style;
  if (strcmp(name,"pair_style") == 0) return (void *) lmp->force->pair_style;
  if (strcmp(name,"bond_style") == 0) return (void *) lmp->force->bond_style;
  if (strcmp(name,"angle_style") == 0) return (void *) lmp->force->angle_style;
  if (strcmp(name,"dihedral_style") == 0) return (void *) lmp->force->dihedral_style;
  if (strcmp(name,"improper_style") == 0) return (void *) lmp->force->improper_style;
  if (strcmp(name,"kspace_style") == 0) return (void *) lmp->force->kspace_style;
  if (strcmp(name,"dt") == 0) return (void *) &lmp->update->dt;
  if (strcmp(name,"ntimestep") == 0) return (void *) &lmp->update->ntimestep;
  // update->atime can be referenced as a pointer
  // thermo "timer" data cannot be, since it is computed on request
  // lammps_get_thermo() can access all thermo keywords by value
  if (strcmp(name,"atime") == 0) return (void *) &lmp->update->atime;
  if (strcmp(name,"atimestep") == 0) return (void *) &lmp->update->atimestep;

  if (utils::strmatch(lmp->update->integrate_style,"^respa")) {
    auto respa = dynamic_cast<Respa *>(lmp->update->integrate);
    if (strcmp(name,"respa_levels") == 0) return (void *) &respa->nlevels;
    if (strcmp(name,"respa_dt") == 0) return (void *) respa->step;
  }
  if (strcmp(name,"boxlo") == 0) return (void *) lmp->domain->boxlo;
  if (strcmp(name,"boxhi") == 0) return (void *) lmp->domain->boxhi;
  if (strcmp(name,"sublo") == 0) return (void *) lmp->domain->sublo;
  if (strcmp(name,"subhi") == 0) return (void *) lmp->domain->subhi;
  // these are only valid for a triclinic cell
  if (lmp->domain->triclinic) {
    if (strcmp(name,"sublo_lambda") == 0)
      return (void *) lmp->domain->sublo_lamda;
    if (strcmp(name,"subhi_lambda") == 0)
      return (void *) lmp->domain->subhi_lamda;
  }
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
  if (strcmp(name,"ntypes") == 0) return (void *) &lmp->atom->ntypes;
  if (strcmp(name,"nbonds") == 0) return (void *) &lmp->atom->nbonds;
  if (strcmp(name,"nangles") == 0) return (void *) &lmp->atom->nangles;
  if (strcmp(name,"ndihedrals") == 0) return (void *) &lmp->atom->ndihedrals;
  if (strcmp(name,"nimpropers") == 0) return (void *) &lmp->atom->nimpropers;
  if (strcmp(name,"nlocal") == 0) return (void *) &lmp->atom->nlocal;
  if (strcmp(name,"nghost") == 0) return (void *) &lmp->atom->nghost;
  if (strcmp(name,"nmax") == 0) return (void *) &lmp->atom->nmax;
  if (strcmp(name,"special_lj") == 0) return (void *) lmp->force->special_lj;
  if (strcmp(name,"special_coul") == 0) return (void *) lmp->force->special_coul;

  if (strcmp(name,"q_flag") == 0) return (void *) &lmp->atom->q_flag;

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

  return nullptr;
}

/* ---------------------------------------------------------------------- */

/** Get data type of a LAMMPS per-atom property
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

This function returns an integer that encodes the data type of the per-atom
property with the specified name. See :cpp:enum:`_LMP_DATATYPE_CONST` for valid
values. Callers of :cpp:func:`lammps_extract_atom` can use this information
to then decide how to cast the ``void *`` pointer and access the data.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  name    string with the name of the extracted property
 * \return         integer constant encoding the data type of the property
 *                 or -1 if not found.
 * */

int lammps_extract_atom_datatype(void *handle, const char *name)
{
  auto lmp = (LAMMPS *) handle;
  return lmp->atom->extract_datatype(name);
}

/* ---------------------------------------------------------------------- */

/** Get pointer to a LAMMPS per-atom property.
 *
\verbatim embed:rst

This function returns a pointer to the location of per-atom properties
(and per-atom-type properties in the case of the 'mass' keyword).
Per-atom data is distributed across sub-domains and thus MPI ranks.  The
returned pointer is cast to ``void *`` and needs to be cast to a pointer
of data type that the entity represents.

A table with supported keywords is included in the documentation
of the :cpp:func:`Atom::extract() <LAMMPS_NS::Atom::extract>` function.

.. warning::

   The pointers returned by this function are generally not persistent
   since per-atom data may be re-distributed, re-allocated, and
   re-ordered at every re-neighboring operation.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  name    string with the name of the extracted property
 * \return         pointer (cast to ``void *``) to the location of the
 *                 requested data or ``NULL`` if not found. */

void *lammps_extract_atom(void *handle, const char *name)
{
  auto lmp = (LAMMPS *) handle;
  return lmp->atom->extract(name);
}

// ----------------------------------------------------------------------
// Library functions to access data from computes, fixes, variables in LAMMPS
// ----------------------------------------------------------------------

/** Get pointer to data from a LAMMPS compute.
 *
\verbatim embed:rst

This function returns a pointer to the location of data provided by a
:doc:`compute` instance identified by the compute-ID.  Computes may
provide global, per-atom, or local data, and those may be a scalar, a
vector, or an array or they may provide the information about the
dimensions of the respective data.  Since computes may provide multiple
kinds of data, it is required to set style and type flags representing
what specific data is desired.  This also determines to what kind of
pointer the returned pointer needs to be cast to access the data
correctly.  The function returns ``NULL`` if the compute ID is not found
or the requested data is not available or current. The following table
lists the available options.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Style (see :cpp:enum:`_LMP_STYLE_CONST`)
     - Type (see :cpp:enum:`_LMP_TYPE_CONST`)
     - Returned type
     - Returned data
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_SCALAR
     - ``double *``
     - Global scalar
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_VECTOR
     - ``double *``
     - Global vector
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_ARRAY
     - ``double **``
     - Global array
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_VECTOR
     - ``int *``
     - Length of global vector
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_ROWS
     - ``int *``
     - Rows of global array
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_COLS
     - ``int *``
     - Columns of global array
   * - LMP_STYLE_ATOM
     - LMP_TYPE_VECTOR
     - ``double *``
     - Per-atom value
   * - LMP_STYLE_ATOM
     - LMP_TYPE_ARRAY
     - ``double **``
     - Per-atom vector
   * - LMP_STYLE_ATOM
     - LMP_SIZE_COLS
     - ``int *``
     - Columns in per-atom array, 0 if vector
   * - LMP_STYLE_LOCAL
     - LMP_TYPE_VECTOR
     - ``double *``
     - Local data vector
   * - LMP_STYLE_LOCAL
     - LMP_TYPE_ARRAY
     - ``double **``
     - Local data array
   * - LMP_STYLE_LOCAL
     - LMP_SIZE_VECTOR
     - ``int *``
     - Alias for using LMP_SIZE_ROWS
   * - LMP_STYLE_LOCAL
     - LMP_SIZE_ROWS
     - ``int *``
     - Number of local array rows or length of vector
   * - LMP_STYLE_LOCAL
     - LMP_SIZE_COLS
     - ``int *``
     - Number of local array columns, 0 if vector

.. warning::

   The pointers returned by this function are generally not persistent
   since the computed data may be re-distributed, re-allocated, and
   re-ordered at every invocation. It is advisable to re-invoke this
   function before the data is accessed, or make a copy if the data shall
   be used after other LAMMPS commands have been issued.

.. note::

   If the compute's data is not computed for the current step, the
   compute will be invoked.  LAMMPS cannot easily check at that time, if
   it is valid to invoke a compute, so it may fail with an error.  The
   caller has to check to avoid such an error.


\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  id      string with ID of the compute
 * \param  style   constant indicating the style of data requested
                   (global, per-atom, or local)
 * \param  type    constant indicating type of data (scalar, vector,
                   or array) or size of rows or columns
 * \return         pointer (cast to ``void *``) to the location of the
 *                 requested data or ``NULL`` if not found. */

void *lammps_extract_compute(void *handle, const char *id, int style, int type)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto compute = lmp->modify->get_compute_by_id(id);
    if (!compute) return nullptr;

    if (style == LMP_STYLE_GLOBAL) {
      if (type == LMP_TYPE_SCALAR) {
        if (!compute->scalar_flag) return nullptr;
        if (compute->invoked_scalar != lmp->update->ntimestep)
          compute->compute_scalar();
        return (void *) &compute->scalar;
      }
      if ((type == LMP_TYPE_VECTOR) || (type == LMP_SIZE_VECTOR)) {
        if (!compute->vector_flag) return nullptr;
        if (compute->invoked_vector != lmp->update->ntimestep)
          compute->compute_vector();
        if (type == LMP_TYPE_VECTOR)
          return (void *) compute->vector;
        else
          return (void *) &compute->size_vector;
      }
      if ((type == LMP_TYPE_ARRAY) || (type == LMP_SIZE_ROWS) || (type == LMP_SIZE_COLS)) {
        if (!compute->array_flag) return nullptr;
        if (compute->invoked_array != lmp->update->ntimestep)
          compute->compute_array();
        if (type == LMP_TYPE_ARRAY)
          return (void *) compute->array;
        else if (type == LMP_SIZE_ROWS)
          return (void *) &compute->size_array_rows;
        else
          return (void *) &compute->size_array_cols;
      }
    }

    if (style == LMP_STYLE_ATOM) {
      if (!compute->peratom_flag) return nullptr;
      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();
      if (type == LMP_TYPE_VECTOR) return (void *) compute->vector_atom;
      if (type == LMP_TYPE_ARRAY) return (void *) compute->array_atom;
      if (type == LMP_SIZE_COLS) return (void *) &compute->size_peratom_cols;
    }

    if (style == LMP_STYLE_LOCAL) {
      if (!compute->local_flag) return nullptr;
      if (compute->invoked_local != lmp->update->ntimestep)
        compute->compute_local();
      if (type == LMP_TYPE_SCALAR) return (void *) &compute->size_local_rows;  /* for backward compatibility */
      if (type == LMP_TYPE_VECTOR) return (void *) compute->vector_local;
      if (type == LMP_TYPE_ARRAY) return (void *) compute->array_local;
      if (type == LMP_SIZE_VECTOR) return (void *) &compute->size_local_rows;  /* alias for LMP_SIZE_ROWS */
      if (type == LMP_SIZE_ROWS) return (void *) &compute->size_local_rows;
      if (type == LMP_SIZE_COLS) return (void *) &compute->size_local_cols;
    }
  }
  END_CAPTURE

  return nullptr;
}

/* ---------------------------------------------------------------------- */

/** Get pointer to data from a LAMMPS fix.
 *
\verbatim embed:rst

This function returns a pointer to data provided by a :doc:`fix`
instance identified by its fix-ID.  Fixes may provide global, per-atom,
or local data, and those may be a scalar, a vector, or an array, or they
may provide the information about the dimensions of the respective data.
Since individual fixes may provide multiple kinds of data, it is
required to set style and type flags representing what specific data is
desired.  This also determines to what kind of pointer the returned
pointer needs to be cast to access the data correctly.  The function
returns ``NULL`` if the fix ID is not found or the requested data is not
available.

.. note::

   When requesting global data, the fix data can only be accessed one
   item at a time without access to the pointer itself.  Thus this
   function will allocate storage for a single double value, copy the
   returned value to it, and returns a pointer to the location of the
   copy.  Therefore the allocated storage needs to be freed after its
   use to avoid a memory leak. Example:

   .. code-block:: c

      double *dptr = (double *) lammps_extract_fix(handle,name,0,1,0,0);
      double value = *dptr;
      lammps_free((void *)dptr);

The following table lists the available options.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Style (see :cpp:enum:`_LMP_STYLE_CONST`)
     - Type (see :cpp:enum:`_LMP_TYPE_CONST`)
     - Returned type
     - Returned data
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_SCALAR
     - ``double *``
     - Copy of global scalar
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_VECTOR
     - ``double *``
     - Copy of global vector element at index nrow
   * - LMP_STYLE_GLOBAL
     - LMP_TYPE_ARRAY
     - ``double *``
     - Copy of global array element at nrow, ncol
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_VECTOR
     - ``int *``
     - Length of global vector
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_ROWS
     - ``int *``
     - Rows in global array
   * - LMP_STYLE_GLOBAL
     - LMP_SIZE_COLS
     - ``int *``
     - Columns in global array
   * - LMP_STYLE_ATOM
     - LMP_TYPE_VECTOR
     - ``double *``
     - Per-atom value
   * - LMP_STYLE_ATOM
     - LMP_TYPE_ARRAY
     - ``double **``
     - Per-atom vector
   * - LMP_STYLE_ATOM
     - LMP_SIZE_COLS
     - ``int *``
     - Columns of per-atom array, 0 if vector
   * - LMP_STYLE_LOCAL
     - LMP_TYPE_VECTOR
     - ``double *``
     - Local data vector
   * - LMP_STYLE_LOCAL
     - LMP_TYPE_ARRAY
     - ``double **``
     - Local data array
   * - LMP_STYLE_LOCAL
     - LMP_SIZE_ROWS
     - ``int *``
     - Number of local data rows
   * - LMP_STYLE_LOCAL
     - LMP_SIZE_COLS
     - ``int *``
     - Number of local data columns

.. warning::

   The pointers returned by this function for per-atom or local data are
   generally not persistent, since the computed data may be re-distributed,
   re-allocated, and re-ordered at every invocation of the fix.  It is thus
   advisable to re-invoke this function before the data is accessed, or
   make a copy, if the data shall be used after other LAMMPS commands have
   been issued.

.. note::

   LAMMPS cannot easily check if it is valid to access the data, so it
   may fail with an error.  The caller has to avoid such an error.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  id      string with ID of the fix
 * \param  style   constant indicating the style of data requested
                   (global, per-atom, or local)
 * \param  type    constant indicating type of data (scalar, vector,
                   or array) or size of rows or columns
 * \param  nrow    row index (only used for global vectors and arrays)
 * \param  ncol    column index (only used for global arrays)
 * \return         pointer (cast to ``void *``) to the location of the
 *                 requested data or ``NULL`` if not found. */

void *lammps_extract_fix(void *handle, const char *id, int style, int type,
                         int nrow, int ncol)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) return nullptr;

    if (style == LMP_STYLE_GLOBAL) {
      if (type == LMP_TYPE_SCALAR) {
        if (!fix->scalar_flag) return nullptr;
        auto dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_scalar();
        return (void *) dptr;
      }
      if (type == LMP_TYPE_VECTOR) {
        if (!fix->vector_flag) return nullptr;
        auto dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_vector(nrow);
        return (void *) dptr;
      }
      if (type == LMP_TYPE_ARRAY) {
        if (!fix->array_flag) return nullptr;
        auto dptr = (double *) malloc(sizeof(double));
        *dptr = fix->compute_array(nrow,ncol);
        return (void *) dptr;
      }
      if (type == LMP_SIZE_VECTOR) {
        if (!fix->vector_flag) return nullptr;
        return (void *) &fix->size_vector;
      }
      if ((type == LMP_SIZE_ROWS) || (type == LMP_SIZE_COLS)) {
        if (!fix->array_flag) return nullptr;
        if (type == LMP_SIZE_ROWS)
          return (void *) &fix->size_array_rows;
        else
          return (void *) &fix->size_array_cols;
      }
    }

    if (style == LMP_STYLE_ATOM) {
      if (!fix->peratom_flag) return nullptr;
      if (type == LMP_TYPE_VECTOR) return (void *) fix->vector_atom;
      if (type == LMP_TYPE_ARRAY) return (void *) fix->array_atom;
      if (type == LMP_SIZE_COLS) return (void *) &fix->size_peratom_cols;
    }

    if (style == LMP_STYLE_LOCAL) {
      if (!fix->local_flag) return nullptr;
      if (type == LMP_TYPE_SCALAR) return (void *) &fix->size_local_rows;
      if (type == LMP_TYPE_VECTOR) return (void *) fix->vector_local;
      if (type == LMP_TYPE_ARRAY) return (void *) fix->array_local;
      if (type == LMP_SIZE_ROWS) return (void *) &fix->size_local_rows;
      if (type == LMP_SIZE_COLS) return (void *) &fix->size_local_cols;
    }
  }
  END_CAPTURE

  return nullptr;
}

/* ---------------------------------------------------------------------- */

/** Get pointer to data from a LAMMPS variable.
 *
\verbatim embed:rst

This function returns a pointer to data from a LAMMPS :doc:`variable`
identified by its name.  When the variable is either an *equal*\ -style
compatible variable, a *vector*\ -style variable, or an *atom*\ -style
variable, the variable is evaluated and the corresponding value(s) returned.
Variables of style *internal* are compatible with *equal*\ -style variables and
so are *python*\ -style variables, if they return a numeric value.  For other
variable styles, their string value is returned.  The function returns
``NULL`` when a variable of the provided *name* is not found or of an
incompatible style.  The *group* argument is only used for *atom*\
-style variables and ignored otherwise, with one exception: for style *vector*,
if *group* is "GET_VECTOR_SIZE", the returned pointer will yield the length
of the vector to be returned when dereferenced. This pointer must be
deallocated after the value is read to avoid a memory leak.
If *group* is set to ``NULL`` when extracting data from an *atom*\ -style
variable, the group is assumed to be "all".

When requesting data from an *equal*\ -style or compatible variable
this function allocates storage for a single double value, copies the
returned value to it, and returns a pointer to the location of the
copy.  Therefore the allocated storage needs to be freed after its
use to avoid a memory leak. Example:

.. code-block:: c

   double *dptr = (double *) lammps_extract_variable(handle,name,NULL);
   double value = *dptr;
   lammps_free((void *)dptr);

For *atom*\ -style variables, the return value is a pointer to an
allocated block of storage of double of the length ``atom->nlocal``.
Since the data returned are a copy, the location will persist, but its
content will not be updated in case the variable is re-evaluated.
To avoid a memory leak, this pointer needs to be freed after use in
the calling program.

For *vector*\ -style variables, the returned pointer is to actual LAMMPS data.
The pointer should not be deallocated. Its length depends on the variable,
compute, or fix data used to construct the *vector*\ -style variable.
This length can be fetched by calling this function with *group* set to the
constant "LMP_SIZE_VECTOR", which returns a ``void *`` pointer that can be
dereferenced to an integer that is the length of the vector. This pointer
needs to be deallocated when finished with it to avoid memory leaks.

For other variable styles the returned pointer needs to be cast to
a char pointer. It should not be deallocated.

.. code-block:: c

   const char *cptr = (const char *) lammps_extract_variable(handle,name,NULL);
   printf("The value of variable %s is %s\n", name, cptr);

.. note::

   LAMMPS cannot easily check if it is valid to access the data
   referenced by the variables (e.g., computes, fixes, or thermodynamic
   info), so it may fail with an error.  The caller has to make certain
   that the data are extracted only when it safe to evaluate the variable
   and thus an error or crash are avoided.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  name    name of the variable
 * \param  group   group-ID for atom style variable or ``NULL``
 * \return         pointer (cast to ``void *``) to the location of the
 *                 requested data or ``NULL`` if not found. */

void *lammps_extract_variable(void *handle, const char *name, const char *group)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    int ivar = lmp->input->variable->find(name);
    if (ivar < 0) return nullptr;

    if (lmp->input->variable->equalstyle(ivar)) {
      auto dptr = (double *) malloc(sizeof(double));
      *dptr = lmp->input->variable->compute_equal(ivar);
      return (void *) dptr;
    } else if (lmp->input->variable->atomstyle(ivar)) {
      if (group == nullptr) group = (char *)"all";
      int igroup = lmp->group->find(group);
      if (igroup < 0) return nullptr;
      int nlocal = lmp->atom->nlocal;
      auto vector = (double *) malloc(nlocal*sizeof(double));
      lmp->input->variable->compute_atom(ivar,igroup,vector,1,0);
      return (void *) vector;
    } else if (lmp->input->variable->vectorstyle(ivar)) {
      double *values = nullptr;
      int nvector = lmp->input->variable->compute_vector(ivar, &values);
      if (group != nullptr && strcmp(group,"LMP_SIZE_VECTOR") == 0) {
          int* nvecptr = (int *) malloc(sizeof(int));
          *nvecptr = nvector;
          return (void *) nvecptr;
      } else
        return (void *) values;
    } else {
      return lmp->input->variable->retrieve(name);
    }
  }
  END_CAPTURE
  return nullptr;
}

/* ---------------------------------------------------------------------- */

/** Get data type of a LAMMPS variable.
 *
\verbatim embed:rst

.. versionadded:: 3Nov2022

This function returns an integer that encodes the data type of the variable
with the specified name. See :cpp:enum:`_LMP_VAR_CONST` for valid values.
Callers of :cpp:func:`lammps_extract_variable` can use this information to
decide how to cast the ``void *`` pointer and access the data.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  name    string with the name of the extracted variable
 * \return         integer constant encoding the data type of the property
 *                 or -1 if not found.
 **/

int lammps_extract_variable_datatype(void *handle, const char *name)
{
  auto lmp = (LAMMPS*) handle;

  BEGIN_CAPTURE
  {
    int ivar = lmp->input->variable->find(name);
    if (ivar < 0) return -1;

    if (lmp->input->variable->equalstyle(ivar))
      return LMP_VAR_EQUAL;
    else if (lmp->input->variable->atomstyle(ivar))
      return LMP_VAR_ATOM;
    else if (lmp->input->variable->vectorstyle(ivar))
      return LMP_VAR_VECTOR;
    else
      return LMP_VAR_STRING;
  }
  END_CAPTURE
  return -1;
}

/* ---------------------------------------------------------------------- */

/** Set the value of a string-style variable.
 *
 * This function assigns a new value from the string str to the
 * string-style variable name. Returns -1 if a variable of that
 * name does not exist or is not a string-style variable, otherwise 0.
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  name    name of the variable
 * \param  str     new value of the variable
 * \return         0 on success or -1 on failure
 */
int lammps_set_variable(void *handle, char *name, char *str)
{
  auto lmp = (LAMMPS *) handle;
  int err = -1;

  BEGIN_CAPTURE
  {
    err = lmp->input->variable->set_string(name,str);
  }
  END_CAPTURE

  return err;
}

// ----------------------------------------------------------------------
// Library functions for scatter/gather operations of data
// ----------------------------------------------------------------------

/** Gather the named atom-based entity for all atoms across all processes,
 * in order.
 *
\verbatim embed:rst

This subroutine gathers data for all atoms and stores them in a
one-dimensional array allocated by the user. The data will be ordered by
atom ID, which requires consecutive atom IDs (1 to *natoms*\ ). If you need
a similar array but have non-consecutive atom IDs, see
:cpp:func:`lammps_gather_atoms_concat`; for a similar array but for a subset
of atoms, see :cpp:func:`lammps_gather_atoms_subset`.

The *data* array will be ordered in groups of *count* values, sorted by atom ID
(e.g., if *name* is *x* and *count* = 3, then *data* = x[0][0], x[0][1],
x[0][2], x[1][0], x[1][1], x[1][2], x[2][0], :math:`\dots`);
*data* must be pre-allocated by the caller to length (*count* :math:`\times`
*natoms*), as queried by :cpp:func:`lammps_get_natoms`,
:cpp:func:`lammps_extract_global`, or :cpp:func:`lammps_extract_setting`.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle  pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., *x* or *charge*)
 * \param type    0 for ``int`` values, 1 for ``double`` values
 * \param count   number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with *image* if you want
 *                a single image flag unpacked into (*x*,*y*,*z*) components.
 * \param data    per-atom values packed in a 1-dimensional array of length
 *                *natoms* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
   method:
     alloc and zero count*Natom length vector
     loop over Nlocal to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_gather_atoms(void *handle, const char *name, int type, int count,
                         void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
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
    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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

    } else if (type == 1) {
      double *vector = nullptr;
      double **array = nullptr;
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
    } else {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_atoms: unsupported data type");
      return;
    }
#endif
  }
  END_CAPTURE
}

/** Gather the named atom-based entity for all atoms across all processes,
 * unordered.
 *
\verbatim embed:rst

This subroutine gathers data for all atoms and stores them in a
one-dimensional array allocated by the user. The data will be a concatenation
of chunks from each processor's owned atoms, in whatever order the atoms are
in on each processor. This process has no requirement that the atom IDs be
consecutive. If you need the ID of each atom, you can do another
:cpp:func:`lammps_gather_atoms_concat` call with *name* set to ``id``.
If you have consecutive IDs and want the data to be in order, use
:cpp:func:`lammps_gather_atoms`; for a similar array but for a subset
of atoms, use :cpp:func:`lammps_gather_atoms_subset`.

The *data* array will be in groups of *count* values, with *natoms*
groups total, but not in order by atom ID (e.g., if *name* is *x* and *count*
is 3, then *data* might be something like x[10][0], x[10][1], x[10][2],
x[2][0], x[2][1], x[2][2], x[4][0], :math:`\dots`); *data* must be
pre-allocated by the caller to length (*count* :math:`\times` *natoms*), as
queried by :cpp:func:`lammps_get_natoms`,
:cpp:func:`lammps_extract_global`, or :cpp:func:`lammps_extract_setting`.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name:   desired quantity (e.g., *x* or *charge*\ )
 * \param type:   0 for ``int`` values, 1 for ``double`` values
 * \param count:  number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with "image" if you want
 *                single image flags unpacked into (*x*,*y*,*z*)
 * \param data:   per-atom values packed in a 1-dimensional array of length
 *                *natoms* \* *count*.
 *
 */

/* ----------------------------------------------------------------------
   method:
     Allgather Nlocal atoms from each proc into data
------------------------------------------------------------------------- */

void lammps_gather_atoms_concat(void *handle, const char *name, int type,
                                int count, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_concat() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
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
    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // perform MPI_Allgatherv on each proc's chunk of Nlocal atoms

    int nprocs = lmp->comm->nprocs;

    int *recvcounts,*displs;
    lmp->memory->create(recvcounts,nprocs,"lib/gather:recvcounts");
    lmp->memory->create(displs,nprocs,"lib/gather:displs");

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
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
#endif
  }
  END_CAPTURE
}

/** Gather the named atom-based entity for a subset of atoms.
 *
\verbatim embed:rst

This subroutine gathers data for the requested atom IDs and stores them in a
one-dimensional array allocated by the user. The data will be ordered by atom
ID, but there is no requirement that the IDs be consecutive. If you wish to
return a similar array for *all* the atoms, use :cpp:func:`lammps_gather_atoms`
or :cpp:func:`lammps_gather_atoms_concat`.

The *data* array will be in groups of *count* values, sorted by atom ID
in the same order as the array *ids* (e.g., if *name* is *x*, *count* = 3, and
*ids* is {100, 57, 210}, then *data* might look like {x[100][0], x[100][1],
x[100][2], x[57][0], x[57][1], x[57][2], x[210][0], :math:`\dots`);
*ids* must be provided by the user with length *ndata*, and
*data* must be pre-allocated by the caller to length
(*count* :math:`\times` *ndata*).

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name:   desired quantity (e.g., *x* or *charge*)
 * \param type:   0 for ``int`` values, 1 for ``double`` values
 * \param count:  number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with "image" if you want
 *                single image flags unpacked into (*x*,*y*,*z*)
 * \param ndata:  number of atoms for which to return data (can be all of them)
 * \param ids:    list of *ndata* atom IDs for which to return data
 * \param data:   per-atom values packed in a 1-dimensional array of length
 *                *ndata* \* *count*.
 *
 */

/* ----------------------------------------------------------------------
   method:
     alloc and zero count*Ndata length vector
     loop over Ndata to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_gather_atoms_subset(void *handle, const char *name, int type,
                                int count, int ndata, int *ids, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms_subset: atoms must have mappable ids");
      return;
    }

    void *vptr = lmp->atom->extract(name);
    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_atoms_subset: "
                            "unknown property name");
      return;
    }

    // copy = Ndata length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
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
#endif
  }
  END_CAPTURE
}

/** Scatter the named atom-based entities in *data* to all processes.
 *
\verbatim embed:rst

This subroutine takes data stored in a one-dimensional array supplied by the
user and scatters them to all atoms on all processes. The data must be
ordered by atom ID, with the requirement that the IDs be consecutive.
Use :cpp:func:`lammps_scatter_atoms_subset` to scatter data for some (or all)
atoms, unordered.

The *data* array needs to be ordered in groups of *count* values, sorted by
atom ID (e.g., if *name* is *x* and *count* = 3, then
*data* = {x[0][0], x[0][1], x[0][2], x[1][0], x[1][1], x[1][2], x[2][0],
:math:`\dots`}); *data* must be of length (*count* :math:`\times` *natoms*).

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle  pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., *x* or *charge*)
 * \param type    0 for ``int`` values, 1 for ``double`` values
 * \param count   number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with *image* if you have
 *                a single image flag packed into (*x*,*y*,*z*) components.
 * \param data    per-atom values packed in a one-dimensional array of length
 *                *natoms* \* *count*.
 *
 */

/* ----------------------------------------------------------------------
   method:
     loop over Natoms, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

void lammps_scatter_atoms(void *handle, const char *name, int type, int count,
                          void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset;

    // error if tags are not defined or not consecutive or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms: ids must exist, be consecutive, and be mapped");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lmp->atom->extract(name);
    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,
                            "lammps_scatter_atoms: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      auto dptr = (double *) data;

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
#endif
  }
  END_CAPTURE
}

/** Scatter the named atom-based entities in *data* from a subset of atoms
 *  to all processes.
 *
\verbatim embed:rst

This subroutine takes data stored in a one-dimensional array supplied by the
user and scatters them to a subset of atoms on all processes. The array
*data* contains data associated with atom IDs, but there is no requirement that
the IDs be consecutive, as they are provided in a separate array.
Use :cpp:func:`lammps_scatter_atoms` to scatter data for all atoms, in order.

The *data* array needs to be organized in groups of *count* values, with the
groups in the same order as the array *ids*. For example, if you want *data*
to be the array {x[1][0], x[1][1], x[1][2], x[100][0], x[100][1], x[100][2],
x[57][0], x[57][1], x[57][2]}, then *count* = 3, *ndata* = 3, and *ids* would
be {1, 100, 57}.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name:   desired quantity (e.g., *x* or *charge*)
 * \param type:   0 for ``int`` values, 1 for ``double`` values
 * \param count:  number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with "image" if you have
 *                all the image flags packed into (*xyz*)
 * \param ndata:  number of atoms listed in *ids* and *data* arrays
 * \param ids:    list of *ndata* atom IDs to scatter data to
 * \param data    per-atom values packed in a 1-dimensional array of length
 *                *ndata* \* *count*.
 *
 */

/* ----------------------------------------------------------------------
   scatter the named atom-based entity in data to a subset of atoms
   data is ordered by provided atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see scatter_atoms() to scatter data for all atoms, ordered by consecutive IDs
   name = desired quantity (e.g., x or charge)
   type = 0 for integer values, 1 for double values
   count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
     use count = 3 with "image" for xyz to be packed into single image flag
   ndata = # of atoms in ids and data (could be all atoms)
   ids = list of Ndata atom IDs to scatter data to
   data = atom-based values in 1d data, ordered by count, then by atom ID
     (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
     data must be correct length = count*Ndata
   method:
     loop over Ndata, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

void lammps_scatter_atoms_subset(void *handle, const char *name, int type,
                                 int count, int ndata, int *ids, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms_subset: atoms must have mapped ids");
      return;
    }

    void *vptr = lmp->atom->extract(name);
    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,
                            "lammps_scatter_atoms_subset: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      auto dptr = (double *) data;

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
#endif
  }
  END_CAPTURE
}

/** Gather type and constituent atom info for all bonds
 *
\verbatim embed:rst

.. versionadded:: 28Jul2021

This function copies the list of all bonds into a buffer provided by
the calling code. The buffer will be filled with bond type, bond atom 1,
bond atom 2 for each bond. Thus the buffer has to be allocated to the
dimension of 3 times the **total** number of bonds times the size of
the LAMMPS "tagint" type, which is either 4 or 8 bytes depending on
whether they are stored in 32-bit or 64-bit integers, respectively.
This size depends on the compile time settings used when compiling
the LAMMPS library and can be queried by calling
:cpp:func:`lammps_extract_setting()` with the keyword "tagint".

When running in parallel, the data buffer must be allocated on **all**
MPI ranks and will be filled with the information for **all** bonds
in the system.

Below is a brief C code demonstrating accessing this collected bond information.

.. code-block:: c

   #include "library.h"

   #include <stdint.h>
   #include <stdio.h>
   #include <stdlib.h>

   int main(int argc, char **argv)
   {
       int tagintsize;
       int64_t i, nbonds;
       void *handle, *bonds;

       handle = lammps_open_no_mpi(0, NULL, NULL);
       lammps_file(handle, "in.some_input");

       tagintsize = lammps_extract_setting(handle, "tagint");
       if (tagintsize == 4)
           nbonds = *(int32_t *)lammps_extract_global(handle, "nbonds");
        else
           nbonds = *(int64_t *)lammps_extract_global(handle, "nbonds");
       bonds = malloc(nbonds * 3 * tagintsize);

       lammps_gather_bonds(handle, bonds);

       if (lammps_extract_setting(handle, "world_rank") == 0) {
           if (tagintsize == 4) {
               int32_t *bonds_real = (int32_t *)bonds;
               for (i = 0; i < nbonds; ++i) {
                   printf("bond % 4ld: type = %d, atoms: % 4d  % 4d\n",i,
                          bonds_real[3*i], bonds_real[3*i+1], bonds_real[3*i+2]);
               }
           } else {
               int64_t *bonds_real = (int64_t *)bonds;
               for (i = 0; i < nbonds; ++i) {
                   printf("bond % 4ld: type = %ld, atoms: % 4ld  % 4ld\n",i,
                          bonds_real[3*i], bonds_real[3*i+1], bonds_real[3*i+2]);
               }
           }
       }

       lammps_close(handle);
       lammps_mpi_finalize();
       free(bonds);
       return 0;
   }

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  data    pointer to data to copy the result to */

void lammps_gather_bonds(void *handle, void *data)
{
  auto lmp = (LAMMPS *) handle;
  BEGIN_CAPTURE {
    void *val = lammps_extract_global(handle,"nbonds");
    bigint nbonds = *(bigint *)val;

    // no bonds
    if (nbonds == 0) return;

    // count per MPI rank bonds, determine offsets and allocate local buffers
    int localbonds = lmp->atom->avec->pack_bond(nullptr);
    int nprocs = lmp->comm->nprocs;
    int *bufsizes = new int[nprocs];
    int *bufoffsets = new int[nprocs];
    MPI_Allgather(&localbonds, 1, MPI_INT, bufsizes, 1, MPI_INT, lmp->world);
    bufoffsets[0] = 0;
    bufsizes[0] *= 3;           // 3 items per bond: type, atom1, atom2
    for (int i = 1; i < nprocs; ++i) {
      bufoffsets[i] = bufoffsets[i-1] + bufsizes[i-1];
      bufsizes[i] *= 3;         // 3 items per bond: type, atom1, atom2
    }

    tagint **bonds;
    // add 1 to localbonds, so "bonds" does not become a NULL pointer
    lmp->memory->create(bonds, localbonds+1, 3, "library:gather_bonds:localbonds");
    lmp->atom->avec->pack_bond(bonds);
    MPI_Allgatherv(&bonds[0][0], 3*localbonds, MPI_LMP_TAGINT, data, bufsizes,
                   bufoffsets, MPI_LMP_TAGINT, lmp->world);
    lmp->memory->destroy(bonds);
    delete[] bufsizes;
    delete[] bufoffsets;
  }
  END_CAPTURE
}

/** Gather type and constituent atom info for all angles
 *
\verbatim embed:rst

.. versionadded:: 8Feb2023

This function copies the list of all angles into a buffer provided by
the calling code. The buffer will be filled with angle type, angle atom 1,
angle atom 2, angle atom 3 for each angle. Thus the buffer has to be allocated to the
dimension of 4 times the **total** number of angles times the size of
the LAMMPS "tagint" type, which is either 4 or 8 bytes depending on
whether they are stored in 32-bit or 64-bit integers, respectively.
This size depends on the compile time settings used when compiling
the LAMMPS library and can be queried by calling
:cpp:func:`lammps_extract_setting()` with the keyword "tagint".

When running in parallel, the data buffer must be allocated on **all**
MPI ranks and will be filled with the information for **all** angles
in the system.

Below is a brief C code demonstrating accessing this collected angle information.

.. code-block:: c

   #include "library.h"

   #include <stdint.h>
   #include <stdio.h>
   #include <stdlib.h>

   int main(int argc, char **argv)
   {
       int tagintsize;
       int64_t i, nangles;
       void *handle, *angles;

       handle = lammps_open_no_mpi(0, NULL, NULL);
       lammps_file(handle, "in.some_input");

       tagintsize = lammps_extract_setting(handle, "tagint");
       if (tagintsize == 4)
           nangles = *(int32_t *)lammps_extract_global(handle, "nangles");
        else
           nangles = *(int64_t *)lammps_extract_global(handle, "nangles");
       angles = malloc(nangles * 4 * tagintsize);

       lammps_gather_angles(handle, angles);

       if (lammps_extract_setting(handle, "world_rank") == 0) {
           if (tagintsize == 4) {
               int32_t *angles_real = (int32_t *)angles;
               for (i = 0; i < nangles; ++i) {
                   printf("angle % 4ld: type = %d, atoms: % 4d  % 4d  % 4d\n",i,
                          angles_real[4*i], angles_real[4*i+1], angles_real[4*i+2], angles_real[4*i+3]);
               }
           } else {
               int64_t *angles_real = (int64_t *)angles;
               for (i = 0; i < nangles; ++i) {
                   printf("angle % 4ld: type = %ld, atoms: % 4ld  % 4ld  % 4ld\n",i,
                          angles_real[4*i], angles_real[4*i+1], angles_real[4*i+2], angles_real[4*i+3]);
               }
           }
       }

       lammps_close(handle);
       lammps_mpi_finalize();
       free(angles);
       return 0;
   }

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  data    pointer to data to copy the result to */

void lammps_gather_angles(void *handle, void *data)
{
  auto lmp = (LAMMPS *) handle;
  BEGIN_CAPTURE {
    void *val = lammps_extract_global(handle,"nangles");
    bigint nangles = *(bigint *)val;

    // no angles
    if (nangles == 0) return;

    // count per MPI rank angles, determine offsets and allocate local buffers
    int localangles = lmp->atom->avec->pack_angle(nullptr);
    int nprocs = lmp->comm->nprocs;
    int *bufsizes = new int[nprocs];
    int *bufoffsets = new int[nprocs];
    MPI_Allgather(&localangles, 1, MPI_INT, bufsizes, 1, MPI_INT, lmp->world);
    bufoffsets[0] = 0;
    bufsizes[0] *= 4;           // 4 items per angle: type, atom1, atom2, atom3
    for (int i = 1; i < nprocs; ++i) {
      bufoffsets[i] = bufoffsets[i-1] + bufsizes[i-1];
      bufsizes[i] *= 4;         // 4 items per angle: type, atom1, atom2, atom3
    }

    tagint **angles;
    // add 1 to localangles, so "angles" does not become a NULL pointer
    lmp->memory->create(angles, localangles+1, 4, "library:gather_angles:localangles");
    lmp->atom->avec->pack_angle(angles);
    MPI_Allgatherv(&angles[0][0], 4*localangles, MPI_LMP_TAGINT, data, bufsizes,
                   bufoffsets, MPI_LMP_TAGINT, lmp->world);
    lmp->memory->destroy(angles);
    delete[] bufsizes;
    delete[] bufoffsets;
  }
  END_CAPTURE
}

/** Gather type and constituent atom info for all dihedrals
 *
\verbatim embed:rst

.. versionadded:: 8Feb2023

This function copies the list of all dihedrals into a buffer provided by
the calling code. The buffer will be filled with dihedral type, dihedral atom 1,
dihedral atom 2, dihedral atom 3, dihedral atom 4 for each dihedral.
Thus the buffer has to be allocated to the
dimension of 5 times the **total** number of dihedrals times the size of
the LAMMPS "tagint" type, which is either 4 or 8 bytes depending on
whether they are stored in 32-bit or 64-bit integers, respectively.
This size depends on the compile time settings used when compiling
the LAMMPS library and can be queried by calling
:cpp:func:`lammps_extract_setting()` with the keyword "tagint".

When running in parallel, the data buffer must be allocated on **all**
MPI ranks and will be filled with the information for **all** dihedrals
in the system.

Below is a brief C code demonstrating accessing this collected dihedral information.

.. code-block:: c

   #include "library.h"

   #include <stdint.h>
   #include <stdio.h>
   #include <stdlib.h>

   int main(int argc, char **argv)
   {
       int tagintsize;
       int64_t i, ndihedrals;
       void *handle, *dihedrals;

       handle = lammps_open_no_mpi(0, NULL, NULL);
       lammps_file(handle, "in.some_input");

       tagintsize = lammps_extract_setting(handle, "tagint");
       if (tagintsize == 4)
           ndihedrals = *(int32_t *)lammps_extract_global(handle, "ndihedrals");
        else
           ndihedrals = *(int64_t *)lammps_extract_global(handle, "ndihedrals");
       dihedrals = malloc(ndihedrals * 5 * tagintsize);

       lammps_gather_dihedrals(handle, dihedrals);

       if (lammps_extract_setting(handle, "world_rank") == 0) {
           if (tagintsize == 4) {
               int32_t *dihedrals_real = (int32_t *)dihedrals;
               for (i = 0; i < ndihedrals; ++i) {
                   printf("dihedral % 4ld: type = %d, atoms: % 4d  % 4d  % 4d  % 4d\n",i,
                          dihedrals_real[5*i], dihedrals_real[5*i+1], dihedrals_real[5*i+2], dihedrals_real[5*i+3], dihedrals_real[5*i+4]);
               }
           } else {
               int64_t *dihedrals_real = (int64_t *)dihedrals;
               for (i = 0; i < ndihedrals; ++i) {
                   printf("dihedral % 4ld: type = %ld, atoms: % 4ld  % 4ld  % 4ld  % 4ld\n",i,
                          dihedrals_real[5*i], dihedrals_real[5*i+1], dihedrals_real[5*i+2], dihedrals_real[5*i+3], dihedrals_real[5*i+4]);
               }
           }
       }

       lammps_close(handle);
       lammps_mpi_finalize();
       free(dihedrals);
       return 0;
   }

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  data    pointer to data to copy the result to */

void lammps_gather_dihedrals(void *handle, void *data)
{
  auto lmp = (LAMMPS *) handle;
  BEGIN_CAPTURE {
    void *val = lammps_extract_global(handle,"ndihedrals");
    bigint ndihedrals = *(bigint *)val;

    // no dihedrals
    if (ndihedrals == 0) return;

    // count per MPI rank dihedrals, determine offsets and allocate local buffers
    int localdihedrals = lmp->atom->avec->pack_dihedral(nullptr);
    int nprocs = lmp->comm->nprocs;
    int *bufsizes = new int[nprocs];
    int *bufoffsets = new int[nprocs];
    MPI_Allgather(&localdihedrals, 1, MPI_INT, bufsizes, 1, MPI_INT, lmp->world);
    bufoffsets[0] = 0;
    bufsizes[0] *= 5;           // 5 items per dihedral: type, atom1, atom2, atom3, atom4
    for (int i = 1; i < nprocs; ++i) {
      bufoffsets[i] = bufoffsets[i-1] + bufsizes[i-1];
      bufsizes[i] *= 5;         // 5 items per dihedral: type, atom1, atom2, atom3, atom4
    }

    tagint **dihedrals;
    // add 1 to localdihedrals, so "dihedrals" does not become a NULL pointer
    lmp->memory->create(dihedrals, localdihedrals+1, 5, "library:gather_dihedrals:localdihedrals");
    lmp->atom->avec->pack_dihedral(dihedrals);
    MPI_Allgatherv(&dihedrals[0][0], 5*localdihedrals, MPI_LMP_TAGINT, data, bufsizes,
                   bufoffsets, MPI_LMP_TAGINT, lmp->world);
    lmp->memory->destroy(dihedrals);
    delete[] bufsizes;
    delete[] bufoffsets;
  }
  END_CAPTURE
}

/** Gather type and constituent atom info for all impropers
 *
\verbatim embed:rst

.. versionadded:: 8Feb2023

This function copies the list of all impropers into a buffer provided by
the calling code. The buffer will be filled with improper type, improper atom 1,
improper atom 2, improper atom 3, improper atom 4 for each improper.
Thus the buffer has to be allocated to the
dimension of 5 times the **total** number of impropers times the size of
the LAMMPS "tagint" type, which is either 4 or 8 bytes depending on
whether they are stored in 32-bit or 64-bit integers, respectively.
This size depends on the compile time settings used when compiling
the LAMMPS library and can be queried by calling
:cpp:func:`lammps_extract_setting()` with the keyword "tagint".

When running in parallel, the data buffer must be allocated on **all**
MPI ranks and will be filled with the information for **all** impropers
in the system.

Below is a brief C code demonstrating accessing this collected improper information.

.. code-block:: c

   #include "library.h"

   #include <stdint.h>
   #include <stdio.h>
   #include <stdlib.h>

   int main(int argc, char **argv)
   {
       int tagintsize;
       int64_t i, nimpropers;
       void *handle, *impropers;

       handle = lammps_open_no_mpi(0, NULL, NULL);
       lammps_file(handle, "in.some_input");

       tagintsize = lammps_extract_setting(handle, "tagint");
       if (tagintsize == 4)
           nimpropers = *(int32_t *)lammps_extract_global(handle, "nimpropers");
        else
           nimpropers = *(int64_t *)lammps_extract_global(handle, "nimpropers");
       impropers = malloc(nimpropers * 5 * tagintsize);

       lammps_gather_impropers(handle, impropers);

       if (lammps_extract_setting(handle, "world_rank") == 0) {
           if (tagintsize == 4) {
               int32_t *impropers_real = (int32_t *)impropers;
               for (i = 0; i < nimpropers; ++i) {
                   printf("improper % 4ld: type = %d, atoms: % 4d  % 4d  % 4d  % 4d\n",i,
                          impropers_real[5*i], impropers_real[5*i+1], impropers_real[5*i+2], impropers_real[5*i+3], impropers_real[5*i+4]);
               }
           } else {
               int64_t *impropers_real = (int64_t *)impropers;
               for (i = 0; i < nimpropers; ++i) {
                   printf("improper % 4ld: type = %ld, atoms: % 4ld  % 4ld  % 4ld  % 4ld\n",i,
                          impropers_real[5*i], impropers_real[5*i+1], impropers_real[5*i+2], impropers_real[5*i+3], impropers_real[5*i+4]);
               }
           }
       }

       lammps_close(handle);
       lammps_mpi_finalize();
       free(impropers);
       return 0;
   }

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \param  data    pointer to data to copy the result to */

void lammps_gather_impropers(void *handle, void *data)
{
  auto lmp = (LAMMPS *) handle;
  BEGIN_CAPTURE {
    void *val = lammps_extract_global(handle,"nimpropers");
    bigint nimpropers = *(bigint *)val;

    // no impropers
    if (nimpropers == 0) return;

    // count per MPI rank impropers, determine offsets and allocate local buffers
    int localimpropers = lmp->atom->avec->pack_improper(nullptr);
    int nprocs = lmp->comm->nprocs;
    int *bufsizes = new int[nprocs];
    int *bufoffsets = new int[nprocs];
    MPI_Allgather(&localimpropers, 1, MPI_INT, bufsizes, 1, MPI_INT, lmp->world);
    bufoffsets[0] = 0;
    bufsizes[0] *= 5;           // 5 items per improper: type, atom1, atom2, atom3, atom4
    for (int i = 1; i < nprocs; ++i) {
      bufoffsets[i] = bufoffsets[i-1] + bufsizes[i-1];
      bufsizes[i] *= 5;         // 5 items per improper: type, atom1, atom2, atom3, atom4
    }

    tagint **impropers;
    // add 1 to localimpropers, so "impropers" does not become a NULL pointer
    lmp->memory->create(impropers, localimpropers+1, 5, "library:gather_impropers:localimpropers");
    lmp->atom->avec->pack_improper(impropers);
    MPI_Allgatherv(&impropers[0][0], 5*localimpropers, MPI_LMP_TAGINT, data, bufsizes,
                   bufoffsets, MPI_LMP_TAGINT, lmp->world);
    lmp->memory->destroy(impropers);
    delete[] bufsizes;
    delete[] bufoffsets;
  }
  END_CAPTURE
}

/** Gather the named per-atom, per-atom fix, per-atom compute, or fix property/atom-based entities
 *  from all processes, in order by atom ID.
 *
\verbatim embed:rst

This subroutine gathers data from all processes and stores them in a one-dimensional array
allocated by the user. The array *data* will be ordered by atom ID, which requires consecutive IDs
(1 to *natoms*\ ). If you need a similar array but for non-consecutive atom IDs, see
:cpp:func:`lammps_gather_concat`; for a similar array but for a subset of atoms, see
:cpp:func:`lammps_gather_subset`.

The *data* array will be ordered in groups of *count* values, sorted by atom ID (e.g., if *name* is
*x*, then *data* is {x[0][0], x[0][1], x[0][2], x[1][0], x[1][1], x[1][2], x[2][0],
:math:`\dots`}); *data* must be pre-allocated by the caller to the correct length
(*count*\ :math:`{}\times{}`\ *natoms*), as queried by :cpp:func:`lammps_get_natoms`,
:cpp:func:`lammps_extract_global`, or :cpp:func:`lammps_extract_setting`.

This function will return an error if fix or compute data are requested and the fix or compute ID
given does not have per-atom data.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle  pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., "x" or "f" for atom properties, "f_id" for per-atom fix
 *                data, "c_id" for per-atom compute data, "d_name" or "i_name" for fix
 *                property/atom vectors with *count* = 1, "d2_name" or "i2_name" for fix
 *                property/atom vectors with *count* > 1)
 * \param type    0 for ``int`` values, 1 for ``double`` values
 * \param count   number of per-atom values (e.g., 1 for *type* or *charge*, 3 for *x* or *f*);
 *                use *count* = 3 with *image* if you want the image flags unpacked into
 *                (*x*,*y*,*z*) components.
 * \param data    per-atom values packed into a one-dimensional array of length
 *                *natoms* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
  Contributing author: Thomas Swinburne (CNRS & CINaM, Marseille, France)
  gather the named atom-based entity for all atoms
    return it in user-allocated data
  data will be ordered by atom ID
    requirement for consecutive atom IDs (1 to N)
  see gather_concat() to return data for all atoms, unordered
  see gather_subset() to return data for only a subset of atoms
  name = "x" , "f" or other atom properties
         "f_fix", "c_compute" for fixes / computes
         "d_name" or "i_name" for fix property/atom vectors with count = 1
         "d2_name" or "i2_name" for fix property/atom arrays with count > 1
         will return error if fix/compute isn't atom-based
  type = 0 for integer values, 1 for double values
  count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
    use count = 3 with "image" if want single image flag unpacked into xyz
  return atom-based values in 1d data, ordered by count, then by atom ID
    (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
    data must be pre-allocated by caller to correct length
    correct length = count*Natoms, as queried by get_natoms()
  method:
    alloc and zero count*Natom length vector
    loop over Nlocal to fill vector with my values
    Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_gather(void *handle, const char *name, int type, int count, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,offset,ltype;

    // error if tags are not defined or not consecutive

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);
    void *vptr = lmp->atom->extract(name);

    // fix

    if (vptr==nullptr && utils::strmatch(name,"^f_")) {

      auto fix = lmp->modify->get_fix_by_id(&name[2]);
      if (!fix) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: unknown fix id");
        return;
      }

      if (fix->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: fix does not return peratom data");
        return;
      }
      if ((count > 1) && (fix->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: count != values peratom for fix");
        return;
      }

      if (lmp->update->ntimestep % fix->peratom_freq) {
        if (lmp->comm->me == 0)
          lmp->error->all(FLERR,"lammps_gather: fix not computed at compatible time");
        return;
      }

      if (count==1) vptr = (void *) fix->vector_atom;
      else vptr = (void *) fix->array_atom;
    }

    // compute

    if (vptr==nullptr && utils::strmatch(name,"^c_")) {

      auto compute = lmp->modify->get_compute_by_id(&name[2]);
      if (!compute) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: unknown compute id");
        return;
      }

      if (compute->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: compute does not return peratom data");
        return;
      }
      if ((count > 1) && (compute->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: count != values peratom for compute");
        return;
      }

      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();

      if (count==1) vptr = (void *) compute->vector_atom;
      else vptr = (void *) compute->array_atom;
    }

    // custom fix property/atom vector or array

    if ((vptr == nullptr) && utils::strmatch(name,"^[id]2?_")) {

      int idx,icol;
      if (utils::strmatch(name,"^[id]_")) idx = lmp->atom->find_custom(&name[2],ltype,icol);
      else idx = lmp->atom->find_custom(&name[3],ltype,icol);

      if (idx < 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: unknown property/atom id");
        return;
      }

      if (ltype != type) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: mismatch property/atom type");
        return;
      }
      if ((count == 1) && (icol != 0)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: mismatch property/atom count");
        return;
      }
      if ((count > 1) && (icol != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: mismatch property/atom count");
        return;
      }

      if (count == 1) {
        if (ltype==0) vptr = (void *) lmp->atom->ivector[idx];
        else vptr = (void *) lmp->atom->dvector[idx];
      } else {
        if (ltype==0) vptr = (void *) lmp->atom->iarray[idx];
        else vptr = (void *) lmp->atom->darray[idx];
      }
    }

    // no match

    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather: undefined property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type==0) {
      int *vector = nullptr;
      int **array = nullptr;

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
      double *vector = nullptr;
      double **array = nullptr;
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
#endif
  }
  END_CAPTURE
}

/** Gather the named per-atom, per-atom fix, per-atom compute, or fix property/atom-based entities
 *  from all processes, unordered.
 *
\verbatim embed:rst

This subroutine gathers data for all atoms and stores them in a one-dimensional array allocated by
the user. The data will be a concatenation of chunks from each processor's owned atoms, in
whatever order the atoms are in on each processor. This process has no requirement that the atom
IDs be consecutive. If you need the ID of each atom, you can do another call to either
:cpp:func:`lammps_gather_atoms_concat` or :cpp:func:`lammps_gather_concat` with *name* set to
``id``. If you have consecutive IDs and want the data to be in order, use
:cpp:func:`lammps_gather`; for a similar array but for a subset of atoms, use
:cpp:func:`lammps_gather_subset`.

The *data* array will be in groups of *count* values, with *natoms* groups total, but not in order
by atom ID (e.g., if *name* is *x* and *count* is 3, then *data* might be something like
{x[10][0], x[10][1], x[10][2], x[2][0], x[2][1], x[2][2], x[4][0], :math:`\dots`}); *data* must be
pre-allocated by the caller to length (*count* :math:`\times` *natoms*), as queried by
:cpp:func:`lammps_get_natoms`, :cpp:func:`lammps_extract_global`, or
:cpp:func:`lammps_extract_setting`.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name:   desired quantity (e.g., "x" or "f" for atom properties, "f_id" for per-atom fix
 *                data, "c_id" for per-atom compute data, "d_name" or "i_name" for fix
 *                property/atom vectors with count = 1, "d2_name" or "i2_name" for fix
 *                property/atom vectors with count > 1)
 * \param type:   0 for ``int`` values, 1 for ``double`` values
 * \param count:  number of per-atom values (e.g., 1 for *type* or *charge*, 3 for *x* or *f*);
 *                use *count* = 3 with *image* if you want the image flags unpacked into
 *                (*x*,*y*,*z*) components.
 * \param data:   per-atom values packed into a one-dimensional array of length
 *                *natoms* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
  Contributing author: Thomas Swinburne (CNRS & CINaM, Marseille, France)
  gather the named atom-based entity for all atoms
    return it in user-allocated data
  data will be ordered by atom ID
    requirement for consecutive atom IDs (1 to N)
  see gather() to return data ordered by consecutive atom IDs
  see gather_subset() to return data for only a subset of atoms
  name = "x" , "f" or other atom properties
         "f_fix", "c_compute" for fixes / computes
         "d_name" or "i_name" for fix property/atom vectors with count = 1
         "d2_name" or "i2_name" for fix property/atom arrays with count > 1
         will return error if fix/compute isn't atom-based
  type = 0 for integer values, 1 for double values
  count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
    use count = 3 with "image" if want single image flag unpacked into xyz
  return atom-based values in 1d data, ordered by count, then by atom ID
    (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
    data must be pre-allocated by caller to correct length
    correct length = count*Natoms, as queried by get_natoms()
  method:
    alloc and zero count*Natom length vector
    loop over Nlocal to fill vector with my values
    Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_gather_concat(void *handle, const char *name, int type, int count,
                          void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather_concat"
                          " not compatible with -DLAMMPS_BIGBIG");
#else
    int i,offset,ltype;

    // error if tags are not defined or not consecutive

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_concat");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);
    void *vptr = lmp->atom->extract(name);

    // fix

    if (vptr==nullptr && utils::strmatch(name,"^f_")) {

      auto fix = lmp->modify->get_fix_by_id(&name[2]);
      if (!fix) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: unknown fix id");
        return;
      }

      if (fix->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: fix does not return peratom data");
        return;
      }
      if ((count > 1) && (fix->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: count != values peratom for fix");
        return;
      }
      if (lmp->update->ntimestep % fix->peratom_freq) {
        if (lmp->comm->me == 0)
          lmp->error->all(FLERR,"lammps_gather_concat: fix not computed at compatible time");
        return;
      }

      if (count==1) vptr = (void *) fix->vector_atom;
      else vptr = (void *) fix->array_atom;
    }

    // compute

    if (vptr==nullptr && utils::strmatch(name,"^c_")) {

      auto compute = lmp->modify->get_compute_by_id(&name[2]);
      if (!compute) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: unknown compute id");
        return;
      }

      if (compute->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: compute does not return peratom data");
        return;
      }
      if ((count > 1) && (compute->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: count != values peratom for compute");
        return;
      }

      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();

      if (count==1) vptr = (void *) compute->vector_atom;
      else vptr = (void *) compute->array_atom;
    }

    // custom per-atom vector or array

    if ((vptr==nullptr) && utils::strmatch(name,"^[id]2?_")) {

      int idx,icol;
      if (utils::strmatch(name,"^[id]_")) idx = lmp->atom->find_custom(&name[2],ltype,icol);
      else idx = lmp->atom->find_custom(&name[3],ltype,icol);

      if (idx < 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: unknown property/atom id");
        return;
      }

      if (ltype != type) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: mismatch property/atom type");
        return;
      }
      if ((count == 1) && (icol != 0)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: mismatch property/atom count");
        return;
      }
      if ((count > 1) && (icol != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_concat: mismatch property/atom count");
        return;
      }

      if (count == 1) {
        if (ltype==0) vptr = (void *) lmp->atom->ivector[idx];
        else vptr = (void *) lmp->atom->dvector[idx];
      } else {
        if (ltype==0) vptr = (void *) lmp->atom->iarray[idx];
        else vptr = (void *) lmp->atom->darray[idx];
      }
    }

    // no match

    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_concat: undefined property name");
      return;
    }

    // perform MPI_Allgatherv on each proc's chunk of Nlocal atoms

    int nprocs = lmp->comm->nprocs;

    int *recvcounts,*displs;
    lmp->memory->create(recvcounts,nprocs,"lib/gather:recvcounts");
    lmp->memory->create(displs,nprocs,"lib/gather:displs");

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;

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
      double *vector = nullptr;
      double **array = nullptr;
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
#endif
  }
  END_CAPTURE
}

/** Gather the named per-atom, per-atom fix, per-atom compute, or fix property/atom-based entities
 *  from all processes for a subset of atoms.
 *
\verbatim embed:rst

This subroutine gathers data for the requested atom IDs and stores them in a one-dimensional array
allocated by the user. The data will be ordered by atom ID, but there is no requirement that the
IDs be consecutive. If you wish to return a similar array for *all* the atoms, use
:cpp:func:`lammps_gather` or :cpp:func:`lammps_gather_concat`.

The *data* array will be in groups of *count* values, sorted by atom ID in the same order as the
array *ids* (e.g., if *name* is *x*, *count* = 3, and *ids* is {100, 57, 210}, then *data* might
look like {x[100][0], x[100][1], x[100][2], x[57][0], x[57][1], x[57][2], x[210][0],
:math:`\dots`}); *ids* must be provided by the user with length *ndata*, and *data* must be
pre-allocated by the caller to length (*count*\ :math:`{}\times{}`\ *ndata*).

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., "x" or "f" for atom properties, "f_id" for per-atom fix
 *                data, "c_id" for per-atom compute data, "d_name" or "i_name" for fix
 *                property/atom vectors with *count* = 1, "d2_name" or "i2_name" for fix
 *                property/atom vectors with *count* > 1)
 * \param type    0 for ``int`` values, 1 for ``double`` values
 * \param count   number of per-atom values (e.g., 1 for *type* or *charge*, 3 for *x* or *f*);
 *                use *count* = 3 with *image* if you want the image flags unpacked into
 *                (*x*,*y*,*z*) components.
 * \param ndata:  number of atoms for which to return data (can be all of them)
 * \param ids:    list of *ndata* atom IDs for which to return data
 * \param data    per-atom values packed into a one-dimensional array of length
 *                *ndata* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
  Contributing author: Thomas Swinburne (CNRS & CINaM, Marseille, France)
  gather the named atom-based entity for all atoms
    return it in user-allocated data
  data will be ordered by atom ID
    requirement for consecutive atom IDs (1 to N)
  see gather() to return data ordered by consecutive atom IDs
  see gather_concat() to return data for all atoms, unordered
  name = "x" , "f" or other atom properties
         "f_fix", "c_compute" for fixes / computes
         "d_name" or "i_name" for fix property/atom vectors with count = 1
         "d2_name" or "i2_name" for fix property/atom arrays with count > 1
         will return error if fix/compute isn't atom-based
  type = 0 for integer values, 1 for double values
  count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
    use count = 3 with "image" if want single image flag unpacked into xyz
  return atom-based values in 1d data, ordered by count, then by atom ID
    (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
    data must be pre-allocated by caller to correct length
    correct length = count*Natoms, as queried by get_natoms()
  method:
    alloc and zero count*Natom length vector
    loop over Nlocal to fill vector with my values
    Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_gather_subset(void *handle, const char *name, int type, int count,
                          int ndata, int *ids, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_gather_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset,ltype;
    tagint id;

    // error if tags are not defined or no atom map

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_subset");
      return;
    }

    void *vptr = lmp->atom->extract(name);

    // fix

    if (vptr==nullptr && utils::strmatch(name,"^f_")) {

      auto fix = lmp->modify->get_fix_by_id(&name[2]);
      if (!fix) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: unknown fix id");
        return;
      }

      if (fix->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: fix does not return peratom data");
        return;
      }
      if ((count > 1) && (fix->size_peratom_cols != count)) {
        lmp->error->warning(FLERR,"lammps_gather_subset: count != values peratom for fix");
        return;
      }
      if (lmp->update->ntimestep % fix->peratom_freq) {
        if (lmp->comm->me == 0)
          lmp->error->all(FLERR,"lammps_gather_subset: fix not computed at compatible time");
        return;
      }

      if (count==1) vptr = (void *) fix->vector_atom;
      else vptr = (void *) fix->array_atom;
    }

    // compute

    if (vptr==nullptr && utils::strmatch(name,"^c_")) {

      auto compute = lmp->modify->get_compute_by_id(&name[2]);
      if (!compute) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: unknown compute id");
        return;
      }

      if (compute->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: compute does not return peratom data");
        return;
      }
      if ((count > 1) && (compute->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: count != values peratom for compute");
        return;
      }

      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();

      if (count==1) vptr = (void *) compute->vector_atom;
      else vptr = (void *) compute->array_atom;
    }

    // custom fix property/atom vector or array

    if ((vptr == nullptr) && utils::strmatch(name,"^[id]2?_")) {

      int idx,icol;
      if (utils::strmatch(name,"^[id]_")) idx = lmp->atom->find_custom(&name[2],ltype,icol);
      else idx = lmp->atom->find_custom(&name[3],ltype,icol);

      if (idx < 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: unknown property/atom id");
        return;
      }

      if (ltype != type) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: mismatch property/atom type");
        return;
      }
      if (count == 1 && icol != 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: mismatch property/atom count");
        return;
      }
      if (count > 1 && icol != count) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather_subset: mismatch property/atom count");
        return;
      }

      if (count == 1) {
        if (ltype==0) vptr = (void *) lmp->atom->ivector[idx];
        else vptr = (void *) lmp->atom->dvector[idx];
      } else {
        if (ltype==0) vptr = (void *) lmp->atom->iarray[idx];
        else vptr = (void *) lmp->atom->darray[idx];
      }
    }

    // no match

    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_gather_subset: undefined property name");
      return;
    }

    // copy = Ndata length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;

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
#endif
  }
  END_CAPTURE
}

/** Scatter the named per-atom, per-atom fix, per-atom compute, or fix property/atom-based
 *  entity in *data* to all processes.
 *
\verbatim embed:rst

This subroutine takes data stored in a one-dimensional array supplied by the user and scatters
them to all atoms on all processes. The data must be ordered by atom ID, with the requirement that
the IDs be consecutive. Use :cpp:func:`lammps_scatter_subset` to scatter data for some (or all)
atoms, unordered.

The *data* array needs to be ordered in groups of *count* values, sorted by atom ID (e.g., if
*name* is *x* and *count* = 3, then *data* = {x[0][0], x[0][1], x[0][2], x[1][0], x[1][1],
x[1][2], x[2][0], :math:`\dots`}); *data* must be of length (*count* :math:`\times` *natoms*).

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle  pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., "x" or "f" for atom properties, "f_id" for per-atom fix
 *                data, "c_id" for per-atom compute data, "d_name" or "i_name" for fix
 *                property/atom vectors with *count* = 1, "d2_name" or "i2_name" for fix
 *                property/atom vectors with *count* > 1)
 * \param type    0 for ``int`` values, 1 for ``double`` values
 * \param count   number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with *image* if you have
 *                a single image flag packed into (*x*,*y*,*z*) components.
 * \param data    per-atom values packed in a one-dimensional array of length
 *                *natoms* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
  Contributing author: Thomas Swinburne (CNRS & CINaM, Marseille, France)
  scatter the named atom-based entity in data to all atoms
  data will be ordered by atom ID
    requirement for consecutive atom IDs (1 to N)
  see scatter_subset() to scatter data for some (or all) atoms, unordered
  name = "x" , "f" or other atom properties
         "f_fix", "c_compute" for fixes / computes
         "d_name" or "i_name" for fix property/atom vectors with count = 1
         "d2_name" or "i2_name" for fix property/atom arrays with count > 1
         will return error if fix/compute isn't atom-based
  type = 0 for integer values, 1 for double values
  count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
    use count = 3 with "image" if want single image flag unpacked into xyz
  return atom-based values in 1d data, ordered by count, then by atom ID
    (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
    data must be pre-allocated by caller to correct length
    correct length = count*Natoms, as queried by get_natoms()
  method:
    alloc and zero count*Natom length vector
    loop over Nlocal to fill vector with my values
    Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

void lammps_scatter(void *handle, const char *name, int type, int count,
                    void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_scatter() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset,ltype;

    // error if tags are not defined or not consecutive or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);
    void *vptr = lmp->atom->extract(name);

    // fix

    if (vptr==nullptr && utils::strmatch(name,"^f_")) {

      auto fix = lmp->modify->get_fix_by_id(&name[2]);
      if (!fix) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: unknown fix id");
        return;
      }

      if (fix->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: fix does not return peratom data");
        return;
      }
      if ((count > 1) && (fix->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: count != values peratom for fix");
        return;
      }

      if (count==1) vptr = (void *) fix->vector_atom;
      else vptr = (void *) fix->array_atom;
    }

    // compute

    if (vptr==nullptr && utils::strmatch(name,"^c_")) {

      auto compute = lmp->modify->get_compute_by_id(&name[2]);
      if (!compute) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: unknown compute id");
        return;
      }

      if (compute->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: compute does not return peratom data");
        return;
      }
      if ((count > 1) && (compute->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: count != values peratom for compute");
        return;
      }

      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();

      if (count==1) vptr = (void *) compute->vector_atom;
      else vptr = (void *) compute->array_atom;
    }

    // custom fix property/atom vector or array

    if ((vptr == nullptr) && utils::strmatch(name,"^[id]2?_")) {

      int idx,icol;
      if (utils::strmatch(name,"^[id]_")) idx = lmp->atom->find_custom(&name[2],ltype,icol);
      else idx = lmp->atom->find_custom(&name[3],ltype,icol);

      if (idx < 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: unknown property/atom id");
        return;
      }

      if (ltype != type) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: mismatch property/atom type");
        return;
      }
      if (count == 1 && icol != 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: mismatch property/atom count");
        return;
      }
      if (count > 1 && icol != count) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter: mismatch property/atom count");
        return;
      }

      if (count == 1) {
        if (ltype==0) vptr = (void *) lmp->atom->ivector[idx];
        else vptr = (void *) lmp->atom->dvector[idx];
      } else {
        if (ltype==0) vptr = (void *) lmp->atom->iarray[idx];
        else vptr = (void *) lmp->atom->darray[idx];
      }
    }

    // no match

    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_scatter: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      auto dptr = (double *) data;

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
#endif
  }
  END_CAPTURE
}

/** Scatter the named per-atom, per-atom fix, per-atom compute, or fix property/atom-based
 * entities in *data* from a subset of atoms to all processes.
 *
\verbatim embed:rst

This subroutine takes data stored in a one-dimensional array supplied by the
user and scatters them to a subset of atoms on all processes. The array
*data* contains data associated with atom IDs, but there is no requirement that
the IDs be consecutive, as they are provided in a separate array.
Use :cpp:func:`lammps_scatter` to scatter data for all atoms, in order.

The *data* array needs to be organized in groups of *count* values, with the
groups in the same order as the array *ids*. For example, if you want *data*
to be the array {x[1][0], x[1][1], x[1][2], x[100][0], x[100][1], x[100][2],
x[57][0], x[57][1], x[57][2]}, then *count* = 3, *ndata* = 3, and *ids* would
be {1, 100, 57}.

This function is not compatible with ``-DLAMMPS_BIGBIG``.

\endverbatim
 *
 * \param handle: pointer to a previously created LAMMPS instance
 * \param name    desired quantity (e.g., "x" or "f" for atom properties, "f_id" for per-atom fix
 *                data, "c_id" for per-atom compute data, "d_name" or "i_name" for fix
 *                property/atom vectors with *count* = 1, "d2_name" or "i2_name" for fix
 *                property/atom vectors with *count* > 1)
 * \param type:   0 for ``int`` values, 1 for ``double`` values
 * \param count:  number of per-atom values (e.g., 1 for *type* or *charge*,
 *                3 for *x* or *f*); use *count* = 3 with "image" if you want
 *                single image flags unpacked into (*x*,*y*,*z*)
 * \param ndata:  number of atoms listed in *ids* and *data* arrays
 * \param ids:    list of *ndata* atom IDs to scatter data to
 * \param data    per-atom values packed in a 1-dimensional array of length
 *                *ndata* \* *count*.
 *
 */
/* ----------------------------------------------------------------------
  Contributing author: Thomas Swinburne (CNRS & CINaM, Marseille, France)
   scatter the named atom-based entity in data to a subset of atoms
   data is ordered by provided atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see scatter_atoms() to scatter data for all atoms, ordered by consecutive IDs
   name = "x" , "f" or other atom properties
          "d_name" or "i_name" for fix property/atom quantities
          "f_fix", "c_compute" for fixes / computes
          will return error if fix/compute doesn't isn't atom-based
   type = 0 for integer values, 1 for double values
   count = # of per-atom values (e.g., 1 for type or charge, 3 for x or f)
     use count = 3 with "image" for xyz to be packed into single image flag
   ndata = # of atoms in ids and data (could be all atoms)
   ids = list of Ndata atom IDs to scatter data to
   data = atom-based values in 1d data, ordered by count, then by atom ID
     (e.g., x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...)
     data must be correct length = count*Ndata
   method:
     loop over Ndata, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

void lammps_scatter_subset(void *handle, const char *name,int type, int count,
                                 int ndata, int *ids, void *data)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
#if defined(LAMMPS_BIGBIG)
    lmp->error->all(FLERR,"Library function lammps_scatter_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
#else
    int i,j,m,offset,ltype;
    tagint id;

    // error if tags are not defined or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == Atom::MAP_NONE) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms_subset");
      return;
    }

    void *vptr = lmp->atom->extract(name);

    // fix

    if (vptr==nullptr && utils::strmatch(name,"^f_")) {

      auto fix = lmp->modify->get_fix_by_id(&name[2]);
      if (!fix) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: unknown fix id");
        return;
      }

      if (fix->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: fix does not return peratom data");
        return;
      }
      if ((count > 1) && (fix->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: count != values peratom for fix");
        return;
      }

      if (count==1) vptr = (void *) fix->vector_atom;
      else vptr = (void *) fix->array_atom;
    }

    // compute

    if (vptr==nullptr && utils::strmatch(name,"^c_")) {

      auto compute = lmp->modify->get_compute_by_id(&name[2]);
      if (!compute) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: unknown compute id");
        return;
      }

      if (compute->peratom_flag == 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: compute does not return peratom data");
        return;
      }
      if ((count > 1) && (compute->size_peratom_cols != count)) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: count != values peratom for compute");
        return;
      }

      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();

      if (count==1) vptr = (void *) compute->vector_atom;
      else vptr = (void *) compute->array_atom;
    }

    // custom fix property/atom vector or array

    if ((vptr == nullptr) && utils::strmatch(name,"^[id]2?_")) {

      int idx,icol;
      if (utils::strmatch(name,"^[id]_")) idx = lmp->atom->find_custom(&name[2],ltype,icol);
      else idx = lmp->atom->find_custom(&name[3],ltype,icol);

      if (idx < 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: unknown property/atom id");
        return;
      }

      if (ltype != type) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_scatter_subset: mismatch property/atom type");
        return;
      }
      if (count == 1 && icol != 0) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: mismatch property/atom count");
        return;
      }
      if (count > 1 && icol != count) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,"lammps_gather: mismatch property/atom count");
        return;
      }

      if (count == 1) {
        if (ltype==0) vptr = (void *) lmp->atom->ivector[idx];
        else vptr = (void *) lmp->atom->dvector[idx];
      } else {
        if (ltype==0) vptr = (void *) lmp->atom->iarray[idx];
        else vptr = (void *) lmp->atom->darray[idx];
      }
    }

    // no match

    if (vptr == nullptr) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"lammps_scatter_atoms_subset: "
                                  "unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = nullptr;
      int **array = nullptr;
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
      double *vector = nullptr;
      double **array = nullptr;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      auto dptr = (double *) data;

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
#endif
  }
  END_CAPTURE
}

/* ---------------------------------------------------------------------- */

/** Create N atoms from list of coordinates
 *
\verbatim embed:rst

The prototype for this function when compiling with ``-DLAMMPS_BIGBIG``
is:

.. code-block:: c

   int lammps_create_atoms(void *handle, int n, int64_t *id, int *type, double *x, double *v, int64_t *image, int bexpand);

This function creates additional atoms from a given list of coordinates
and a list of atom types.  Additionally the atom-IDs, velocities, and
image flags may be provided.  If atom-IDs are not provided, they will be
automatically created as a sequence following the largest existing
atom-ID.

This function is useful to add atoms to a simulation or - in tandem with
:cpp:func:`lammps_reset_box` - to restore a previously extracted and
saved state of a simulation.  Additional properties for the new atoms
can then be assigned via the :cpp:func:`lammps_scatter_atoms`
:cpp:func:`lammps_extract_atom` functions.

For non-periodic boundaries, atoms will **not** be created that have
coordinates outside the box unless it is a shrink-wrap boundary and the
shrinkexceed flag has been set to a non-zero value.  For periodic
boundaries atoms will be wrapped back into the simulation cell and its
image flags adjusted accordingly, unless explicit image flags are
provided.

The function returns the number of atoms created or -1 on failure (e.g.,
when called before as box has been created).

Coordinates and velocities have to be given in a 1d-array in the order
X(1),Y(1),Z(1),X(2),Y(2),Z(2),...,X(N),Y(N),Z(N).

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance
 * \param  n        number of atoms, N, to be added to the system
 * \param  id       pointer to N atom IDs; ``NULL`` will generate IDs
 * \param  type     pointer to N atom types (required)
 * \param  x        pointer to 3N doubles with x-,y-,z- positions
                    of the new atoms (required)
 * \param  v        pointer to 3N doubles with x-,y-,z- velocities
                    of the new atoms (set to 0.0 if ``NULL``)
 * \param  image    pointer to N imageint sets of image flags, or ``NULL``
 * \param  bexpand  if 1, atoms outside of shrink-wrap boundaries will
                    still be created and not dropped and the box extended
 * \return          number of atoms created on success;
                    -1 on failure (no box, no atom IDs, etc.) */

int lammps_create_atoms(void *handle, int n, const tagint *id, const int *type,
                        const double *x, const double *v, const imageint *image,
                        int bexpand)
{
  auto lmp = (LAMMPS *) handle;
  bigint natoms_prev = lmp->atom->natoms;

  BEGIN_CAPTURE
  {
    // error if box does not exist or tags not defined

    int flag = 0;
    std::string msg("Failure in lammps_create_atoms: ");
    if (lmp->domain->box_exist == 0) {
      flag = 1;
      msg += "trying to create atoms before before simulation box is defined";
    }
    if (lmp->atom->tag_enable == 0) {
      flag = 1;
      msg += "must have atom IDs to use this function";
    }

    if (flag) {
      if (lmp->comm->me == 0) lmp->error->warning(FLERR,msg);
      return -1;
    }

    // loop over all N atoms on all MPI ranks
    // if this proc would own it based on its coordinates, invoke create_atom()
    // optionally set atom tags and velocities

    Atom *atom = lmp->atom;
    Domain *domain = lmp->domain;
    int nlocal = atom->nlocal;

    int nlocal_prev = nlocal;
    double xdata[3];
    imageint idata, *img;

    for (int i = 0; i < n; i++) {
      xdata[0] = x[3*i];
      xdata[1] = x[3*i+1];
      xdata[2] = x[3*i+2];
      if (image) {
        idata = image[i];
        img = &idata;
      } else img = nullptr;
      const tagint tag = id ? id[i] : 0;

      // create atom only on MPI rank that would own it

      if (!domain->ownatom(tag, xdata, img, bexpand)) continue;

      atom->avec->create_atom(type[i],xdata);
      if (id) atom->tag[nlocal] = id[i];
      else atom->tag[nlocal] = 0;
      if (v) {
        atom->v[nlocal][0] = v[3*i];
        atom->v[nlocal][1] = v[3*i+1];
        atom->v[nlocal][2] = v[3*i+2];
      }
      if (image) atom->image[nlocal] = image[i];
      nlocal++;
    }

    // if no tags are given explicitly, create new and unique tags

    if (id == nullptr) atom->tag_extend();

    // reset box info, if extended when adding atoms.

    if (bexpand) domain->reset_box();

    // need to reset atom->natoms inside LAMMPS

    bigint ncurrent = nlocal;
    MPI_Allreduce(&ncurrent,&lmp->atom->natoms,1,MPI_LMP_BIGINT,
                  MPI_SUM,lmp->world);

    // init per-atom fix/compute/variable values for created atoms

    atom->data_fix_compute_variable(nlocal_prev,nlocal);

    // if global map exists, reset it
    // invoke map_init() b/c atom count has grown

    if (lmp->atom->map_style != Atom::MAP_NONE) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }
  }
  END_CAPTURE;
  return (int) lmp->atom->natoms - natoms_prev;
}

// ----------------------------------------------------------------------
// Library functions for accessing neighbor lists
// ----------------------------------------------------------------------

/** Find index of a neighbor list requested by a pair style
 *
 * This function determines which of the available neighbor lists for
 * pair styles matches the given conditions.  It first matches the style
 * name. If exact is 1 the name must match exactly, if exact is 0, a
 * regular expression or sub-string match is done.  If the pair style is
 * hybrid or hybrid/overlay the style is matched against the sub styles
 * instead.
 * If a the same pair style is used multiple times as a sub-style, the
 * nsub argument must be > 0 and represents the nth instance of the sub-style
 * (same as for the pair_coeff command, for example).  In that case
 * nsub=0 will not produce a match and this function will return -1.
 *
 * The final condition to be checked is the request ID (reqid).  This
 * will normally be 0, but some pair styles request multiple neighbor
 * lists and set the request ID to a value > 0.
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  style    String used to search for pair style instance
 * \param  exact    Flag to control whether style should match exactly or only
 *                  a regular expression / sub-string match is applied.
 * \param  nsub     match nsub-th hybrid sub-style instance of the same style
 * \param  reqid    request id to identify neighbor list in case there are
 *                  multiple requests from the same pair style instance
 * \return          return neighbor list index if found, otherwise -1 */

int lammps_find_pair_neighlist(void *handle, const char *style, int exact, int nsub, int reqid) {
  auto lmp = (LAMMPS *) handle;
  Pair *pair = lmp->force->pair_match(style, exact, nsub);

  if (pair != nullptr) {
    // find neigh list
    for (int i = 0; i < lmp->neighbor->nlist; i++) {
      NeighList *list = lmp->neighbor->lists[i];
      if ((list->requestor_type == NeighList::PAIR)
           && (pair == list->requestor)
           && (list->id == reqid) ) return i;
    }
  }
  return -1;
}

/* ---------------------------------------------------------------------- */

/** Find index of a neighbor list requested by a fix
 *
 * The neighbor list request from a fix is identified by the fix ID and
 * the request ID.  The request ID is typically 0, but will be > 0 in
 * case a fix has multiple neighbor list requests.
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param id       Identifier of fix instance
 * \param reqid    request id to identify neighbor list in case there are
 *                 multiple requests from the same fix
 * \return         return neighbor list index if found, otherwise -1  */

int lammps_find_fix_neighlist(void *handle, const char *id, int reqid) {
  auto lmp = (LAMMPS *) handle;
  auto fix = lmp->modify->get_fix_by_id(id);
  if (!fix) return -1;

  // find neigh list
  for (int i = 0; i < lmp->neighbor->nlist; i++) {
    NeighList *list = lmp->neighbor->lists[i];
    if ((list->requestor_type == NeighList::FIX)
         && (fix == list->requestor)
         && (list->id == reqid) ) return i;
  }
  return -1;
}

/* ---------------------------------------------------------------------- */

/** Find index of a neighbor list requested by a compute
 *
 * The neighbor list request from a compute is identified by the compute
 * ID and the request ID.  The request ID is typically 0, but will be
 * > 0 in case a compute has multiple neighbor list requests.
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param id       Identifier of compute instance
 * \param reqid    request id to identify neighbor list in case there are
 *                 multiple requests from the same compute
 * \return         return neighbor list index if found, otherwise -1 */

int lammps_find_compute_neighlist(void *handle, const char *id, int reqid) {
  auto lmp = (LAMMPS *) handle;
  auto compute = lmp->modify->get_compute_by_id(id);
  if (!compute) return -1;

  // find neigh list
  for (int i = 0; i < lmp->neighbor->nlist; i++) {
    NeighList * list = lmp->neighbor->lists[i];
    if ((list->requestor_type == NeighList::COMPUTE)
         && (compute == list->requestor)
         && (list->id == reqid) ) return i;
  }
  return -1;
}

/* ---------------------------------------------------------------------- */

/** Return the number of entries in the neighbor list with given index
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param idx      neighbor list index
 * \return         return number of entries in neighbor list, -1 if idx is
 *                 not a valid index
 */
int lammps_neighlist_num_elements(void *handle, int idx) {
  auto   lmp = (LAMMPS *) handle;
  Neighbor * neighbor = lmp->neighbor;

  if (idx < 0 || idx >= neighbor->nlist) {
    return -1;
  }

  NeighList * list = neighbor->lists[idx];
  return list->inum;
}

/* ---------------------------------------------------------------------- */

/** Return atom local index, number of neighbors, and array of neighbor local
 * atom indices of neighbor list entry
 *
 * \param handle          pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param idx             index of this neighbor list in the list of all neighbor lists
 * \param element         index of this neighbor list entry
 * \param[out] iatom      local atom index (i.e. in the range [0, nlocal + nghost), -1 if
                          invalid idx or element value
 * \param[out] numneigh   number of neighbors of atom iatom or 0
 * \param[out] neighbors  pointer to array of neighbor atom local indices or NULL */

void lammps_neighlist_element_neighbors(void *handle, int idx, int element, int *iatom, int *numneigh, int **neighbors) {
  auto   lmp = (LAMMPS *) handle;
  Neighbor * neighbor = lmp->neighbor;
  *iatom = -1;
  *numneigh = 0;
  *neighbors = nullptr;

  if (idx < 0 || idx >= neighbor->nlist) {
    return;
  }

  NeighList * list = neighbor->lists[idx];

  if (element < 0 || element >= list->inum) {
    return;
  }

  int i = list->ilist[element];
  *iatom     = i;
  *numneigh  = list->numneigh[i];
  *neighbors = list->firstneigh[i];
}

// ----------------------------------------------------------------------
// Library functions for accessing LAMMPS configuration
// ----------------------------------------------------------------------

/** Get numerical representation of the LAMMPS version date.
 *
\verbatim embed:rst

The :cpp:func:`lammps_version` function returns an integer representing
the version of the LAMMPS code in the format YYYYMMDD.  This can be used
to implement backward compatibility in software using the LAMMPS library
interface.  The specific format guarantees, that this version number is
growing with every new LAMMPS release.

\endverbatim
 *
 * \param  handle  pointer to a previously created LAMMPS instance
 * \return         an integer representing the version data in the
 *                 format YYYYMMDD */

int lammps_version(void *handle)
{
  auto lmp = (LAMMPS *) handle;
  return lmp->num_ver;
}

/** Get operating system and architecture information
 *
\verbatim embed:rst

.. versionadded:: 9Oct2020

The :cpp:func:`lammps_get_os_info` function can be used to retrieve
detailed information about the hosting operating system and
compiler/runtime.

A suitable buffer for a C-style string has to be provided and its length.
The assembled text will be truncated to not overflow this buffer. The
string is typically a few hundred bytes long.

\endverbatim
 *
 * \param  buffer    string buffer to copy the information to
 * \param  buf_size  size of the provided string buffer */

/* ---------------------------------------------------------------------- */

void lammps_get_os_info(char *buffer, int buf_size)
{
  if (buf_size <=0) return;
  buffer[0] = buffer[buf_size-1] = '\0';
  std::string txt = platform::os_info() + "\n";
  txt += platform::compiler_info();
  txt += " with " + platform::openmp_standard() + "\n";
  strncpy(buffer, txt.c_str(), buf_size-1);
}

/* ---------------------------------------------------------------------- */

/** This function is used to query whether LAMMPS was compiled with
 *  a real MPI library or in serial. For the real MPI library it
 *  reports the size of the MPI communicator in bytes (4 or 8),
 *  which allows to check for compatibility with a hosting code.
 *
 * \return 0 when compiled with MPI STUBS, otherwise the MPI_Comm size in bytes */

int lammps_config_has_mpi_support()
{
#ifdef MPI_STUBS
  return 0;
#else
  return sizeof(MPI_Comm);
#endif
}

/* ---------------------------------------------------------------------- */

/** Check if the LAMMPS library supports reading or writing compressed
 * files via a pipe to gzip or similar compression programs

\verbatim embed:rst
Several LAMMPS commands (e.g., :doc:`read_data`, :doc:`write_data`,
:doc:`dump styles atom, custom, and xyz <dump>`) support reading and
writing compressed files via creating a pipe to the ``gzip`` program.
This function checks whether this feature was :ref:`enabled at compile
time <gzip>`. It does **not** check whether``gzip`` or any other
supported compression programs themselves are installed and usable.
\endverbatim
 *
 * \return 1 if yes, otherwise 0
 */
int lammps_config_has_gzip_support() {
  return Info::has_gzip_support() ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Check if the LAMMPS library supports writing PNG format images

\verbatim embed:rst
The LAMMPS :doc:`dump style image <dump_image>` supports writing multiple
image file formats.  Most of them, however, need support from an external
library, and using that has to be :ref:`enabled at compile time <graphics>`.
This function checks whether support for the `PNG image file format
<https://en.wikipedia.org/wiki/Portable_Network_Graphics>`_ is available
in the current LAMMPS library.
\endverbatim
 *
 * \return 1 if yes, otherwise 0
 */
int lammps_config_has_png_support() {
  return Info::has_png_support() ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Check if the LAMMPS library supports writing JPEG format images

\verbatim embed:rst
The LAMMPS :doc:`dump style image <dump_image>` supports writing multiple
image file formats.  Most of them, however, need support from an external
library, and using that has to be :ref:`enabled at compile time <graphics>`.
This function checks whether support for the `JPEG image file format
<https://jpeg.org/jpeg/>`_ is available in the current LAMMPS library.
\endverbatim
 *
 * \return 1 if yes, otherwise 0
 */
int lammps_config_has_jpeg_support() {
  return Info::has_jpeg_support() ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Check if the LAMMPS library supports creating movie files via a pipe to ffmpeg

\verbatim embed:rst
The LAMMPS :doc:`dump style movie <dump_image>` supports generating movies
from images on-the-fly via creating a pipe to the
`ffmpeg <https://ffmpeg.org/>`_ program.
This function checks whether this feature was :ref:`enabled at compile time <graphics>`.
It does **not** check whether the ``ffmpeg`` itself is installed and usable.
\endverbatim
 *
 * \return 1 if yes, otherwise 0
 */
int lammps_config_has_ffmpeg_support() {
  return Info::has_ffmpeg_support() ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Check whether LAMMPS errors will throw C++ exceptions.
 *
\verbatim embed:rst

.. deprecated:: TBD

   LAMMPS has now exceptions always enabled, so this function
   will now always return 1 and can be removed from applications
   using the library interface.

In case of an error, LAMMPS will either abort or throw a C++ exception.
The latter has to be :ref:`enabled at compile time <exceptions>`.
This function checks if exceptions were enabled.

When using the library interface with C++ exceptions enabled,
the library interface functions will "catch" them and the
error status can then be checked by calling
:cpp:func:`lammps_has_error` and the most recent error message
can be retrieved via :cpp:func:`lammps_get_last_error_message`.
This can allow to restart a calculation or delete and recreate
the LAMMPS instance when C++ exceptions are enabled.  One application
of using exceptions this way is the :ref:`lammps_shell`.  If C++
exceptions are disabled and an error happens during a call to
LAMMPS, the application will terminate.
\endverbatim
 * \return 1 if yes, otherwise 0
 */
int lammps_config_has_exceptions() {
  return Info::has_exceptions() ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Check whether a specific package has been included in LAMMPS
 *
\verbatim embed:rst
This function checks whether the LAMMPS library in use includes the
specific :doc:`LAMMPS package <Packages>` provided as argument.
\endverbatim
 *
 * \param name string with the name of the package
 * \return 1 if included, 0 if not.
 */
int lammps_config_has_package(const char *name) {
  return Info::has_package(name) ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Count the number of installed packages in the LAMMPS library.
 *
\verbatim embed:rst
This function counts how many :doc:`LAMMPS packages <Packages>` are
included in the LAMMPS library in use.
\endverbatim
 *
 * \return number of packages included
 */
int lammps_config_package_count() {
  int i = 0;
  while (LAMMPS::installed_packages[i] != nullptr) {
    ++i;
  }
  return i;
}

/* ---------------------------------------------------------------------- */

/** Get the name of a package in the list of installed packages in the LAMMPS library.
 *
\verbatim embed:rst
This function copies the name of the package with the index *idx* into the
provided C-style string buffer.  The length of the buffer must be provided
as *buf_size* argument.  If the name of the package exceeds the length of the
buffer, it will be truncated accordingly.  If the index is out of range,
the function returns 0 and *buffer* is set to an empty string, otherwise 1;
\endverbatim
 *
 * \param idx index of the package in the list of included packages (0 <= idx < package count)
 * \param buffer string buffer to copy the name of the package to
 * \param buf_size size of the provided string buffer
 * \return 1 if successful, otherwise 0
 */
int lammps_config_package_name(int idx, char *buffer, int buf_size) {
  int maxidx = lammps_config_package_count();
  if ((idx < 0) || (idx >= maxidx)) {
      buffer[0] = '\0';
      return 0;
  }

  strncpy(buffer, LAMMPS::installed_packages[idx], buf_size);
  return 1;
}

/** Check for compile time settings in accelerator packages included in LAMMPS.
 *
\verbatim embed:rst
This function checks availability of compile time settings of included
:doc:`accelerator packages <Speed_packages>` in LAMMPS.
Supported packages names are "GPU", "KOKKOS", "INTEL", and "OPENMP".
Supported categories are "api" with possible settings "cuda", "hip", "phi",
"pthreads", "opencl", "openmp", and "serial", and "precision" with
possible settings "double", "mixed", and "single".  If the combination
of package, category, and setting is available, the function returns 1,
otherwise 0.
\endverbatim
 *
 * \param  package   string with the name of the accelerator package
 * \param  category  string with the category name of the setting
 * \param  setting   string with the name of the specific setting
 * \return 1 if available, 0 if not.
 */
int lammps_config_accelerator(const char *package,
                              const char *category,
                              const char *setting)
{
  return Info::has_accelerator_feature(package,category,setting) ? 1 : 0;
}

/** Check for presence of a viable GPU package device
 *
\verbatim embed:rst

.. versionadded:: 14May2021

The :cpp:func:`lammps_has_gpu_device` function checks at runtime if
an accelerator device is present that can be used with the
:doc:`GPU package <Speed_gpu>`. If at least one suitable device is
present the function will return 1, otherwise 0.

More detailed information about the available device or devices can
be obtained by calling the
:cpp:func:`lammps_get_gpu_device_info` function.

\endverbatim
 *
 * \return  1 if viable device is available, 0 if not.  */

int lammps_has_gpu_device()
{
  return Info::has_gpu_device() ? 1: 0;
}

/** Get GPU package device information
 *
\verbatim embed:rst

.. versionadded:: 14May2021

The :cpp:func:`lammps_get_gpu_device_info` function can be used to retrieve
detailed information about any accelerator devices that are viable for use
with the :doc:`GPU package <Speed_gpu>`.  It will produce a string that is
equivalent to the output of the ``nvc_get_device`` or ``ocl_get_device`` or
``hip_get_device`` tools that are compiled alongside LAMMPS if the GPU
package is enabled.

A suitable buffer for a C-style string has to be provided and its length.
The assembled text will be truncated to not overflow this buffer.  This
string can be several kilobytes long, if multiple devices are present.

\endverbatim
 *
 * \param  buffer    string buffer to copy the information to
 * \param  buf_size  size of the provided string buffer */

void lammps_get_gpu_device_info(char *buffer, int buf_size)
{
  if (buf_size <= 0) return;
  buffer[0] = buffer[buf_size-1] = '\0';
  std::string devinfo = Info::get_gpu_device_info();
  strncpy(buffer, devinfo.c_str(), buf_size-1);
}

/* ---------------------------------------------------------------------- */

/** Check if a specific style has been included in LAMMPS
 *
\verbatim embed:rst
This function checks if the LAMMPS library in use includes the
specific *style* of a specific *category* provided as an argument.
Valid categories are: *atom*\ , *integrate*\ , *minimize*\ ,
*pair*\ , *bond*\ , *angle*\ , *dihedral*\ , *improper*\ , *kspace*\ ,
*compute*\ , *fix*\ , *region*\ , *dump*\ , and *command*\ .
\endverbatim
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  category  category of the style
 * \param  name      name of the style
 * \return           1 if included, 0 if not.
 */
int lammps_has_style(void *handle, const char *category, const char *name) {
  auto lmp = (LAMMPS *) handle;
  Info info(lmp);
  return info.has_style(category, name) ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

/** Count the number of styles of category in the LAMMPS library.
 *
\verbatim embed:rst
This function counts how many styles in the provided *category*
are included in the LAMMPS library in use.
Please see :cpp:func:`lammps_has_style` for a list of valid
categories.
\endverbatim
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param category category of styles
 * \return number of styles in category
 */
int lammps_style_count(void *handle, const char *category) {
  auto lmp = (LAMMPS *) handle;
  Info info(lmp);
  return info.get_available_styles(category).size();
}

/* ---------------------------------------------------------------------- */

/** Look up the name of a style by index in the list of style of a given category in the LAMMPS library.
 *
 *
 * This function copies the name of the *category* style with the index
 * *idx* into the provided C-style string buffer.  The length of the buffer
 * must be provided as *buf_size* argument.  If the name of the style
 * exceeds the length of the buffer, it will be truncated accordingly.
 * If the index is out of range, the function returns 0 and *buffer* is
 * set to an empty string, otherwise 1.
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param category category of styles
 * \param idx      index of the style in the list of *category* styles (0 <= idx < style count)
 * \param buffer   string buffer to copy the name of the style to
 * \param buf_size size of the provided string buffer
 * \return 1 if successful, otherwise 0
 */
int lammps_style_name(void *handle, const char *category, int idx,
                      char *buffer, int buf_size) {
  auto lmp = (LAMMPS *) handle;
  Info info(lmp);
  auto styles = info.get_available_styles(category);

  if ((idx >= 0) && (idx < (int) styles.size())) {
    strncpy(buffer, styles[idx].c_str(), buf_size);
    return 1;
  }

  buffer[0] = '\0';
  return 0;
}

/* ---------------------------------------------------------------------- */

/** Check if a specific ID exists in the current LAMMPS instance
 *
\verbatim embed:rst

.. versionadded:: 9Oct2020

This function checks if the current LAMMPS instance a *category* ID of
the given *name* exists.  Valid categories are: *compute*\ , *dump*\ ,
*fix*\ , *group*\ , *molecule*\ , *region*\ , and *variable*\ .

\endverbatim
 *
 * \param  handle    pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  category  category of the id
 * \param  name      name of the id
 * \return           1 if included, 0 if not.
 */
int lammps_has_id(void *handle, const char *category, const char *name) {
  auto lmp = (LAMMPS *) handle;

  if (strcmp(category,"compute") == 0) {
    if (lmp->modify->get_compute_by_id(name)) return 1;
  } else if (strcmp(category,"dump") == 0) {
    if (lmp->output->get_dump_by_id(name)) return 1;
  } else if (strcmp(category,"fix") == 0) {
    if (lmp->modify->get_fix_by_id(name)) return 1;
  } else if (strcmp(category,"group") == 0) {
    if (lmp->group->find(name) >= 0) return 1;
  } else if (strcmp(category,"molecule") == 0) {
    if (lmp->atom->find_molecule(name) >= 0) return 1;
  } else if (strcmp(category,"region") == 0) {
    if (lmp->domain->get_region_by_id(name)) return 1;
  } else if (strcmp(category,"variable") == 0) {
    if (lmp->input->variable->find(name) >= 0) return 1;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

/** Count the number of IDs of a category.
 *
\verbatim embed:rst

.. versionadded:: 9Oct2020

This function counts how many IDs in the provided *category*
are defined in the current LAMMPS instance.
Please see :cpp:func:`lammps_has_id` for a list of valid
categories.

\endverbatim
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param category category of IDs
 * \return number of IDs in category
 */
int lammps_id_count(void *handle, const char *category) {
  auto lmp = (LAMMPS *) handle;
  if (strcmp(category,"compute") == 0) {
    return lmp->modify->get_compute_list().size();
  } else if (strcmp(category,"dump") == 0) {
    return lmp->output->get_dump_list().size();
  } else if (strcmp(category,"fix") == 0) {
    return lmp->modify->get_fix_list().size();
  } else if (strcmp(category,"group") == 0) {
    return lmp->group->ngroup;
  } else if (strcmp(category,"molecule") == 0) {
    return lmp->atom->nmolecule;
  } else if (strcmp(category,"region") == 0) {
    return lmp->domain->get_region_list().size();
  } else if (strcmp(category,"variable") == 0) {
    return lmp->input->variable->nvar;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

/** Look up the name of an ID by index in the list of IDs of a given category.
 *
\verbatim embed:rst

.. versionadded:: 9Oct2020

This function copies the name of the *category* ID with the index
*idx* into the provided C-style string buffer.  The length of the buffer
must be provided as *buf_size* argument.  If the name of the style
exceeds the length of the buffer, it will be truncated accordingly.
If the index is out of range, the function returns 0 and *buffer* is
set to an empty string, otherwise 1.

\endverbatim
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param category category of IDs
 * \param idx      index of the ID in the list of *category* styles (0 <= idx < count)
 * \param buffer   string buffer to copy the name of the style to
 * \param buf_size size of the provided string buffer
 * \return 1 if successful, otherwise 0
 */
int lammps_id_name(void *handle, const char *category, int idx, char *buffer, int buf_size) {
  auto lmp = (LAMMPS *) handle;
  if (idx < 0) return 0;

  if (strcmp(category,"compute") == 0) {
    auto icompute = lmp->modify->get_compute_by_index(idx);
    if (icompute) {
      strncpy(buffer, icompute->id, buf_size);
      return 1;
    }
  } else if (strcmp(category,"dump") == 0) {
    auto idump = lmp->output->get_dump_by_index(idx);
    if (idump) {
      strncpy(buffer, idump->id, buf_size);
      return 1;
    }
  } else if (strcmp(category,"fix") == 0) {
    auto ifix = lmp->modify->get_fix_by_index(idx);
    if (ifix) {
      strncpy(buffer, ifix->id, buf_size);
      return 1;
    }
  } else if (strcmp(category,"group") == 0) {
    if ((idx >= 0) && (idx < lmp->group->ngroup)) {
      strncpy(buffer, lmp->group->names[idx], buf_size);
      return 1;
    }
  } else if (strcmp(category,"molecule") == 0) {
    if ((idx >= 0) && (idx < lmp->atom->nmolecule)) {
      strncpy(buffer, lmp->atom->molecules[idx]->id, buf_size);
      return 1;
    }
  } else if (strcmp(category,"region") == 0) {
    auto regions = lmp->domain->get_region_list();
    if ((idx >= 0) && (idx < (int) regions.size())) {
      strncpy(buffer, regions[idx]->id, buf_size);
      return 1;
    }
  } else if (strcmp(category,"variable") == 0) {
    if ((idx >= 0) && (idx < lmp->input->variable->nvar)) {
      strncpy(buffer, lmp->input->variable->names[idx], buf_size);
      return 1;
    }
  }
  buffer[0] = '\0';
  return 0;
}

/* ---------------------------------------------------------------------- */

/** Count the number of loaded plugins
 *
\verbatim embed:rst

.. versionadded:: 10Mar2021

This function counts how many plugins are currently loaded.

\endverbatim
 *
 * \return number of loaded plugins
 */
int lammps_plugin_count()
{
#if defined(LMP_PLUGIN)
  return plugin_get_num_plugins();
#else
  return 0;
#endif
}

/* ---------------------------------------------------------------------- */

/** Look up the info of a loaded plugin by its index in the list of plugins
 *
\verbatim embed:rst

.. versionadded:: 10Mar2021

This function copies the name of the *style* plugin with the index
*idx* into the provided C-style string buffer.  The length of the buffer
must be provided as *buf_size* argument.  If the name of the style
exceeds the length of the buffer, it will be truncated accordingly.
If the index is out of range, the function returns 0 and *buffer* is
set to an empty string, otherwise 1.

\endverbatim
 *
 * \param  idx       index of the plugin in the list all or *style* plugins
 * \param  stylebuf  string buffer to copy the style of the plugin to
 * \param  namebuf   string buffer to copy the name of the plugin to
 * \param  buf_size  size of the provided string buffers
 * \return 1 if successful, otherwise 0
 */
int lammps_plugin_name(int idx, char *stylebuf, char *namebuf, int buf_size)
{
#if defined(LMP_PLUGIN)
  stylebuf[0] = namebuf[0] = '\0';

  const lammpsplugin_t *plugin = plugin_get_info(idx);
  if (plugin) {
    strncpy(stylebuf,plugin->style,buf_size);
    strncpy(namebuf,plugin->name,buf_size);
    return 1;
  }
#endif
  return 0;
}

// ----------------------------------------------------------------------
// utility functions
// ----------------------------------------------------------------------

/** Encode three integer image flags into a single imageint.
 *
\verbatim embed:rst

The prototype for this function when compiling with ``-DLAMMPS_BIGBIG``
is:

.. code-block:: c

   int64_t lammps_encode_image_flags(int ix, int iy, int iz);

This function performs the bit-shift, addition, and bit-wise OR
operations necessary to combine the values of three integers
representing the image flags in x-, y-, and z-direction.  Unless
LAMMPS is compiled with -DLAMMPS_BIGBIG, those integers are
limited 10-bit signed integers [-512, 511].  Otherwise the return
type changes from ``int`` to ``int64_t`` and the valid range for
the individual image flags becomes [-1048576,1048575],
i.e. that of a 21-bit signed integer.  There is no check on whether
the arguments conform to these requirements.

\endverbatim
 *
 * \param  ix  image flag value in x
 * \param  iy  image flag value in y
 * \param  iz  image flag value in z
 * \return     encoded image flag integer */

imageint lammps_encode_image_flags(int ix, int iy, int iz)
{
  imageint image = ((imageint) (ix + IMGMAX) & IMGMASK) |
    (((imageint) (iy + IMGMAX) & IMGMASK) << IMGBITS) |
    (((imageint) (iz + IMGMAX) & IMGMASK) << IMG2BITS);
  return image;
}

/* ---------------------------------------------------------------------- */

/** Decode a single image flag integer into three regular integers
 *
\verbatim embed:rst

The prototype for this function when compiling with ``-DLAMMPS_BIGBIG``
is:

.. code-block:: c

   void lammps_decode_image_flags(int64_t image, int *flags);

This function does the reverse operation of
:cpp:func:`lammps_encode_image_flags` and takes an image flag integer
does the bit-shift and bit-masking operations to decode it and stores
the resulting three regular integers into the buffer pointed to by
*flags*.

\endverbatim
 *
 * \param  image  encoded image flag integer
 * \param  flags  pointer to storage where the decoded image flags are stored. */

void lammps_decode_image_flags(imageint image, int *flags)
{
  flags[0] = (image & IMGMASK) - IMGMAX;
  flags[1] = (image >> IMGBITS & IMGMASK) - IMGMAX;
  flags[2] = (image >> IMG2BITS) - IMGMAX;
}

/* ---------------------------------------------------------------------- */

/** Set up the callback function for a fix external instance with the given ID.

\verbatim embed:rst

Fix :doc:`external <fix_external>` allows programs that are running LAMMPS through
its library interface to modify certain LAMMPS properties on specific
timesteps, similar to the way other fixes do.

This function sets the callback function for use with the "pf/callback"
mode. The function has to have C language bindings with the prototype:

.. code-block:: c

   void func(void *ptr, bigint timestep, int nlocal, tagint *ids, double **x, double **fexternal);

The argument *ptr* to this function will be stored in fix external and
the passed as the first argument calling the callback function `func()`.
This would usually be a pointer to the active LAMMPS instance, i.e. the same
pointer as the *handle* argument.  This would be needed to call
functions that set the global or per-atom energy or virial contributions
from within the callback function.

The callback mechanism is one of the two modes of how forces and can be
applied to a simulation with the help of fix external. The alternative
is the array mode where you call :cpp:func:`lammps_fix_external_get_force`.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

.. versionchanged:: 28Jul2021

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  funcptr  pointer to callback function
 * \param  ptr      pointer to object in calling code, passed to callback function as first argument */

void lammps_set_fix_external_callback(void *handle, const char *id, FixExternalFnPtr funcptr, void *ptr)
{
  auto lmp = (LAMMPS *) handle;
  auto  callback = (FixExternal::FnPtr) funcptr;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Cannot find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style 'external'", id);

    auto fext = dynamic_cast<FixExternal *>(fix);
    fext->set_callback(callback, ptr);
  }
  END_CAPTURE
}

/** Get pointer to the force array storage in a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

Fix :doc:`external <fix_external>` allows programs that are running
LAMMPS through its library interface to add or modify certain LAMMPS
properties on specific timesteps, similar to the way other fixes do.

This function provides access to the per-atom force storage in a fix
external instance with the given fix-ID to be added to the individual
atoms when using the "pf/array" mode.  The *fexternal* array can be
accessed like other "native" per-atom arrays accessible via the
:cpp:func:`lammps_extract_atom` function.  Please note that the array
stores holds the forces for *local* atoms for each MPI ranks, in the
order determined by the neighbor list build.  Because the underlying
data structures can change as well as the order of atom as they migrate
between MPI processes because of the domain decomposition
parallelization, this function should be always called immediately
before the forces are going to be set to get an up-to-date pointer.
You can use, for example, :cpp:func:`lammps_extract_setting` to obtain
the number of local atoms `nlocal` and then assume the dimensions of
the returned force array as ``double force[nlocal][3]``.

This is an alternative to the callback mechanism in fix external set up
by :cpp:func:`lammps_set_fix_external_callback`. The main difference is
that this mechanism can be used when forces are be pre-computed and the
control alternates between LAMMPS and the external code, while the
callback mechanism can call the external code to compute the force when
the fix is triggered and needs them.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle     pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id         fix ID of fix external instance
 * \return            a pointer to the per-atom force array allocated by the fix */

double **lammps_fix_external_get_force(void *handle, const char *id)
{
  auto lmp = (LAMMPS *) handle;
  double **fexternal = nullptr;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    int tmp;
    fexternal = (double **)fix->extract("fexternal",tmp);
  }
  END_CAPTURE
  return fexternal;
}

/** Set the global energy contribution for a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback` and
:cpp:func:`lammps_fix_external_get_force` to also set the contribution
to the global energy from the external code.  The value of the *eng*
argument will be stored in the fix and applied on the current and all
following timesteps until changed by another call to this function.
The energy is in energy units as determined by the current :doc:`units <units>`
settings and is the **total** energy of the contribution.  Thus when
running in parallel all MPI processes have to call this function with
the **same** value and this will be returned as scalar property of the
fix external instance when accessed in LAMMPS input commands or from
variables.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  eng      total energy to be added to the global energy */

void lammps_fix_external_set_energy_global(void *handle, const char *id, double eng)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto fext = dynamic_cast<FixExternal*>(fix);
    fext->set_energy_global(eng);
  }
  END_CAPTURE
}

/** Set the global virial contribution for a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback`
and :cpp:func:`lammps_fix_external_get_force` to set the contribution to
the global virial from the external code.

The 6 values of the *virial* array will be stored in the fix and applied
on the current and all following timesteps until changed by another call
to this function. The components of the virial need to be stored in the
order: *xx*, *yy*, *zz*, *xy*, *xz*, *yz*.  In LAMMPS the virial is
stored internally as `stress*volume` in units of `pressure*volume` as
determined by the current :doc:`units <units>` settings and is the
**total** contribution.  Thus when running in parallel all MPI processes
have to call this function with the **same** value and this will then
be added by fix external.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  virial   the 6 global stress tensor components to be added to the global virial */

void lammps_fix_external_set_virial_global(void *handle, const char *id, double *virial)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto  fext = dynamic_cast<FixExternal*>(fix);
    fext->set_virial_global(virial);
  }
  END_CAPTURE
}

/** Set the per-atom energy contribution for a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback`
to set the per-atom energy contribution due to the fix from the external code
as part of the callback function.  For this to work, the handle to the
LAMMPS object must be passed as the *ptr* argument when registering the
callback function.

.. note::

   This function is fully independent from :cpp:func:`lammps_fix_external_set_energy_global`
   and will **NOT** add any contributions to the global energy tally
   and **NOT** check whether the sum of the contributions added here are
   consistent with the global added energy.


Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  eng      pointer to array of length nlocal with the energy to be added to the per-atom energy */

void lammps_fix_external_set_energy_peratom(void *handle, const char *id, double *eng)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto fext = dynamic_cast<FixExternal*>(fix);
    fext->set_energy_peratom(eng);
  }
  END_CAPTURE
}

/** Set the per-atom virial contribution for a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback`
to set the per-atom virial contribution due to the fix from the external code
as part of the callback function.  For this to work, the handle to the
LAMMPS object must be passed as the *ptr* argument when registering the
callback function.

.. note::

   This function is fully independent from :cpp:func:`lammps_fix_external_set_virial_global`
   and will **NOT** add any contributions to the global virial tally
   and **NOT** check whether the sum of the contributions added here are
   consistent with the global added virial.

The order and units of the per-atom stress tensor elements are the same
as for the global virial.  The code in fix external assumes the
dimensions of the per-atom virial array is ``double virial[nlocal][6]``.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  virial   a list of nlocal entries with the 6 per-atom stress tensor components to be added to the per-atom virial */

void lammps_fix_external_set_virial_peratom(void *handle, const char *id, double **virial)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto  fext = dynamic_cast<FixExternal*>(fix);
    fext->set_virial_peratom(virial);
  }
  END_CAPTURE
}

/** Set the vector length for a global vector stored with fix external for analysis

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback` and
:cpp:func:`lammps_fix_external_get_force` to set the length of a global vector of
properties that will be stored with the fix via
:cpp:func:`lammps_fix_external_set_vector`.

This function needs to be called **before** a call to
:cpp:func:`lammps_fix_external_set_vector` and **before** a run or minimize
command. When running in parallel it must be called from **all** MPI
processes and with the same length parameter.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  len      length of the global vector to be stored with the fix */

void lammps_fix_external_set_vector_length(void *handle, const char *id, int len)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto fext = dynamic_cast<FixExternal*>(fix);
    fext->set_vector_length(len);
  }
  END_CAPTURE
}

/** Store a global vector value for a fix external instance with the given ID.

\verbatim embed:rst

.. versionadded:: 28Jul2021

This is a companion function to :cpp:func:`lammps_set_fix_external_callback` and
:cpp:func:`lammps_fix_external_get_force` to set the values of a global vector of
properties that will be stored with the fix.  And can be accessed from
within LAMMPS input commands (e.g., fix ave/time or variables) when used
in a vector context.

This function needs to be called **after** a call to
:cpp:func:`lammps_fix_external_set_vector_length` and the  and **before** a run or minimize
command.  When running in parallel it must be called from **all** MPI
processes and with the **same** index and value parameters.  The value
is assumed to be extensive.

.. note::

   The index in the *idx* parameter is 1-based, i.e. the first element
   is set with idx = 1 and the last element of the vector with idx = N,
   where N is the value of the *len* parameter of the call to
   :cpp:func:`lammps_fix_external_set_vector_length`.

Please see the documentation for :doc:`fix external <fix_external>` for
more information about how to use the fix and how to couple it with an
external code.

\endverbatim
 *
 * \param  handle   pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param  id       fix ID of fix external instance
 * \param  idx      1-based index of in global vector
 * \param  val      value to be stored in global vector */

void lammps_fix_external_set_vector(void *handle, const char *id, int idx, double val)
{
  auto lmp = (LAMMPS *) handle;

  BEGIN_CAPTURE
  {
    auto fix = lmp->modify->get_fix_by_id(id);
    if (!fix) lmp->error->all(FLERR,"Can not find fix with ID '{}'!", id);

    if (strcmp("external",fix->style) != 0)
      lmp->error->all(FLERR,"Fix '{}' is not of style external!", id);

    auto  fext = dynamic_cast<FixExternal*>(fix);
    fext->set_vector(idx, val);
  }
  END_CAPTURE
}

/* ---------------------------------------------------------------------- */

/** Flush output buffers

\verbatim embed:rst
This function can be used to flush buffered output to be written to screen
and logfile pointers to simplify capturing output from LAMMPS library calls.
\endverbatim
 *
 * \param  handle    pointer to a previously created LAMMPS instance cast to ``void *``.
 */
void lammps_flush_buffers(void *handle) {
  utils::flush_buffers((LAMMPS *) handle);
}

/* ---------------------------------------------------------------------- */

/** Free memory buffer allocated by LAMMPS.
 *
\verbatim embed:rst

Some of the LAMMPS C library interface functions return data as pointer
to a buffer that has been allocated by LAMMPS or the library interface.
This function can be used to delete those in order to avoid memory
leaks.

\endverbatim
 *
 * \param  ptr  pointer to data allocated by LAMMPS */

void lammps_free(void *ptr)
{
  free(ptr);
}

/* ---------------------------------------------------------------------- */

/** Check if LAMMPS is currently inside a run or minimization
 *
 * This function can be used from signal handlers or multi-threaded
 * applications to determine if the LAMMPS instance is currently active.
 *
 * \param  handle pointer to a previously created LAMMPS instance cast to ``void *``.
 * \return        0 if idle or >0 if active */

int lammps_is_running(void *handle)
{
  auto   lmp = (LAMMPS *) handle;
  return lmp->update->whichflag;
}

/** Force a timeout to stop an ongoing run cleanly.
 *
 * This function can be used from signal handlers or multi-threaded
 * applications to cleanly terminate an ongoing run.
 *
 * \param  handle pointer to a previously created LAMMPS instance cast to ``void *`` */

void lammps_force_timeout(void *handle)
{
  auto   lmp = (LAMMPS *) handle;
  return lmp->timer->force_timeout();
}

// ----------------------------------------------------------------------
// Library functions for error handling with exceptions enabled
// ----------------------------------------------------------------------

/** Check if there is a (new) error message available

\verbatim embed:rst
This function can be used to query if an error inside of LAMMPS
has thrown a :ref:`C++ exception <exceptions>`.

.. note::

   .. versionchanged: 2Aug2023

   The *handle* pointer may be ``NULL`` for this function, as would be
   the case when a call to create a LAMMPS instance has failed.  Then
   this function will not check the error status inside the LAMMPS
   instance, but instead would check the global error buffer of the
   library interface.

.. note::

   This function will always report "no error" when the LAMMPS library
   has been compiled without ``-DLAMMPS_EXCEPTIONS``, which turns fatal
   errors aborting LAMMPS into C++ exceptions. You can use the library
   function :cpp:func:`lammps_config_has_exceptions` to check whether this is
   the case.
\endverbatim
 *
 * \param handle   pointer to a previously created LAMMPS instance cast to ``void *`` or NULL
 * \return 0 on no error, 1 on error.
 */
int lammps_has_error(void *handle)
{
  if (handle) {
    LAMMPS *lmp = (LAMMPS *) handle;
    Error *error = lmp->error;
    return (error->get_last_error().empty()) ? 0 : 1;
  } else {
    return lammps_last_global_errormessage.empty() ? 0 : 1;
  }
}

/* ---------------------------------------------------------------------- */

/** Copy the last error message into the provided buffer

\verbatim embed:rst
This function can be used to retrieve the error message that was set
in the event of an error inside of LAMMPS which resulted in a
:ref:`C++ exception <exceptions>`.  A suitable buffer for a C-style
string has to be provided and its length.  If the internally stored
error message is longer, it will be truncated accordingly.  The return
value of the function corresponds to the kind of error: a "1" indicates
an error that occurred on all MPI ranks and is often recoverable, while
a "2" indicates an abort that would happen only in a single MPI rank
and thus may not be recoverable, as other MPI ranks may be waiting on
the failing MPI ranks to send messages.

.. note::

   .. versionchanged: 2Aug2023

   The *handle* pointer may be ``NULL`` for this function, as would be
   the case when a call to create a LAMMPS instance has failed.  Then
   this function will not check the error buffer inside the LAMMPS
   instance, but instead would check the global error buffer of the
   library interface.

.. note::

   This function will do nothing when the LAMMPS library has been
   compiled without ``-DLAMMPS_EXCEPTIONS``, which turns errors aborting
   LAMMPS into C++ exceptions.  You can use the library function
   :cpp:func:`lammps_config_has_exceptions` to check whether this is the case.
\endverbatim
 *
 * \param  handle    pointer to a previously created LAMMPS instance cast to ``void *`` or NULL.
 * \param  buffer    string buffer to copy the error message to
 * \param  buf_size  size of the provided string buffer
 * \return           1 when all ranks had the error, 2 on a single rank error. */

int lammps_get_last_error_message(void *handle, char *buffer, int buf_size)
{
  if (handle) {
    LAMMPS *lmp = (LAMMPS *) handle;
    Error *error = lmp->error;
    buffer[0] = buffer[buf_size-1] = '\0';

    if (!error->get_last_error().empty()) {
      int error_type = error->get_last_error_type();
      strncpy(buffer, error->get_last_error().c_str(), buf_size-1);
      error->set_last_error("", ERROR_NONE);
      return error_type;
    }
  } else {
    buffer[0] = buffer[buf_size-1] = '\0';

    if (!lammps_last_global_errormessage.empty()) {
      strncpy(buffer, lammps_last_global_errormessage.c_str(), buf_size-1);
      lammps_last_global_errormessage.clear();
      return 1;
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

/** Return API version of embedded Python interpreter

\verbatim embed:rst

.. versionadded:: 3Nov2022

This function is used by the ML-IAP python code (mliappy) to verify
the API version of the embedded python interpreter of the PYTHON
package.  It returns -1 if the PYTHON package is not enabled.

\endverbatim
 *
 * \return   PYTHON_API_VERSION constant of the python interpreter or -1 */

int lammps_python_api_version() {
#if defined(LMP_PYTHON)
  return PYTHON_API_VERSION;
#else
  return -1;
#endif
}

// Local Variables:
// fill-column: 72
// End:
