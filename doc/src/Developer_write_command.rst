Writing a new command style
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Command styles allow to do system manipulations or interfaces to the
operating system.

In the text below, we will discuss the implementation of one example.  As
shown on the page for :doc:`writing or extending command styles
<Modify_command>`, in order to implement a new command style, a new class
must be written that is either directly or indirectly derived from the
``Command`` class.  There is just one method that must be implemented:
``Command::command()``.  In addition, a custom constructor is needed to get
access to the members of the ``LAMMPS`` class like the ``Error`` class to
print out error messages.  The ``Command::command()`` method processes the
arguments passed to the command in the input and executes it.  Any other
methods would be for the convenience of implementation of the new command.

In general, new command styles should be added to the :ref:`EXTRA-COMMAND
package <PKG-EXTRA-COMMAND>`.  If you feel that your contribution should be
added to a different package, please consult with the :doc:`LAMMPS
developers <Intro_authors>` first.  The contributed code needs to support
the :doc:`traditional GNU make build process <Build_make>` **and** the
:doc:`CMake build process <Build_cmake>`.

----

Case 1: Implementing the geturl command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, we will describe the procedure of adding a simple command
style to LAMMPS: the :doc:`geturl command <geturl>` that allows to download
files directly without having to rely on an external program like "wget" or
"curl".  The complete implementation can be found in the files
``src/EXTRA-COMMAND/geturl.cpp`` and ``src/EXTRA-COMMAND/geturl.h`` of the
LAMMPS source code.

Interfacing the *libcurl* library
"""""""""""""""""""""""""""""""""

Rather than implementing the various protocols for downloading files, we
rely on an external library: `libcurl library <https:://curl.se/libcurl/>`_.
This requires that the library and its headers are installed.  For the
traditional GNU make build system, this simply requires edits to the machine
makefile to add compilation flags like for other libraries.  For the CMake
based build system, we need to add some lines to the file
``cmake/Modules/Packages/EXTRA-COMMAND.cmake``:

.. code-block:: cmake

   find_package(CURL QUIET COMPONENTS HTTP HTTPS)
   option(WITH_CURL "Enable libcurl support" ${CURL_FOUND})
   if(WITH_CURL)
     find_package(CURL REQUIRED COMPONENTS HTTP HTTPS)
     target_compile_definitions(lammps PRIVATE -DLAMMPS_CURL)
     target_link_libraries(lammps PRIVATE CURL::libcurl)
   endif()

The first ``find_package()`` command uses a built-in CMake module to find
an existing *libcurl* installation with development headers and support for
using the HTTP and HTTPS protocols.  The "QUIET" flag ensures that there is
no screen output and no error if the search fails.  The status of the search
is recorded in the "${CURL_FOUND}" variable.  That variable sets the default
of the WITH_CURL option, which toggles whether support for *libcurl* is included
or not.

The second ``find_package()`` uses the "REQUIRED" flag to produce an error
if the WITH_CURL option was set to ``True``, but no suitable *libcurl*
implementation with development support was found.  This construct is used
so that the CMake script code inside the ``if(WITH_CURL)`` and ``endif()``
block can be expanded later to download and compile *libcurl* as part of the
LAMMPS build process, if it is not found locally.  The
``target_compile_definitions()`` function added the define ``-DLAMMPS_CURL``
to the compilation flags when compiling objects for the LAMMPS library.
This allows to always compile the :doc:`geturl command <geturl>`, but use
pre-processing to compile in the interface to *libcurl* only when it is
present and usable and otherwise stop with an error message about the
unavailability of *libcurl* to execute the functionality of the command.

Header file
"""""""""""

The first segment of any LAMMPS source should be the copyright and
license statement.  Note the marker in the first line to indicate to
editors like emacs that this file is a C++ source, even though the .h
extension suggests a C source (this is a convention inherited from the
very beginning of the C++ version of LAMMPS).

.. code-block:: c++

   /* -*- c++ -*- ----------------------------------------------------------
      LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
      https://www.lammps.org/, Sandia National Laboratories
      LAMMPS development team: developers@lammps.org

      Copyright (2003) Sandia Corporation.  Under the terms of Contract
      DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
      certain rights in this software.  This software is distributed under
      the GNU General Public License.

      See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

Every command style must be registered in LAMMPS by including the following
lines of code in the second part of the header after the copyright
message and before the include guards for the class definition:

.. code-block:: c++

   #ifdef COMMAND_CLASS
   // clang-format off
   CommandStyle(geturl,GetURL);
   // clang-format on
   #else

This block between ``#ifdef COMMAND_CLASS`` and ``#else`` will be
included by the ``Input`` class in ``input.cpp`` to build a map of
"factory functions" that will create an instance of a Command class
and call its ``command()`` method.  The map connects the name of the
command ``geturl`` with the name of the class ``GetURL``.  During
compilation, LAMMPS constructs a file ``style_command.h`` that contains
``#include`` statements for all "installed" command styles.  Before
including ``style_command.h`` into ``input.cpp``, the ``COMMAND_CLASS``
define is set and the ``CommandStyle(name,class)`` macro defined.  The
code of the macro adds the installed command styles to the "factory map"
which enables the ``Input`` to execute the command.

The list of header files to include in ``style_command.h`` is automatically
updated by the build system if there are new files, so the presence of the
new header file in the ``src/EXTRA-COMMAND`` folder and the enabling of the
EXTRA-COMMAND package will trigger LAMMPS to include the new command style
when it is (re-)compiled.  The "// clang-format" format comments are needed
so that running :ref:`clang-format <clang-format>` on the file will not
insert unwanted blanks which would break the ``CommandStyle`` macro.

The third part of the header file is the actual class definition of the
``GetURL`` class.  This has the custom constructor and the ``command()``
method implemented by this command style.  For the constructor there is
nothing to do but to pass the ``lmp`` pointer to the base class.  Since the
``command()`` method is labeled "virtual" in the base class, it must be
given the "override" property.

.. code-block:: c++

   #ifndef LMP_GETURL_H
   #define LMP_GETURL_H

   #include "command.h"

   namespace LAMMPS_NS {

   class GetURL : public Command {
    public:
     GetURL(class LAMMPS *lmp) : Command(lmp) {};
     void command(int, char **) override;
   };
   }    // namespace LAMMPS_NS
   #endif
   #endif

The "override" property helps to detect unexpected mismatches because
compilation will stop with an error in case the signature of a function
is changed in the base class without also changing it in all derived
classes.

Implementation file
"""""""""""""""""""

We move on to the implementation of the ``GetURL`` class in the
``geturl.cpp`` file.  This file also starts with a LAMMPS copyright and
license header.  Below that notice is typically the space where comments may
be added with additional information about this specific file, the
author(s), affiliation(s), and email address(es).  This way the contributing
author(s) can be easily contacted, when there are questions about the
implementation later.  Since the file(s) may be around for a long time, it
is beneficial to use some kind of "permanent" email address, if possible.

.. code-block:: c++

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
      Contributing authors:  Axel Kohlmeyer (Temple U),
   ------------------------------------------------------------------------- */

   #include "geturl.h"

   #include "comm.h"
   #include "error.h"

   #if defined(LAMMPS_CURL)
   #include <curl/curl.h>
   #endif

   using namespace LAMMPS_NS;

The second section of the implementation file has various include
statements.  The include file for the class header has to come first, then a
couple of LAMMPS classes (sorted alphabetically) followed by the header for
the *libcurl* interface.  This is wrapped into an ``#ifdef`` block so that
LAMMPS will compile this file without error when the *libcurl* header is not
available and thus the define not set.  The final statement of this segment
imports the ``LAMMPS_NS::`` namespace globally for this file.  This way, all
LAMMPS specific functions and classes do not have to be prefixed with
``LAMMPS_NS::``.

The command() function (required)
"""""""""""""""""""""""""""""""""

Since the required custom constructor is trivial and implemented in the
header, there is only one function that must be implemented for a command
style and that is the ``command()`` function.

.. code-block:: c++

   void GetURL::command(int narg, char **arg)
   {
   #if !defined(LAMMPS_CURL)
     error->all(FLERR, "LAMMPS has not been compiled with libcurl support");
   #else
     if (narg < 1) utils::missing_cmd_args(FLERR, "geturl", error);
     int verify = 1;
     int overwrite = 1;
     int verbose = 0;

This first part also has the ``#ifdef`` block depending on the LAMMPS_CURL
define.  This way the command will simply print an error, if *libcurl* is
not available but will not fail to compile.  Furthermore, it sets the
defaults for the following optional arguments.

.. code-block:: c++

     // process arguments

     std::string url = arg[0];

     // sanity check

     if ((url.find(':') == std::string::npos) || (url.find('/') == std::string::npos))
       error->all(FLERR, "URL '{}' is not a supported URL", url);

     std::string output = url.substr(url.find_last_of('/') + 1);
     if (output.empty()) error->all(FLERR, "URL '{}' must end in a file string", url);

This block stores the positional, i.e. non-optional argument of the URL to
be downloaded and adds a couple of sanity checks on the string to make sure it is
a valid URL.  Also it derives the default name of the output file from the URL.

.. code-block:: c++

     int iarg = 1;
     while (iarg < narg) {
       if (strcmp(arg[iarg], "output") == 0) {
         if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl output", error);
         output = arg[iarg + 1];
         ++iarg;
       } else if (strcmp(arg[iarg], "overwrite") == 0) {
         if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl overwrite", error);
         overwrite = utils::logical(FLERR, arg[iarg + 1], false, lmp);
         ++iarg;
       } else if (strcmp(arg[iarg], "verify") == 0) {
         if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl verify", error);
         verify = utils::logical(FLERR, arg[iarg + 1], false, lmp);
         ++iarg;
       } else if (strcmp(arg[iarg], "verbose") == 0) {
         if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl verbose", error);
         verbose = utils::logical(FLERR, arg[iarg + 1], false, lmp);
         ++iarg;
       } else {
         error->all(FLERR, "Unknown geturl keyword: {}", arg[iarg]);
       }
       ++iarg;
     }

This block parses the optional arguments following the URL and stops with an
error if there are arguments missing or an unknown argument is encountered.

.. code-block:: c++

     // only download files from rank 0

     if (comm->me != 0) return;

     if (!overwrite && platform::file_is_readable(output)) return;

     // open output file for writing

     FILE *out = fopen(output.c_str(), "wb");
     if (!out)
       error->all(FLERR, "Cannot open output file {} for writing: {}", output, utils::getsyserror());

Here all MPI ranks other than 0 will return, so that the URL download will
only happen from a single MPI rank. For that rank the output file is opened
for writing using the C library function ``fopen()``.

.. code-block:: c++

     // initialize curl and perform download

     CURL *curl;
     curl_global_init(CURL_GLOBAL_DEFAULT);
     curl = curl_easy_init();
     if (curl) {
       (void) curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
       (void) curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) out);
       (void) curl_easy_setopt(curl, CURLOPT_FILETIME, 1L);
       (void) curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
       if (verbose && screen) {
         (void) curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
         (void) curl_easy_setopt(curl, CURLOPT_STDERR, (void *) screen);
       }
       if (!verify) {
         (void) curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
         (void) curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);
       }
       auto res = curl_easy_perform(curl);
       if (res != CURLE_OK) {
         long response = 0L;
         curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response);
         error->one(FLERR, "Download of {} failed with: {} {}", output, curl_easy_strerror(res),
                    response);
       }
       curl_easy_cleanup(curl);

This block now implements the actual URL download with the selected options
via the "easy" interface of *libcurl*.  For the details of what these
function calls do, please have a look at the `*libcurl documentation
<https://curl.se/libcurl/c/allfuncs.html>`_.

 .. code-block:: c++

     }
     curl_global_cleanup();
     fclose(out);
   #endif
   }

Finally, the previously opened file is closed and the command is complete.
