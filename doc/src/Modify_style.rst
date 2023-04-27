LAMMPS programming style and requirements for contributions
===========================================================

The following is a summary of the current requirements and
recommendations for including contributed source code or documentation
into the LAMMPS software distribution.

Motivation
----------

The LAMMPS developers are committed to providing a software package that
is versatile, reliable, high-quality, efficient, portable, and easy to
maintain and modify.  Achieving all of these goals is challenging since
a large part of LAMMPS consists of contributed code from many different
authors and not many of them are professionally trained programmers and
familiar with the idiosyncrasies of maintaining a large software
package.  In addition, changes that interfere with the parallel
efficiency of the core code must be avoided.  As LAMMPS continues to
grow and more features and functionality are added, it becomes a
necessity to be more discriminating with new contributions while also
working at the same time to improve the existing code.

The following requirements and recommendations are provided to help
maintaining or improving that status.  Where possible we utilize
available continuous integration tools to search for common programming
mistakes, portability limitations, incompatible formatting, and
undesired side effects.  It is indicated which requirements are strict,
and which represent a preference and thus are negotiable or optional.

Please feel free to contact the LAMMPS core developers in case you need
additional explanations or clarifications or in case you need assistance
in realizing the (strict) requirements for your contributions.

Licensing requirements (strict)
-------------------------------

Contributing authors agree when submitting a pull request that their
contributions can be distributed under the LAMMPS license
conditions. This is the GNU public license in version 2 (not 3 or later)
for the publicly distributed versions, e.g. on the LAMMPS homepage or on
GitHub.  On request we also make a version of LAMMPS available under
LGPL 2.1 terms; this will usually be the latest available or a previous
stable version with a few LGPL 2.1 incompatible files removed. More details
are found on the :doc:`LAMMPS open-source license page <Intro_opensource>`.

Your new source files should have the LAMMPS copyright, GPL notice, and
your name and email address at the top, like other user-contributed
LAMMPS source files.

Contributions may be under a different license as long as that
license does not conflict with the aforementioned terms.  Contributions
that use code with a conflicting license can be split into two parts:

1. the core parts (i.e. parts that must be in the `src` tree) that are
   licensed under compatible terms and bundled with the LAMMPS sources
2. an external library that must be downloaded and compiled (either
   separately or as part of the LAMMPS compilation)

Please note, that this split licensed mode may complicate including the
contribution in binary packages.

Using pull requests on GitHub (preferred)
-----------------------------------------

All contributions to LAMMPS are processed as pull requests on GitHub
(this also applies to the work of the core LAMMPS developers).  A
:doc:`tutorial for submitting pull requests on GitHub <Howto_github>` is
provided.  If this is still problematic, contributors may contact any of
the core LAMMPS developers for help or to create a pull request on their
behalf.  This latter way of submission may delay the integration as it
depends on the amount of time required to prepare the pull request and
free time available by the LAMMPS developer in question to spend on this
task.

Integration testing (strict)
----------------------------

Contributed code, like all pull requests, must pass the automated
tests on GitHub before it can be merged with the LAMMPS distribution.
These tests compile LAMMPS in a variety of environments and settings and
run the bundled unit tests.  At the discretion of the LAMMPS developer
managing the pull request, additional tests may be activated that test
for "side effects" on running a collection of input decks and create
consistent results.  Also, the translation of the documentation to HTML
and PDF is tested for.

More specifically, this means that contributed source code **must**
compile with the most current version of LAMMPS with ``-DLAMMPS_BIGBIG``
in addition to the default setting of ``-DLAMMPS_SMALLBIG``.  The code
needs to work correctly in both cases and also in serial and parallel
using MPI.

Some "disruptive" changes may break tests and require updates to the
testing tools or scripts or tests themselves.  This is rare.  If in
doubt, contact the LAMMPS developer that is assigned to the pull request
for further details and explanations and suggestions of what needs to be
done.

Documentation (strict)
----------------------

Contributions that add new styles or commands or augment existing ones
must include the corresponding new or modified documentation in
`ReStructuredText format <rst_>`_ (.rst files in the ``doc/src/``
folder). The documentation shall be written in American English and the
.rst file must use only ASCII characters so it can be cleanly translated
to PDF files (via `sphinx <https://www.sphinx-doc.org>`_ and PDFLaTeX).
Special characters may be included via embedded math expression typeset
in a LaTeX subset.

.. _rst: https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html

When adding new commands, they need to be integrated into the sphinx
documentation system, and the corresponding command tables and lists
updated. When translating the documentation into html files there should
be no warnings. When adding a new package also some lists describing
packages must be updated as well as a package specific description added
and, if necessary, some package specific build instructions included.

As appropriate, the text files with the documentation can include inline
mathematical expression or figures (see ``doc/JPG`` for examples).
Additional PDF files with further details (see ``doc/PDF`` for examples) may
also be included.  The page should also include literature citations as
appropriate; see the bottom of ``doc/fix_nh.rst`` for examples and the
earlier part of the same file for how to format the cite itself.
Citation labels must be unique across **all** .rst files.  The
"Restrictions" section of the page should indicate if your command is
only available if LAMMPS is built with the appropriate FOO package.  See
other package doc files for examples of how to do this.

Please run at least "make html" and "make spelling" and carefully
inspect and proofread the resulting HTML format doc page before
submitting your code.  Upon submission of a pull request, checks for
error free completion of the HTML and PDF build will be performed and
also a spell check, a check for correct anchors and labels, and a check
for completeness of references all styles in their corresponding tables
and lists is run.  In case the spell check reports false positives they
can be added to the file ``doc/utils/sphinx-config/false_positives.txt``

Contributions that add or modify the library interface or "public" APIs
from the C++ code or the Fortran module must include suitable doxygen
comments in the source and corresponding changes to the documentation
sources for the "Programmer Guide" guide section of the LAMMPS manual.

Examples (preferred)
--------------------

In most cases, it is preferred that example scripts (simple, small, fast
to complete on 1 CPU) are included that demonstrate the use of new or
extended functionality. These are typically under the examples or
examples/PACKAGES directory.  A few guidelines for such example input
decks.

- commands that generate output should be commented out (except when the
  output is the sole purpose or the feature, e.g. for a new compute).

- commands like :doc:`log <log>`, :doc:`echo <echo>`, :doc:`package
  <package>`, :doc:`processors <processors>`, :doc:`suffix <suffix>` may
  **not** be used in the input file (exception: "processors * * 1" or
  similar is acceptable when used to avoid unwanted domain decomposition
  of empty volumes).

- outside of the log files, no generated output should be included

- custom thermo_style settings may not include output measuring CPU or other time
  as that makes comparing the thermo output between different runs more complicated.

- input files should be named ``in.name``, data files should be named
  ``data.name`` and log files should be named ``log.version.name.<compiler>.<ncpu>``

- the total file size of all the inputs and outputs should be small

- where possible, potential files from the "potentials" folder or data
  file from other folders should be re-used through symbolic links

Howto document (optional)
-------------------------

If your feature requires some more complex steps and explanations to be
used correctly or some external or bundled tools or scripts, we
recommend that you also contribute a :doc:`Howto document <Howto>`
providing some more background information and some tutorial material.
This can also be used to provide more in-depth explanations for bundled
examples.

As a general rule-of-thumb, the more clear and self-explanatory you make
your documentation, README files and examples, and the easier you make
it for people to get started, the more likely it is that users will try
out your new feature.

Programming style requirements (varied)
---------------------------------------

The LAMMPS developers aim to employ a consistent programming style and
naming conventions across the entire code base, as this helps with
maintenance, debugging, and understanding the code, both for developers
and users.

The files `pair_lj_cut.h`, `pair_lj_cut.cpp`, `utils.h`, and `utils.cpp`
may serve as representative examples.

Command or Style names, file names, and keywords (mostly strict)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All user-visible command or style names should be all lower case and
should only use letters, numbers, or forward slashes.  They should be
descriptive and initialisms should be avoided unless they are well
established (e.g. lj for Lennard-Jones).  For a compute style
"some/name" the source files must be called `compute_some_name.h` and
`compute_some_name.cpp`. The "include guard" would then be
`LMP_COMPUTE_SOME_NAME_H` and the class name `ComputeSomeName`.

Whitespace and permissions (preferred)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Source files should not contain TAB characters unless required by the
syntax (e.g. in makefiles) and no trailing whitespace.  Text files
should be added with Unix-style line endings (LF-only). Git will
automatically convert those in both directions when running on Windows;
use dos2unix on Linux machines to convert files.  Text files should have
a line ending on the last line.

All files should have 0644 permissions, i.e writable to the user only
and readable by all and no executable permissions.  Executable
permissions (0755) should only be on shell scripts or python or similar
scripts for interpreted script languages.

You can check for these issues with the python scripts in the
:ref:`"tools/coding_standard" <coding_standard>` folder.  When run
normally with a source file or a source folder as argument, they will
list all non-conforming lines.  By adding the `-f` flag to the command
line, they will modify the flagged files to try removing the detected
issues.

Indentation and placement of braces (strongly preferred)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS uses 2 characters per indentation level and lines should be
kept within 100 characters wide.

For new files added to the "src" tree, a `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_ configuration file is
provided under the name `.clang-format`.  This file is compatible with
clang-format version 8 and later. With that file present, files can be
reformatted according to the configuration with a command like:
`clang-format -i new-file.cpp`.  Ideally, this is done while writing the
code or before a pull request is submitted.  Blocks of code where the
reformatting from clang-format yields undesirable output may be
protected with placing a pair `// clang-format off` and `// clang-format
on` comments around that block.

Error or warning messages and explanations (preferred)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionchanged:: 4May2022

Starting with LAMMPS version 4 May 2022 the LAMMPS developers have
agreed on a new policy for error and warning messages.

Previously, all error and warning strings were supposed to be listed in
the class header files with an explanation.  Those would then be
regularly "harvested" and transferred to alphabetically sorted lists in
the manual.  To avoid excessively long lists and to reduce effort, this
came with a requirement to have rather generic error messages (e.g.
"Illegal ... command").  To identify the specific cause, the name of the
source file and the line number of the error location would be printed,
so that one could look up the cause by reading the source code.

The new policy encourages more specific error messages that ideally
indicate the cause directly and no further lookup would be needed.
This is aided by using the `{fmt} library <https://fmt.dev>`_ to convert
the Error class commands so that they take a variable number of arguments
and error text will be treated like a {fmt} syntax format string.
Error messages should still kept to a single line or two lines at the most.

For more complex explanations or errors that have multiple possible
reasons, a paragraph should be added to the `Error_details` page with an
error code reference (e.g. ``.. _err0001:``) then the utility function
:cpp:func:`utils::errorurl() <LAMMPS_NS::utils::errorurl>` can be used
to generate an URL that will directly lead to that paragraph.  An error
for missing arguments can be easily generated using the
:cpp:func:`utils::missing_cmd_args()
<LAMMPS_NS::utils::missing_cmd_args>` convenience function.

The transformation of existing LAMMPS code to this new scheme is ongoing
and - given the size of the LAMMPS source code - will take a significant
amount of time until completion.  However, for new code following the
new approach is strongly preferred.  The expectation is that the new
scheme will make it easier for LAMMPS users, developers, and
maintainers.

An example for this approach would be the
``src/read_data.cpp`` and ``src/atom.cpp`` files that implement the
:doc:`read_data <read_data>` and :doc:`atom_modify <atom_modify>`
commands and that may create :ref:`"Unknown identifier in data file" <err0001>`
errors that seem difficult to debug for users because they may have
one of multiple possible reasons, and thus require some additional explanations.

Programming language standards (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The core of LAMMPS is written in C++11 in a style that can be mostly
described as "C with classes".  Advanced C++ features like operator
overloading or excessive use of templates are avoided with the intent to
keep the code readable to programmers that have limited C++ programming
experience.  C++ constructs are acceptable when they help improve the
readability and reliability of the code, e.g. when using the
`std::string` class instead of manipulating pointers and calling the
string functions of the C library.  In addition a collection of
convenient :doc:`utility functions and classes <Developer_utils>` for
recurring tasks and a collection of
:doc:`platform neutral functions <Developer_platform>` for improved
portability are provided.

Included Fortran code has to be compatible with the Fortran 2003
standard.  Python code must be compatible with Python 3.5.  Large parts
of LAMMPS (including the :ref:`PYTHON package <PKG-PYTHON>`) are also
compatible with Python 2.7.  Compatibility with Python 2.7 is
desirable, but compatibility with Python 3.5 is **required**.

Compatibility with these older programming language standards is very
important to maintain portability and availability of LAMMPS on many
platforms.  This applies especially to HPC cluster environments, which
tend to be running older software stacks and LAMMPS users may be
required to use those older tools for access to advanced hardware
features or not have the option to install newer compilers or libraries.

Programming conventions (varied)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following is a collection of conventions that should be applied when
writing code for LAMMPS.  Following these steps will make it much easier
to integrate your contribution. Please have a look at the existing files
in packages in the src directory for examples.  As a demonstration for
how can be adapted to these conventions you may compare the REAXFF
package with the what it looked like when it was called USER-REAXC.  If
you are uncertain, please ask.

- system headers or from installed libraries are include with angular
  brackets (example: ``#include <vector>``), while local include file
  use double quotes (example: ``#include "atom.h"``).

- when including system header files from the C library use the
    C++-style names (``<cstdlib>`` or ``<cstring>``) instead of the
    C-style names (``<stdlib.h>`` or ``<string.h>``)

- the order of ``#include`` statements in a file ``some_name.cpp`` that
  implements a class ``SomeName`` defined in a header file
  ``some_name.h`` should be as follows:

  - ``#include "some_name.h"`` followed by an empty line

  - LAMMPS include files e.g. ``#include "comm.h"`` or ``#include
    "modify.h"`` in alphabetical order followed by an empty line

  - system header files from the C++ or C standard library followed by
    an empty line

  - ``using namespace LAMMPS_NS`` or other namespace imports.

- I/O is done via the C-style stdio library and **not** iostreams.

- Do not use so-called "alternative tokens" like ``and``, ``or``,
  ``not`` and similar, but rather use the corresponding operators
  ``&&``, ``||``, and ``!``.  The alternative tokens are not available
  by default on all compilers, and also we want to maintain a consistent
  programming style.

- Output to the screen and the logfile should be using the corresponding
  FILE pointers and only be done on MPI rank 0.  Use the :cpp:func:`utils::logmesg`
  convenience function where possible.

- Usage of C++11 `virtual`, `override`, `final` keywords: Please follow the
  `C++ Core Guideline C.128 <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rh-override>`_.
  That means, you should only use `virtual` to declare a new virtual
  function, `override` to indicate you are overriding an existing virtual
  function, and `final` to prevent any further overriding.

- Trivial destructors: Prefer not writing destructors when they are empty and `default`.

  .. code-block:: c++

     // don't write destructors for A or B like this
     class A : protected Pointers {
      public:
        A();
        ~A() override {}
     };

     class B : protected Pointers {
      public:
        B();
        ~B() override = default;
     };

     // instead, let the compiler create the implicit default destructor by not writing it
     class A : protected Pointers {
      public:
        A();
     };

     class B : protected Pointers {
      public:
        B();
     };

- Header files, especially those defining a "style", should only use
  the absolute minimum number of include files and **must not** contain
  any ``using`` statements. Typically that would be only the header for
  the base class. Instead any include statements should be put into the
  corresponding implementation files and forward declarations be used.
  For implementation files, the "include what you use" principle should
  be employed.  However, there is the notable exception that when the
  ``pointers.h`` header is included (or one of the base classes derived
  from it) certain headers will always be included and thus do not need
  to be explicitly specified.
  These are: `mpi.h`, `cstddef`, `cstdio`, `cstdlib`, `string`, `utils.h`,
  `vector`, `fmt/format.h`, `climits`, `cinttypes`.
  This also means any such file can assume that `FILE`, `NULL`, and
  `INT_MAX` are defined.

- Header files that define a new LAMMPS style (i.e. that have a
  ``SomeStyle(some/name,SomeName);`` macro in them) should only use the
  include file for the base class and otherwise use forward declarations
  and pointers; when interfacing to a library use the PIMPL (pointer
  to implementation) approach where you have a pointer to a struct
  that contains all library specific data (and thus requires the library
  header) but use a forward declaration and define the struct only in
  the implementation file. This is a **strict** requirement since this
  is where type clashes between packages and hard to find bugs have
  regularly manifested in the past.

- Please use clang-format only to reformat files that you have
  contributed.  For header files containing a ``SomeStyle(keyword,
  ClassName)`` macros it is required to have this macro embedded with a
  pair of ``// clang-format off``, ``// clang-format on`` comments and
  the line must be terminated with a semi-colon (;).  Example:

  .. code-block:: c++

     #ifdef COMMAND_CLASS
     // clang-format off
     CommandStyle(run,Run);
     // clang-format on
     #else

     #ifndef LMP_RUN_H
     [...]

  You may also use ``// clang-format on/off`` throughout your files
  to protect individual sections from being reformatted.

- We rarely accept new styles in the core src folder.  Thus please
  review the list of :doc:`available Packages <Packages_details>` to see
  if your contribution could be added to be added to one of them.  It
  should fit into the general purposed of that package.  If it does not
  fit well, it may be added to one of the EXTRA- packages or the MISC
  package.


Contributing a package
----------------------

If your contribution has several related features that are not covered
by one of the existing packages or is dependent on a library (bundled or
external), it is best to make it a package directory with a name like
FOO.  In addition to your new files, the directory should contain a
README text file.  The README should contain your name and contact
information and a brief description of what your new package does.


Build system (strongly preferred)
---------------------------------

LAMMPS currently supports two build systems: one that is based on
:doc:`traditional Makefiles <Build_make>` and one that is based on
:doc:`CMake <Build_cmake>`.  Thus your contribution must be compatible
with and support both.

For a single pair of header and implementation files that are an
independent feature, it is usually only required to add them to
`src/.gitignore``.

For traditional make, if your contributed files or package depend on
other LAMMPS style files or packages also being installed (e.g. because
your file is a derived class from the other LAMMPS class), then an
Install.sh file is also needed to check for those dependencies and
modifications to src/Depend.sh to trigger the checks.  See other README
and Install.sh files in other directories as examples.

Similarly for CMake support, changes may need to be made to
cmake/CMakeLists.txt, some of the files in cmake/presets, and possibly a
file with specific instructions needs to be added to
cmake/Modules/Packages/.  Please check out how this is handled for
existing packages and ask the LAMMPS developers if you need assistance.


Citation reminder (suggested)
-----------------------------

If there is a paper of yours describing your feature (either the
algorithm/science behind the feature itself, or its initial usage, or
its implementation in LAMMPS), you can add the citation to the \*.cpp
source file.  See ``src/DIFFRACTION/compute_saed.cpp`` for an example.
A BibTeX format citation is stored in a string variable at the top
of the file and  a single line of code registering this variable is
added to the constructor of the class.  When your feature is used,
by default, LAMMPS will print the brief info and the DOI
in the first line to the screen and the full citation to the log file.

If there is additional functionality (which may have been added later)
described in a different publication, additional citation descriptions
may be added for as long as they are only registered when the
corresponding keyword activating this functionality is used.  With these
options it is possible to have LAMMPS output a specific citation
reminder whenever a user invokes your feature from their input script.
Please note that you should *only* use this for the *most* relevant
paper for a feature and a publication that you or your group authored.
E.g. adding a citation in the code for a paper by Nose and Hoover if you
write a fix that implements their integrator is not the intended usage.
That latter kind of citation should just be included in the
documentation page you provide describing your contribution.  If you are
not sure what the best option would be, please contact the LAMMPS
developers for advice.


Testing (optional)
------------------

If your contribution contains new utility functions or a supporting class
(i.e. anything that does not depend on a LAMMPS object), new unit tests
should be added to a suitable folder in the ``unittest`` tree.
When adding a new LAMMPS style computing forces or selected fixes,
a ``.yaml`` file with a test configuration and reference data should be
added for the styles where a suitable tester program already exists
(e.g. pair styles, bond styles, etc.). Please see
:ref:`this section in the manual <testing>` for more information on
how to enable, run, and expand testing.
