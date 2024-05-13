Requirements for contributions to LAMMPS
========================================

The following is a summary of the current requirements and
recommendations for including contributed source code or documentation
into the LAMMPS software distribution.

Motivation
----------

The LAMMPS developers are committed to provide a software package that
is versatile, reliable, high-quality, efficient, portable, and easy to
maintain and modify.  Achieving all of these goals is challenging
since a large part of LAMMPS consists of contributed code from many
different authors who may not be professionally trained programmers or
familiar with the idiosyncrasies of maintaining a large software
package.  In addition, changes that interfere with the parallel
efficiency of the core code must be avoided.  As LAMMPS continues to
grow and more features and functionality are added, it is necessary to
follow established guidelines when accepting new contributions while
also working at the same time to improve the existing code.

The following requirements and recommendations are provided as a
guide.  They indicate which individual requirements are strict, and
which represent a preference and thus are negotiable or optional.
Please feel free to contact the LAMMPS core developers in case you
need additional explanations or clarifications, or you need assistance
in implementing the (strict) requirements for your contributions.
Requirements include:

* :ref:`Licensing requirements <ReqLicense>` (strict)
* :ref:`Integration testing <ReqIntegrationTesting>` (strict)
* :ref:`Documentation <ReqDocumentation>` (strict)
* :ref:`Programming language standards <ReqProgrammingStandards>` (strict)
* :ref:`Build system <ReqBuildSystem>` (strict)
* :ref:`Command or style names <ReqNaming>` (strict)
* :ref:`Programming style requirements <ReqProgrammingStyle>` (varied)
* :ref:`Examples <ReqExamples>` (preferred)
* :ref:`Error or warning messages and explanations <ReqErrorMessages>` (preferred)
* :ref:`Citation reminder <ReqCitation>` (optional)
* :ref:`Testing <ReqUnitTesting>` (optional)

.. _ReqLicense:

Licensing requirements (strict)
-------------------------------

Contributing authors agree when submitting a pull request that their
contributions can be distributed under the LAMMPS license conditions.
This is the GNU public license in version 2 (not 3 or later) for the
publicly distributed versions, e.g. on the LAMMPS homepage or on
GitHub.  We also have a version of LAMMPS under LGPL 2.1 terms which
is available on request; this will usually be the latest available or
a previous stable version with a few LGPL 2.1 incompatible files
removed.  More details are found on the :doc:`LAMMPS open-source
license page <Intro_opensource>`.

Your new source files should have the LAMMPS copyright and GPL notice,
followed by your name and email address at the top, like other
user-contributed LAMMPS source files.

Contributions may be under a different license as long as that license
does not conflict with the aforementioned terms.  Contributions that
use code with a conflicting license can be split into two parts:

1. the core parts (i.e. parts that must be in the `src` tree) that are
   licensed under compatible terms and bundled with the LAMMPS sources
2. an external library that must be downloaded and compiled (either
   separately or as part of the LAMMPS compilation)

Please note, that this split licensing mode may complicate including
the contribution in binary packages.

.. _ReqIntegrationTesting:

Integration testing (strict)
----------------------------

Where possible we use available continuous integration tools to search
for common programming mistakes, portability limitations, incompatible
formatting, and undesired side effects. Contributed code must pass the
automated tests on GitHub before it can be merged with the LAMMPS
distribution. These tests compile LAMMPS in a variety of environments
and settings and run the bundled unit tests.  At the discretion of the
LAMMPS developer managing the pull request, additional tests may be
activated that test for "side effects" on running a collection of
input decks and create consistent results.  The translation of the
documentation to HTML and PDF is also tested.

This means that contributed source code **must** compile with the most
current version of LAMMPS with ``-DLAMMPS_BIGBIG`` in addition to the
default setting of ``-DLAMMPS_SMALLBIG``.  The code needs to work
correctly in both cases, and also in serial and parallel using MPI.

Some "disruptive" changes may break tests and require updates to the
testing tools or scripts or tests themselves.  This is rare.  If in
doubt, contact the LAMMPS developer that is assigned to the pull
request.

.. _ReqDocumentation:

Documentation (strict)
----------------------

Contributions that add new styles or commands or augment existing ones
must include the corresponding new or modified documentation in
`ReStructuredText format <rst_>`_ (.rst files in the ``doc/src/``
folder). The documentation should be written in American English and the
.rst file must only use ASCII characters, so it can be cleanly
translated to PDF files (via `sphinx <https://www.sphinx-doc.org>`_ and
PDFLaTeX).  Special characters may be included via embedded math
expression typeset in a LaTeX subset.

.. _rst: https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html

When adding new commands, they need to be integrated into the sphinx
documentation system, and the corresponding command tables and lists
updated. When translating the documentation into html files there
should be no warnings. When adding a new package, some lists
describing packages must also be updated as well as a package specific
description added.  Likewise, if necessary, some package specific
build instructions should be included.

As appropriate, the text files with the documentation can include
inline mathematical expressions or figures (see ``doc/JPG`` for
examples).  Additional PDF files with further details may also be
included; see ``doc/PDF`` for examples.  The page should also include
literature citations as appropriate; see the bottom of
``doc/fix_nh.rst`` for examples and the earlier part of the same file
for how to format the cite itself.  Citation labels must be unique
across **all** .rst files.  The "Restrictions" section of the page
should indicate if your command is only available if LAMMPS is built
with the appropriate package.  See other command doc files for
examples of how to do this.

Please run at least "make html" and "make spelling" from within the
doc/src directory, and carefully inspect and proofread the resulting
HTML format doc page before submitting your code.  Upon submission of
a pull request, checks for error free completion of the HTML and PDF
build will be performed and also a spell check, a check for correct
anchors and labels, and a check for completeness of references to all
styles in their corresponding tables and lists is run.  In case the
spell check reports false positives, they can be added to the file
``doc/utils/sphinx-config/false_positives.txt``

Contributions that add or modify the library interface or "public"
APIs from the C++ code or the Fortran module must include suitable
doxygen comments in the source and corresponding changes to the
documentation sources for the "Programmer Guide" guide section of the
LAMMPS manual.

If your feature requires some more complex steps and explanations to
be used correctly or some external or bundled tools or scripts, we
recommend that you also contribute a :doc:`Howto document <Howto>`
providing some more background information and some tutorial material.
This can also be used to provide more in-depth explanations of models
that require use of multiple commands.

As a rule-of-thumb, the more clear and self-explanatory you make your
documentation, README files and examples, and the easier you make it
for people to get started, the more likely it is that users will try
out your new feature.

.. _ReqProgrammingStandards:

Programming language standards (strict)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The core of LAMMPS is written in C++11 in a style that can be mostly
described as "C with classes".  Advanced C++ features like operator
overloading or excessive use of templates are avoided with the intent to
keep the code readable to programmers that have limited C++ programming
experience.  C++ constructs are acceptable when they help improve the
readability and reliability of the code, e.g. when using the
`std::string` class instead of manipulating pointers and calling the
string functions of the C library.  In addition, a collection of
convenient :doc:`utility functions and classes <Developer_utils>` for
recurring tasks and a collection of :doc:`platform neutral functions
<Developer_platform>` for improved portability are provided.
Contributions with code requiring more recent C++ standards are only
accepted as packages with the post C++11 standard code confined to the
package so that it is optional.

Included Fortran code has to be compatible with the Fortran 2003
standard.  Since not all platforms supported by LAMMPS provide good
support for compiling Fortran files, it should be considered to rewrite
these parts as C++ code, if possible and thus allow for a wider adoption
of the contribution.  As of January 2023, all previously included
Fortran code for the LAMMPS executable has been replaced by equivalent
C++ code.

Python code must be compatible with Python 3.5 and later.  Large parts
of LAMMPS (including the :ref:`PYTHON package <PKG-PYTHON>`) are also
compatible with Python 2.7.  Compatibility with Python 2.7 is desirable,
but compatibility with Python 3.5 is **required**.

Compatibility with older programming language standards is very
important to maintain portability and availability of LAMMPS on many
platforms.  This applies especially to HPC cluster environments, which
tend to be running older software stacks and where LAMMPS users may be
required to use those older tools for access to advanced hardware
features or not have the option to install newer compilers or libraries.

.. _ReqBuildSystem:

Build system (strict)
---------------------

LAMMPS currently supports two build systems: one that is based on
:doc:`traditional Makefiles <Build_make>` and one that is based on
:doc:`CMake <Build_cmake>`.  Therefore, your contribution must be
compatible with and support both build systems.

For a single pair of header and implementation files that are an
independent feature, it is usually only required to add them to
``src/.gitignore``.

For traditional make, if your contributed files or package depend on
other LAMMPS style files or packages also being installed
(e.g. because your file is a derived class from the other LAMMPS
class), then an ``Install.sh`` file is also needed to check for those
dependencies and modifications to ``src/Depend.sh`` to trigger the checks.
See other README and Install.sh files in other directories as
examples.

Similarly, for CMake support, changes may need to be made to
``cmake/CMakeLists.txt``, some of the files in ``cmake/presets``, and
possibly a file with specific instructions needs to be added to
``cmake/Modules/Packages/``.  Please check out how this is handled for
existing packages and ask the LAMMPS developers if you need assistance.

.. _ReqNaming:

Command or style names, file names, and keywords (strict)
---------------------------------------------------------

All user-visible command or style names should be all lower case and
should only use letters, numbers, or forward slashes.  They should be
descriptive and initialisms should be avoided unless they are well
established (e.g. lj for Lennard-Jones).  For a compute style
"some/name" the source files must be called ``compute_some_name.h`` and
``compute_some_name.cpp``. The "include guard" in the header file would
then be ``LMP_COMPUTE_SOME_NAME_H`` and the class name
``ComputeSomeName``.

.. _ReqProgrammingStyle:

Programming style requirements (varied)
---------------------------------------

To maintain source code consistency across contributions from many
people, there are various programming style requirements for
contributions to LAMMPS.  Some of these requirements are strict and
must be followed, while others are only preferred and thus may be
skipped.  An in-depth discussion of the style guidelines is provided
in the :doc:`programming style doc page <Modify_style>`.

.. _ReqExamples:

Examples (preferred)
--------------------

For many new features, it is preferred that example scripts (simple,
small, fast to complete on 1 CPU) are included that demonstrate the
use of new or extended functionality. These are typically include
under the examples or examples/PACKAGES directory and are further
described on the :doc:`examples page <Examples>`.  Guidelines for
input scripts include:

- commands that generate output should be commented out (except when the
  output is the sole purpose or the feature, e.g. for a new compute)

- commands like :doc:`log <log>`, :doc:`echo <echo>`, :doc:`package
  <package>`, :doc:`processors <processors>`, :doc:`suffix <suffix>` may
  **not** be used in the input file (exception: "processors * * 1" or
  similar is acceptable when used to avoid unwanted domain decomposition
  of empty volumes)

- outside of the log files, no generated output should be included

- custom thermo_style settings may not include output measuring CPU or other
  time as it complicates comparisons between different runs

- input files should be named ``in.name``, data files should be named
  ``data.name`` and log files should be named ``log.version.name.<compiler>.<ncpu>``

- the total file size of all the inputs and outputs should be small

- where possible, potential files from the "potentials" folder or data
  file from other folders should be re-used through symbolic links

.. _ReqErrorMessages:

Error or warning messages and explanations (preferred)
------------------------------------------------------

.. versionchanged:: 4May2022

Starting with LAMMPS version 4 May 2022, the LAMMPS developers have
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
indicate the cause directly, and requiring no further lookup. This is
aided by the `{fmt} library <https://fmt.dev>`_ enabling Error class
methods that take a variable number of arguments and an error text that
will be treated like a {fmt} syntax format string. Error messages should
still preferably be kept to a single line or two lines at most.

For more complex explanations or errors that have multiple possible
reasons, a paragraph should be added to the `Error_details` page with an
error code reference (e.g. ``.. _err0001:``) then the utility function
:cpp:func:`utils::errorurl() <LAMMPS_NS::utils::errorurl>` can be used
to generate a URL that will directly lead to that paragraph.  An error
for missing arguments can be easily generated using the
:cpp:func:`utils::missing_cmd_args()
<LAMMPS_NS::utils::missing_cmd_args>` convenience function.
An example for this approach would be the
``src/read_data.cpp`` and ``src/atom.cpp`` files that implement the
:doc:`read_data <read_data>` and :doc:`atom_modify <atom_modify>`
commands and that may create :ref:`"Unknown identifier in data file" <err0001>`
errors that may have multiple possible reasons which complicates debugging,
and thus require some additional explanation.

The transformation of existing LAMMPS code to this new scheme is
ongoing.  Given the size of the LAMMPS code base, it will take a
significant amount of time to complete.  For new code, however,
following the new approach is strongly preferred.  The expectation is
that the new scheme will make understanding errors easier for LAMMPS
users, developers, and maintainers.

.. _ReqCitation:

Citation reminder (optional)
-----------------------------

If there is a paper of yours describing your feature (either the
algorithm/science behind the feature itself, or its initial usage, or
its implementation in LAMMPS), you can add the citation to the \*.cpp
source file.  See ``src/DIFFRACTION/compute_saed.cpp`` for an example.
A BibTeX format citation is stored in a string variable at the top of
the file, and a single line of code registering this variable is added
to the constructor of the class.  When your feature is used, then
LAMMPS (by default) will print the brief info and the DOI in the first
line to the screen and the full citation to the log file.

If there is additional functionality (which may have been added later)
described in a different publication, additional citation descriptions
may be added so long as they are only registered when the
corresponding keyword activating this functionality is used.

With these options, it is possible to have LAMMPS output a specific
citation reminder whenever a user invokes your feature from their
input script.  Please note that you should *only* use this for the
*most* relevant paper for a feature and a publication that you or your
group authored.  E.g. adding a citation in the source code for a paper
by Nose and Hoover if you write a fix that implements their integrator
is not the intended usage.  That kind of citation should just be
included in the documentation page you provide describing your
contribution.  If you are not sure what the best option would be,
please contact the LAMMPS developers for advice.

.. _ReqUnitTesting:

Testing (optional)
------------------

If your contribution contains new utility functions or a supporting
class (i.e. anything that does not depend on a LAMMPS object), new
unit tests should be added to a suitable folder in the ``unittest``
tree.  When adding a new LAMMPS style computing forces or selected
fixes, a ``.yaml`` file with a test configuration and reference data
should be added for the styles where a suitable tester program already
exists (e.g. pair styles, bond styles, etc.). Please see :ref:`this
section in the manual <testing>` for more information on how to
enable, run, and expand testing.
