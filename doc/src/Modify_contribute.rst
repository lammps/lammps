Submitting new features for inclusion in LAMMPS
===============================================

We encourage users to submit new features or modifications for LAMMPS to
`the core developers <https://www.lammps.org/authors.html>`_ so they
can be added to the LAMMPS distribution. The preferred way to manage and
coordinate this is via the LAMMPS project on `GitHub
<https://github.com/lammps/lammps>`_.  Please see the :doc:`GitHub
Tutorial <Howto_github>` for a demonstration on how to do that.  An
alternative is to contact the LAMMPS developers or the indicated
developer of a package or feature directly and send in your contribution
via e-mail, but that can add a significant delay on getting your
contribution included, depending on how busy the respective developer
is, how complex a task it would be to integrate that code, and how
many - if any - changes are required before the code can be included.

For any larger modifications or programming project, you are encouraged
to contact the LAMMPS developers ahead of time in order to discuss
implementation strategies and coding guidelines. That will make it
easier to integrate your contribution and results in less work for
everybody involved.  You are also encouraged to search through the list
of `open issues on GitHub <https://github.com/lammps/lammps/issues>`_
and submit a new issue for a planned feature, so you would not duplicate
the work of others (and possibly get scooped by them) or have your work
duplicated by others.

For informal communication with (some of) the LAMMPS developers you may
ask to join the `LAMMPS developers on Slack
<https://lammps.slack.com>`_.  This slack work space is by invitation
only. Thus for access, please send an e-mail to ``slack@lammps.org``
explaining what part of LAMMPS you are working on.  Only discussions
related to LAMMPS development are tolerated, so this is **NOT** for
people that look for help with compiling, installing, or using
LAMMPS. Please contact the
`lammps-users mailing list <https://www.lammps.org/mail.html>`_ or the
`LAMMPS forum <https://www.lammps.org/forum.html>`_ for those purposes
instead.

How quickly your contribution will be integrated depends largely on how
much effort it will cause to integrate and test it, how many and what
kind of changes it requires to the core codebase, and of how much
interest it is to the larger LAMMPS community.  Please see below for a
checklist of typical requirements.  Once you have prepared everything,
see the :doc:`LAMMPS GitHub Tutorial <Howto_github>` page for
instructions on how to submit your changes or new files through a GitHub
pull request.  If you prefer to submit patches or full files, you should
first make certain, that your code works correctly with the latest
patch-level version of LAMMPS and contains all bug fixes from it.  Then
create a gzipped tar file of all changed or added files or a
corresponding patch file using 'diff -u' or 'diff -c' and compress it
with gzip.  Please only use gzip compression, as this works well and is
available on all platforms.

If the new features/files are broadly useful we may add them as core
files to LAMMPS or as part of a :doc:`package <Packages_list>`.  All
packages are listed and described on the :doc:`Packages details
<Packages_details>` doc page.

Note that by providing us files to release, you are agreeing to make
them open-source, i.e. we can release them under the terms of the GPL
(version 2), used as a license for the rest of LAMMPS.  And as part of
a LGPL (version 2.1) distribution that we make available to developers
on request only and with files that are authorized for that kind of
distribution removed (e.g. interface to FFTW).  See the
:doc:`LAMMPS license <Intro_opensource>` page for details.

.. note::

   If you prefer to actively develop and support your add-on feature
   yourself, then you may wish to make it available for download from
   your own website, as a user package that LAMMPS users can add to
   their copy of LAMMPS.  See the `Offsite LAMMPS packages and tools
   <https://www.lammps.org/offsite.html>`_ page of the LAMMPS web site
   for examples of groups that do this.  We are happy to advertise your
   package and web site from that page.  Simply email the `developers
   <https://www.lammps.org/authors.html>`_ with info about your package
   and we will post it there.  We recommend to name external packages
   USER-\<name\> so they can be easily distinguished from bundled packages
   that do not have the USER- prefix.

.. _lws: https://www.lammps.org

The previous sections of this page describe how to add new "style"
files of various kinds to LAMMPS.  Packages are simply collections of
one or more new class files which are invoked as a new style within a
LAMMPS input script.  If designed correctly, these additions typically
do not require changes to the main core of LAMMPS; they are simply
add-on files.  If you think your new feature requires non-trivial
changes in core LAMMPS files, you should `communicate with the
developers <https://www.lammps.org/authors.html>`_, since we may or
may not want to include those changes for some reason.  An example of a
trivial change is making a parent-class method "virtual" when you derive
a new child class from it.

Here is a checklist of steps you need to follow to submit a single file
or package for our consideration.  Following these steps will save
both you and us time. Please have a look at the existing files in
packages in the src directory for examples. If you are uncertain, please ask.

* All source files you provide must compile with the most current
  version of LAMMPS with multiple configurations. In particular you
  need to test compiling LAMMPS from scratch with -DLAMMPS_BIGBIG
  set in addition to the default -DLAMMPS_SMALLBIG setting. Your code
  will need to work correctly in serial and in parallel using MPI.

* For consistency with the rest of LAMMPS and especially, if you want
  your contribution(s) to be added to main LAMMPS code or one of its
  standard packages, it needs to be written in a style compatible with
  other LAMMPS source files. This means: 2-character indentation per
  level, **no tabs**, no lines over 100 characters. I/O is done via
  the C-style stdio library (mixing of stdio and iostreams is generally
  discouraged), class header files should not import any system headers
  outside of <cstdio>, STL containers should be avoided in headers,
  system header from the C library should use the C++-style names
  (<cstdlib>, <cstdio>, or <cstring>) instead of the C-style names
  <stdlib.h>, <stdio.h>, or <string.h>), and forward declarations
  used where possible or needed to avoid including headers.
  All added code should be placed into the LAMMPS_NS namespace or a
  sub-namespace; global or static variables should be avoided, as they
  conflict with the modular nature of LAMMPS and the C++ class structure.
  Header files must **not** import namespaces with *using*\ .
  This all is so the developers can more easily understand, integrate,
  and maintain your contribution and reduce conflicts with other parts
  of LAMMPS.  This basically means that the code accesses data
  structures, performs its operations, and is formatted similar to other
  LAMMPS source files, including the use of the error class for error
  and warning messages.

* To simplify reformatting contributed code in a way that is compatible
  with the LAMMPS formatting styles, you can use clang-format (version 8
  or later).  The LAMMPS distribution includes a suitable ``.clang-format``
  file which will be applied if you run ``clang-format -i some_file.cpp``
  on your files inside the LAMMPS src tree.  Please only reformat files
  that you have contributed.  For header files containing a
  ``SomeStyle(keyword, ClassName)`` macros it is required to have this
  macro embedded with a pair of ``// clang-format off``, ``// clang-format on``
  commends and the line must be terminated with a semi-colon (;).
  Example:

  .. code-block:: c++

     #ifdef COMMAND_CLASS
     // clang-format off
     CommandStyle(run,Run);
     // clang-format on
     #else

     #ifndef LMP_RUN_H
     [...]

  You may also use ``// clang-format on/off`` throughout your file
  to protect sections of the file from being reformatted.

* Please review the list of :doc:`available Packages <Packages_details>`
  to see if your contribution could be added to be added to one of them.
  It should fit into the general purposed of that package.  If it does not
  fit well, it can be added to one of the EXTRA- packages or the MISC package.

* If your contribution has several related features that are not covered
  by one of the existing packages or is dependent on a library (bundled
  or external), it is best to make it a package directory with a name
  like FOO.  In addition to your new files, the directory should contain
  a README text file.  The README should contain your name and contact
  information and a brief description of what your new package does.  If
  your files depend on other LAMMPS style files also being installed
  (e.g. because your file is a derived class from the other LAMMPS
  class), then an Install.sh file is also needed to check for those
  dependencies and modifications to src/Depend.sh to trigger the checks.
  See other README and Install.sh files in other directories as examples.
  Similarly for CMake support changes need to be made to cmake/CMakeLists.txt,
  the files in cmake/presets, and possibly a file to cmake/Modules/Packages/
  added.  Please check out how this is handled for existing packages and
  ask the LAMMPS developers if you need assistance.  Please submit a pull
  request on GitHub or send us a tarball of this FOO directory and all
  modified files.  Pull requests are strongly encouraged since they greatly
  reduce the effort required to integrate a contribution and simplify the
  process of adjusting the contributed code to cleanly fit into the
  LAMMPS distribution.

* Your new source files need to have the LAMMPS copyright, GPL notice,
  and your name and email address at the top, like other
  user-contributed LAMMPS source files.  They need to create a class
  that is inside the LAMMPS namespace.  To simplify maintenance, we
  may ask to adjust the programming style and formatting style to closer
  match the rest of LAMMPS.  We bundle a clang-format configuration file
  that can help with adjusting the formatting, although this is not a
  strict requirement.

* You **must** also create a **documentation** file for each new command
  or style you are adding to LAMMPS.  For simplicity and convenience,
  the documentation of groups of closely related commands or styles may
  be combined into a single file.  This will be one file for a
  single-file feature.  For a package, it might be several files.  These
  are text files with a .rst extension using the `reStructuredText
  <rst_>`_ markup language, that are then converted to HTML and PDF
  using the `Sphinx <sphinx_>`_ documentation generator tool.  Running
  Sphinx with the included configuration requires Python 3.x.
  Configuration settings and custom extensions for this conversion are
  included in the source distribution, and missing python packages will
  be transparently downloaded into a virtual environment via pip. Thus,
  if your local system is missing required packages, you need access to
  the internet. The translation can be as simple as doing "make html
  pdf" in the doc folder.  As appropriate, the text files can include
  inline mathematical expression or figures (see doc/JPG for examples).
  Additional PDF files with further details (see doc/PDF for examples)
  may also be included.  The page should also include literature
  citations as appropriate; see the bottom of doc/fix_nh.rst for
  examples and the earlier part of the same file for how to format the
  cite itself.  Citation labels must be unique across all .rst files.
  The "Restrictions" section of the page should indicate if your
  command is only available if LAMMPS is built with the appropriate
  FOO package.  See other package doc files for examples of
  how to do this.  Please run at least "make html" and "make spelling"
  and carefully inspect and proofread the resulting HTML format doc page
  before submitting your code.  Upon submission of a pull request,
  checks for error free completion of the HTML and PDF build will be
  performed and also a spell check, a check for correct anchors and
  labels, and a check for completeness of references all styles in their
  corresponding tables and lists is run.  In case the spell check
  reports false positives they can be added to the file
  doc/utils/sphinx-config/false_positives.txt

* For a new package (or even a single command) you should include one or
  more example scripts demonstrating its use.  These should run in no
  more than a couple minutes, even on a single processor, and not require
  large data files as input.  See directories under examples/PACKAGES for
  examples of input scripts other users provided for their packages.
  These example inputs are also required for validating memory accesses
  and testing for memory leaks with valgrind

* If there is a paper of yours describing your feature (either the
  algorithm/science behind the feature itself, or its initial usage, or
  its implementation in LAMMPS), you can add the citation to the \*.cpp
  source file.  See src/EFF/atom_vec_electron.cpp for an example.
  A LaTeX citation is stored in a variable at the top of the file and
  a single line of code registering this variable is added to the
  constructor of the class.  If there is additional functionality (which
  may have been added later) described in a different publication,
  additional citation descriptions may be added for as long as they
  are only registered when the corresponding keyword activating this
  functionality is used.  With these options it is possible to have
  LAMMPS output a specific citation reminder whenever a user invokes
  your feature from their input script.  Note that you should only use
  this for the most relevant paper for a feature and a publication that
  you or your group authored.  E.g. adding a citation in the code for
  a paper by Nose and Hoover if you write a fix that implements their
  integrator is not the intended usage.  That kind of citation should
  just be included in the documentation page you provide describing
  your contribution.  If you are not sure what the best option would
  be, please contact the LAMMPS developers for advice.

Finally, as a general rule-of-thumb, the more clear and
self-explanatory you make your documentation and README files, and the
easier you make it for people to get started, e.g. by providing example
scripts, the more likely it is that users will try out your new feature.

.. _rst: https://docutils.readthedocs.io/en/sphinx-docs/user/rst/quickstart.html
.. _sphinx: https://sphinx-doc.org
