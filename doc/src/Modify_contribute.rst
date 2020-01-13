Submitting new features for inclusion in LAMMPS
===============================================

We encourage users to submit new features or modifications for LAMMPS
to `the core developers <http://lammps.sandia.gov/authors.html>`_ so they
can be added to the LAMMPS distribution. The preferred way to manage
and coordinate this is as of Fall 2016 via the LAMMPS project on
`GitHub <https://github.com/lammps/lammps>`_. An alternative is to
contact the LAMMPS developers or the indicated developer of a package
or feature directly and send in your contribution via e-mail.

For any larger modifications or programming project, you are
encouraged to contact the LAMMPS developers ahead of time, in order to
discuss implementation strategies and coding guidelines, that will
make it easier to integrate your contribution and result in less work
for everybody involved. You are also encouraged to search through the
list of `open issues on GitHub <https://github.com/lammps/lammps/issues>`_ and submit a new issue
for a planned feature, so you would not duplicate the work of others
(and possibly get scooped by them) or have your work duplicated by
others.

How quickly your contribution will be integrated depends largely on
how much effort it will cause to integrate and test it, how much it
requires changes to the core codebase, and of how much interest it is
to the larger LAMMPS community.  Please see below for a checklist of
typical requirements. Once you have prepared everything, see the
:doc:`Using GitHub with LAMMPS Howto <Howto_github>` doc page for instructions on how to
submit your changes or new files through a GitHub pull request. If you
prefer to submit patches or full files, you should first make certain,
that your code works correctly with the latest patch-level version of
LAMMPS and contains all bug fixes from it. Then create a gzipped tar
file of all changed or added files or a corresponding patch file using
'diff -u' or 'diff -c' and compress it with gzip. Please only use gzip
compression, as this works well on all platforms.

If the new features/files are broadly useful we may add them as core
files to LAMMPS or as part of a :doc:`standard package <Packages_standard>`.  Else we will add them as a
user-contributed file or :doc:`user package <Packages_user>`.  Examples
of user packages are in src sub-directories that start with USER.  The
USER-MISC package is simply a collection of (mostly) unrelated single
files, which is the simplest way to have your contribution quickly
added to the LAMMPS distribution.  All the standard and user packages
are listed and described on the :doc:`Packages details <Packages_details>` doc page.

Note that by providing us files to release, you are agreeing to make
them open-source, i.e. we can release them under the terms of the GPL,
used as a license for the rest of LAMMPS.  See the `Open source <http://lammps.sandia.gov/open_source.html>`_ page on the LAMMPS
website for details.

With user packages and files, all we are really providing (aside from
the fame and fortune that accompanies having your name in the source
code and on the `Authors page <http://lammps.sandia.gov/authors.html>`_
of the `LAMMPS WWW site <lws_>`_), is a means for you to distribute your
work to the LAMMPS user community, and a mechanism for others to
easily try out your new feature.  This may help you find bugs or make
contact with new collaborators.  Note that you're also implicitly
agreeing to support your code which means answer questions, fix bugs,
and maintain it if LAMMPS changes in some way that breaks it (an
unusual event).

.. note::

   If you prefer to actively develop and support your add-on
   feature yourself, then you may wish to make it available for download
   from your own website, as a user package that LAMMPS users can add to
   their copy of LAMMPS.  See the `Offsite LAMMPS packages and tools <http://lammps.sandia.gov/offsite.html>`_ page of the LAMMPS web
   site for examples of groups that do this.  We are happy to advertise
   your package and web site from that page.  Simply email the
   `developers <http://lammps.sandia.gov/authors.html>`_ with info about
   your package and we will post it there.

The previous sections of this doc page describe how to add new "style"
files of various kinds to LAMMPS.  Packages are simply collections of
one or more new class files which are invoked as a new style within a
LAMMPS input script.  If designed correctly, these additions typically
do not require changes to the main core of LAMMPS; they are simply
add-on files.  If you think your new feature requires non-trivial
changes in core LAMMPS files, you'll need to `communicate with the developers <http://lammps.sandia.gov/authors.html>`_, since we may or may
not want to make those changes.  An example of a trivial change is
making a parent-class method "virtual" when you derive a new child
class from it.

Here is a checklist of steps you need to follow to submit a single file
or user package for our consideration.  Following these steps will save
both you and us time. See existing files in packages in the src dir for
examples. If you are uncertain, please ask.

* All source files you provide must compile with the most current
  version of LAMMPS with multiple configurations. In particular you
  need to test compiling LAMMPS from scratch with -DLAMMPS\_BIGBIG
  set in addition to the default -DLAMMPS\_SMALLBIG setting. Your code
  will need to work correctly in serial and in parallel using MPI.
* For consistency with the rest of LAMMPS and especially, if you want
  your contribution(s) to be added to main LAMMPS code or one of its
  standard packages, it needs to be written in a style compatible with
  other LAMMPS source files. This means: 2-character indentation per
  level, **no tabs**\ , no lines over 80 characters. I/O is done via
  the C-style stdio library (mixing of stdio and iostreams is generally
  discouraged), class header files should not import any system headers
  outside of <cstdio>, STL containers should be avoided in headers,
  system header from the C library should use the C++-style names
  (<cstdlib>, <cstdio>, or <cstring>) instead of the C-style names
  <stdlib.h>, <stdio.h>, or <string.h>), and forward declarations
  used where possible or needed to avoid including headers.
  All added code should be placed into the LAMMPS\_NS namespace or a
  sub-namespace; global or static variables should be avoided, as they
  conflict with the modular nature of LAMMPS and the C++ class structure.
  Header files must **not** import namespaces with *using*\ .
  This all is so the developers can more easily understand, integrate,
  and maintain your contribution and reduce conflicts with other parts
  of LAMMPS.  This basically means that the code accesses data
  structures, performs its operations, and is formatted similar to other
  LAMMPS source files, including the use of the error class for error
  and warning messages.
* If you want your contribution to be added as a user-contributed
  feature, and it's a single file (actually a \*.cpp and \*.h file) it can
  rapidly be added to the USER-MISC directory.  Send us the one-line
  entry to add to the USER-MISC/README file in that dir, along with the
  2 source files.  You can do this multiple times if you wish to
  contribute several individual features.
* If you want your contribution to be added as a user-contribution and
  it is several related features, it is probably best to make it a user
  package directory with a name like USER-FOO.  In addition to your new
  files, the directory should contain a README text file.  The README
  should contain your name and contact information and a brief
  description of what your new package does.  If your files depend on
  other LAMMPS style files also being installed (e.g. because your file
  is a derived class from the other LAMMPS class), then an Install.sh
  file is also needed to check for those dependencies.  See other README
  and Install.sh files in other USER directories as examples.  Send us a
  tarball of this USER-FOO directory.

  Your new source files need to have the LAMMPS copyright, GPL notice,
  and your name and email address at the top, like other
  user-contributed LAMMPS source files.  They need to create a class
  that is inside the LAMMPS namespace.  If the file is for one of the

* USER packages, including USER-MISC, then we are not as picky about the
  coding style (see above).  I.e. the files do not need to be in the
  same stylistic format and syntax as other LAMMPS files, though that
  would be nice for developers as well as users who try to read your
  code.
* You **must** also create a **documentation** file for each new command or
  style you are adding to LAMMPS. For simplicity and convenience, the
  documentation of groups of closely related commands or styles may be
  combined into a single file.  This will be one file for a single-file
  feature.  For a package, it might be several files.  These are simple
  text files with a specific markup language, that are then auto-converted
  to HTML and PDF. The tools for this conversion are included in the
  source distribution, and the translation can be as simple as doing
  "make html pdf" in the doc folder.
  Thus the documentation source files must be in the same format and
  style as other \*.txt files in the lammps/doc/src directory for similar
  commands and styles; use one or more of them as a starting point.
  A description of the markup can also be found in
  lammps/doc/utils/txt2html/README.html
  As appropriate, the text files can include links to equations
  (see doc/Eqs/\*.tex for examples, we auto-create the associated JPG
  files), or figures (see doc/JPG for examples), or even additional PDF
  files with further details (see doc/PDF for examples).  The doc page
  should also include literature citations as appropriate; see the
  bottom of doc/fix\_nh.txt for examples and the earlier part of the same
  file for how to format the cite itself.  The "Restrictions" section of
  the doc page should indicate that your command is only available if
  LAMMPS is built with the appropriate USER-MISC or USER-FOO package.
  See other user package doc files for examples of how to do this. The
  prerequisite for building the HTML format files are Python 3.x and
  virtualenv, the requirement for generating the PDF format manual
  is the `htmldoc <http://www.htmldoc.org/>`_ software. Please run at least
  "make html" and carefully inspect and proofread the resulting HTML format
  doc page before submitting your code.
* For a new package (or even a single command) you should include one or
  more example scripts demonstrating its use.  These should run in no
  more than a couple minutes, even on a single processor, and not require
  large data files as input.  See directories under examples/USER for
  examples of input scripts other users provided for their packages.
  These example inputs are also required for validating memory accesses
  and testing for memory leaks with valgrind
* If there is a paper of yours describing your feature (either the
  algorithm/science behind the feature itself, or its initial usage, or
  its implementation in LAMMPS), you can add the citation to the \*.cpp
  source file.  See src/USER-EFF/atom\_vec\_electron.cpp for an example.
  A LaTeX citation is stored in a variable at the top of the file and a
  single line of code that references the variable is added to the
  constructor of the class.  Whenever a user invokes your feature from
  their input script, this will cause LAMMPS to output the citation to a
  log.cite file and prompt the user to examine the file.  Note that you
  should only use this for a paper you or your group authored.
  E.g. adding a cite in the code for a paper by Nose and Hoover if you
  write a fix that implements their integrator is not the intended
  usage.  That kind of citation should just be in the doc page you
  provide.


Finally, as a general rule-of-thumb, the more clear and
self-explanatory you make your documentation and README files, and the
easier you make it for people to get started, e.g. by providing example
scripts, the more likely it is that users will try out your new feature.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
