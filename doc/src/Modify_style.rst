LAMMPS programming style and requirements
=========================================

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
