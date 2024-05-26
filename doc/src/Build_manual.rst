Build the LAMMPS documentation
==============================

Depending on how you obtained LAMMPS and whether you have built the
manual yourself, this directory has a number of subdirectories and
files. Here is a list with descriptions:

.. code-block:: bash

   README           # brief info about the documentation
   src              # content files for LAMMPS documentation
   html             # HTML version of the LAMMPS manual (see html/Manual.html)
   utils            # tools and settings for building the documentation
   lammps.1         # man page for the lammps command
   msi2lmp.1        # man page for the msi2lmp command
   Manual.pdf       # large PDF version of entire manual
   LAMMPS.epub      # Manual in ePUB e-book format
   LAMMPS.mobi      # Manual in MOBI e-book format
   docenv           # virtualenv folder for processing the manual sources
   doctrees         # temporary data from processing the manual
   doxygen          # doxygen configuration and output
   .gitignore       # list of files and folders to be ignored by git
   doxygen-warn.log # logfile with warnings from running doxygen
   github-development-workflow.md   # notes on the LAMMPS development workflow

If you downloaded LAMMPS as a tarball from `the LAMMPS website <lws_>`_,
the html folder and the PDF files should be included.

If you downloaded LAMMPS from the public git repository, then the HTML
and PDF files are not included.  You can build the HTML or PDF files yourself,
by typing ``make html``  or ``make pdf`` in the ``doc`` folder.  This requires
various tools and files.  Some of them have to be installed (see below).  For
the rest the build process will attempt to download and install them into
a python virtual environment and local folders.

A current version of the manual (latest feature release, that is the state
of the *release* branch) is is available online at:
`https://docs.lammps.org/ <https://docs.lammps.org/>`_.
A version of the manual corresponding to the ongoing development (that is
the state of the *develop* branch) is available online at:
`https://docs.lammps.org/latest/ <https://docs.lammps.org/latest/>`_
A version of the manual corresponding to the latest stable LAMMPS release
(that is the state of the *stable* branch) is available online at:
`https://docs.lammps.org/stable/ <https://docs.lammps.org/stable/>`_

Build using GNU make
--------------------

The LAMMPS manual is written in `reStructuredText <rst_>`_ format which
can be translated to different output format using the `Sphinx
<https://www.sphinx-doc.org/>`_ document generator tool.  It also
incorporates programmer documentation extracted from the LAMMPS C++
sources through the `Doxygen <https://doxygen.nl/>`_ program.  Currently
the translation to HTML, PDF (via LaTeX), ePUB (for many e-book readers)
and MOBI (for Amazon Kindle readers) are supported.  For that to work a
Python interpreter version 3.8 or later, the ``doxygen`` tools and
internet access to download additional files and tools are required.
This download is usually only required once or after the documentation
folder is returned to a pristine state with ``make clean-all``.

For the documentation build a python virtual environment is set up in
the folder ``doc/docenv`` and various python packages are installed into
that virtual environment via the ``pip`` tool.  For rendering embedded
LaTeX code also the `MathJax <https://www.mathjax.org/>`_ JavaScript
engine needs to be downloaded.  If you need to pass additional options
to the pip commands to work (e.g. to use a web proxy or to point to
additional SSL certificates) you can set them via the ``PIP_OPTIONS``
environment variable or uncomment and edit the ``PIP_OPTIONS`` setting
at beginning of the makefile.

The actual translation is then done via ``make`` commands in the doc
folder.  The following ``make`` commands are available:

.. code-block:: bash

   make html          # generate HTML in html dir using Sphinx
   make pdf           # generate PDF  as Manual.pdf using Sphinx and PDFLaTeX
   make epub          # generate LAMMPS.epub in ePUB format using Sphinx
   make mobi          # generate LAMMPS.mobi in MOBI format using ebook-convert

   make fasthtml      # generate approximate HTML in fasthtml dir using Sphinx
                      # some Sphinx extensions do not work correctly with this

   make clean         # remove intermediate RST files created by HTML build
   make clean-all     # remove entire build folder and any cached data

   make anchor_check  # check for duplicate anchor labels
   make style_check   # check for complete and consistent style lists
   make package_check # check for complete and consistent package lists
   make link_check    # check for broken or outdated URLs
   make spelling      # spell-check the manual

----------

Build using CMake
-----------------

It is also possible to create the HTML version (and **only** the HTML
version) of the manual within the :doc:`CMake build directory
<Build_cmake>`.  The reason for this option is to include the
installation of the HTML manual pages into the "install" step when
installing LAMMPS after the CMake build via ``cmake --build . --target
install``.  The documentation build is included in the default build
target, but can also be requested independently with
``cmake --build . --target doc``.  If you need to pass additional options
to the pip commands to work (e.g. to use a web proxy or to point to
additional SSL certificates) you can set them via the ``PIP_OPTIONS``
environment variable.

.. code-block:: bash

   -D BUILD_DOC=value       # yes or no (default)

----------

Prerequisites for HTML
----------------------

To run the HTML documentation build toolchain, python 3, git, doxygen,
and virtualenv have to be installed locally.  Here are instructions for
common setups:

.. tabs::

   .. tab:: Ubuntu

      .. code-block:: bash

         sudo apt-get install git doxygen

   .. tab:: RHEL or CentOS (Version 7.x)

      .. code-block:: bash

         sudo yum install git doxygen

   .. tab:: Fedora or RHEL/CentOS (8.x or later)

      .. code-block:: bash

         sudo dnf install git doxygen

   .. tab:: macOS

      *Python 3*

      If Python 3 is not available on your macOS system, you can
      download the latest Python 3 macOS package from
      `https://www.python.org <https://www.python.org>`_ and install it.
      This will install both Python 3 and pip3.

Prerequisites for PDF
---------------------

In addition to the tools needed for building the HTML format manual,
a working LaTeX installation with support for PDFLaTeX and a selection
of LaTeX styles/packages are required.  To run the PDFLaTeX translation
the ``latexmk`` script needs to be installed as well.

Prerequisites for ePUB and MOBI
-------------------------------

In addition to the tools needed for building the HTML format manual,
a working LaTeX installation with a few add-on LaTeX packages
as well as the ``dvipng`` tool are required to convert embedded
math expressions transparently into embedded images.

For converting the generated ePUB file to a MOBI format file (for e-book
readers, like Kindle, that cannot read ePUB), you also need to have the
``ebook-convert`` tool from the "calibre" software
installed. `https://calibre-ebook.com/ <https://calibre-ebook.com/>`_
Typing ``make mobi`` will first create the ePUB file and then convert
it.  On the Kindle readers in particular, you also have support for PDF
files, so you could download and view the PDF version as an alternative.


Instructions for Developers
---------------------------

When adding new styles or options to the LAMMPS code, corresponding
documentation is required and either existing files in the ``src``
folder need to be updated or new files added. These files are written in
`reStructuredText <rst_>`_ markup for translation with the Sphinx tool.

Before contributing any documentation, please check that both the HTML
and the PDF format documentation can translate without errors.  During
testing the html translation, you may use the ``make fasthtml`` command
which does an approximate translation (i.e. not all Sphinx features and
extensions will work), but runs very fast because it will only translate
files that have been changed since the last ``make fasthtml`` command.

Please also check the output to the console for any warnings or problems.  There will
be multiple tests run automatically:

- A test for correctness of all anchor labels and their references

- A test that all LAMMPS packages (= folders with sources in
  ``lammps/src``) are documented and listed.  A typical warning shows
  the name of the folder with the suspected new package code and the
  documentation files where they need to be listed:

  .. parsed-literal::

     Found 88 packages
     Package NEWPACKAGE missing in Packages_list.rst
     Package NEWPACKAGE missing in Packages_details.rst

- A test that only standard, printable ASCII text characters are used.
  This runs the command ``env LC_ALL=C grep -n '[^ -~]' src/*.rst`` and
  thus prints all offending lines with filename and line number
  prepended to the screen.  Special characters like Greek letters
  (:math:`\alpha~~\sigma~~\epsilon`), super- or subscripts
  (:math:`x^2~~\mathrm{U}_{LJ}`), mathematical expressions
  (:math:`\frac{1}{2}\mathrm{N}~~x\to\infty`), or the Angstrom symbol
  (:math:`\AA`) should be typeset with embedded LaTeX (like this
  ``:math:`\alpha \sigma \epsilon```, ``:math:`x^2 \mathrm{E}_{LJ}```,
  ``:math:`\frac{1}{2}\mathrm{N} x\to\infty```, or ``:math:`\AA```\ ).

- Embedded LaTeX is rendered in HTML output with `MathJax
  <https://www.mathjax.org/>`_ and in PDF output by passing the embedded
  text to LaTeX.  Some care has to be taken, though, since there are
  limitations which macros and features can be used in either mode, so
  it is recommended to always check whether any new or changed
  documentation does translate and render correctly with either output.

- A test whether all styles are documented and listed in their
  respective overview pages.  A typical output with warnings looks like this:

  .. parsed-literal::

     Parsed style names w/o suffixes from C++ tree in ../src:
        Angle styles:      21    Atom styles:       24
        Body styles:        3    Bond styles:       17
        Command styles:    41    Compute styles:   143
        Dihedral styles:   16    Dump styles:       26
        Fix styles:       223    Improper styles:   13
        Integrate styles:   4    Kspace styles:     15
        Minimize styles:    9    Pair styles:      234
        Reader styles:      4    Region styles:      8
     Compute style entry newcomp is missing or incomplete in Commands_compute.rst
     Compute style entry newcomp is missing or incomplete in compute.rst
     Fix style entry newfix is missing or incomplete in Commands_fix.rst
     Fix style entry newfix is missing or incomplete in fix.rst
     Pair style entry new is missing or incomplete in Commands_pair.rst
     Pair style entry new is missing or incomplete in pair_style.rst
     Found 6 issue(s) with style lists


In addition, there is the option to run a spellcheck on the entire
manual with ``make spelling``.  This requires `a library called enchant
<https://github.com/AbiWord/enchant>`_.  To avoid printing out *false
positives* (e.g. keywords, names, abbreviations) those can be added to
the file ``lammps/doc/utils/sphinx-config/false_positives.txt``.

.. _lws: https://www.lammps.org
.. _rst: https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html
