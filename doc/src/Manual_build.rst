Building the LAMMPS manual
**************************

Depending on how you obtained LAMMPS and whether you have built
the manual yourself, this directory has a varying number of
sub-directories and files. Here is a list with descriptions:

.. code-block:: bash

   README           # brief info about the documentation
   src              # content files for LAMMPS documentation
   html             # HTML version of the LAMMPS manual (see html/Manual.html)
   utils            # tools and settings for building the documentation
   lammps.1         # man page for the lammps command
   msi2lmp.1        # man page for the msi2lmp command
   Manual.pdf       # large PDF version of entire manual
   Developer.pdf    # small PDF with info about how LAMMPS is structured
   LAMMPS.epub      # Manual in ePUB e-book format
   LAMMPS.mobi      # Manual in MOBI e-book format
   docenv           # virtualenv for processing the manual sources
   doctrees         # temporary data from processing the manual
   mathjax          # code and fonts for rendering math in html
   doxygen          # doxygen configuration and output
   .gitignore       # list of files and folders to be ignored by git
   doxygen-warn.log # logfile with warnings from running doxygen
   github-development-workflow.md   # notes on the LAMMPS development workflow
   include-file-conventions.md      # notes on LAMMPS' include file conventions

If you downloaded LAMMPS as a tarball from `the LAMMPS website <lws_>`_,
the html folder and the PDF files should be included.

If you downloaded LAMMPS from the public git repository, then the HTML
and PDF files are not included.  Instead you need to create them, in one
of two ways:

a. You can "fetch" the current HTML and PDF files from the LAMMPS web
   site.  Just type ``make fetch``.  This should download a html_www
   directory and Manual_www.pdf/Developer_www.pdf files.  Note that if
   new LAMMPS features have been added more recently than the date of
   your LAMMPS version, the fetched documentation will include those
   changes (but your source code will not, unless you update your local
   repository).

b. You can build the HTML or PDF files yourself, by typing ``make html``
   or ``make pdf``.  This requires various tools and files.  Some of them
   have to be installed (more on that below). For the rest the build
   process will attempt to download and install them into a python
   virtual environment and local folders.  This download is required
   only once, unless you type ``make clean-all``.  After that, viewing and
   processing of the documentation can be done without internet access.

----------

The generation of all documentation is managed by the Makefile in the
doc directory. The following documentation related build options are
available:

.. code-block:: bash

   make html          # generate HTML in html dir using Sphinx
   make pdf           # generate 2 PDF files (Manual.pdf,Developer.pdf)
                      #   in doc dir via htmldoc and pdflatex
   make fetch         # fetch HTML doc pages and 2 PDF files from web site
                      #   as a tarball and unpack into html dir and 2 PDFs
   make epub          # generate LAMMPS.epub in ePUB format using Sphinx
   make mobi          # generate LAMMPS.mobi in MOBI format using ebook-convert
   make clean         # remove intermediate RST files created by HTML build
   make clean-all     # remove entire build folder and any cached data
   make anchor_check  # check for duplicate anchor labels
   make style_check   # check for complete and consistent style lists
   make package_check # check for complete and consistent package lists
   make spelling      # spell-check the manual

----------

Installing prerequisites for HTML build
=======================================

To run the HTML documentation build toolchain, python 3, git, doxygen,
and virtualenv have to be installed locally.  Here are instructions for
common setups:

Ubuntu
------

.. code-block:: bash

   sudo apt-get install python-virtualenv git doxygen

Fedora (up to version 21) and Red Hat Enterprise Linux or CentOS (up to version 7.x)
------------------------------------------------------------------------------------

.. code-block:: bash

   sudo yum install python3-virtualenv git doxygen

Fedora (since version 22)
-------------------------

.. code-block:: bash

   sudo dnf install python3-virtualenv git doxygen

MacOS X
-------

Python 3
^^^^^^^^

Download the latest Python 3 MacOS X package from
`https://www.python.org <https://www.python.org>`_
and install it.  This will install both Python 3
and pip3.

virtualenv
^^^^^^^^^^

Once Python 3 is installed, open a Terminal and type

.. code-block:: bash

   pip3 install virtualenv

This will install virtualenv from the Python Package Index.

Installing prerequisites for PDF build
======================================

In addition to the tools needed for building the HTML format manual,
a working LaTeX installation with support for PDFLaTeX and a selection
of LaTeX styles/packages are required.

Installing prerequisites for e-book reader builds
=================================================

In addition to the tools needed for building the HTML format manual,
a working LaTeX installation with a few add-on LaTeX packages
as well as the ``dvipng`` tool are required to convert embedded
math expressions transparently into embedded images.

For converting the generated ePUB file to a MOBI format file (for e-book
readers, like Kindle, that cannot read ePUB), you also need to have the
``ebook-convert`` tool from the "calibre" software
installed. `http://calibre-ebook.com/ <http://calibre-ebook.com/>`_
Typing ``make mobi`` will first create the ePUB file and then convert
it.  On the Kindle readers in particular, you also have support for PDF
files, so you could download and view the PDF version as an alternative.


Instructions for Developers
===========================

When adding new styles or options to the LAMMPS code, corresponding
documentation is required and either existing files in the ``src``
folder need to be updated or new files added. These files are written
in `reStructuredText <rst_>`_ markup for translation with the Sphinx tool.

Before contributing any documentation, please check that both the HTML
and the PDF format documentation can translate without errors. Please also
check the output to the console for any warnings or problems.  There will
be multiple tests run automatically:

- A test for correctness of all anchor labels and their references

- A test that all LAMMPS packages (= folders with sources in
  ``lammps/src``) are documented and listed.  A typical warning shows
  the name of the folder with the suspected new package code and the
  documentation files where they need to be listed:

  .. parsed-literal::

     Found 33 standard and 41 user packages
     Standard package NEWPACKAGE missing in Packages_standard.rst
     Standard package NEWPACKAGE missing in Packages_details.rst

- A test that only standard, printable ASCII text characters are used.
  This runs the command ``env LC_ALL=C grep -n '[^ -~]' src/*.rst`` and
  thus prints all offending lines with filename and line number
  prepended to the screen.  Special characters like the Angstrom
  :math:`\mathrm{\mathring{A}}` should be typeset with embedded math
  (like this ``:math:`\mathrm{\mathring{A}}```\ ).

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

.. _rst: https://docutils.readthedocs.io/en/sphinx-docs/user/rst/quickstart.html

.. _lws: https://lammps.sandia.gov
