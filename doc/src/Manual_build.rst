Building the LAMMPS manual
**************************

Depending on how you obtained LAMMPS, the doc directory has up
to 6 sub-directories, 2 Nroff files, and optionally 2 PDF files
plus 2 e-book format files:

.. parsed-literal::

   src             # content files for LAMMPS documentation
   html            # HTML version of the LAMMPS manual (see html/Manual.html)
   utils           # tools and settings for building the documentation
   docenv          # virtualenv for processing the manual sources
   doctrees        # temporary data from processing the manual
   mathjax         # code and fonts for rendering math in html
   Manual.pdf      # large PDF version of entire manual
   Developer.pdf   # small PDF with info about how LAMMPS is structured
   LAMMPS.epub     # Manual in ePUB e-book format
   LAMMPS.mobi     # Manual in MOBI e-book format
   lammps.1        # man page for the lammps command
   msi2lmp.1       # man page for the msi2lmp command

If you downloaded LAMMPS as a tarball from the web site, the html folder
and the PDF files should be included.

If you downloaded LAMMPS from the public git repository, then the HTML
and PDF files are not included.  Instead you need to create them, in one
of two ways:

a. You can "fetch" the current HTML and PDF files from the LAMMPS web
   site.  Just type "make fetch".  This should download a html\_www
   directory and Manual\_www.pdf/Developer\_www.pdf files.  Note that if
   new LAMMPS features have been added more recently than the date of
   your LAMMPS version, the fetched documentation will include those
   changes (but your source code will not, unless you update your local
   repository).

b. You can build the HTML or PDF files yourself, by typing "make html"
   or "make pdf".  This requires various tools including Sphinx, git,
   and the MathJax javascript library, which the build process will attempt
   to download automatically into a virtual environment in the folder
   doc/docenv and the folder mathjax, respectively, if not already available.
   This download is required only once, unless you type "make clean-all".
   After that, viewing and  processing of the documentation can be done
   without internet access.  To generate the PDF version of the manual,
   the PDFLaTeX software and several LaTeX packages are required as well.
   However, those cannot be installed automatically at the moment.

----------

The generation of all documentation is managed by the Makefile in
the doc directory.

.. code-block:: bash

   Documentation Build Options:

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

To run the HTML documentation build toolchain, Python 3 and virtualenv
have to be installed.  Here are instructions for common setups:

Ubuntu
------

.. code-block:: bash

   sudo apt-get install python-virtualenv

Fedora (up to version 21) and Red Hat Enterprise Linux or CentOS (up to version 7.x)
------------------------------------------------------------------------------------

.. code-block:: bash

   sudo yum install python3-virtualenv

Fedora (since version 22)
-------------------------

.. code-block:: bash

   sudo dnf install python3-virtualenv

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

----------

Installing prerequisites for epub build
=======================================

ePUB
----

Same as for HTML. This uses mostly the same tools and configuration
files as the HTML tree. In addition it uses LaTeX to convert embedded
math expressions transparently into embedded images.

For converting the generated ePUB file to a MOBI format file
(for e-book readers, like Kindle, that cannot read ePUB), you
also need to have the 'ebook-convert' tool from the "calibre"
software installed. `http://calibre-ebook.com/ <http://calibre-ebook.com/>`_
You first create the ePUB file and then convert it with 'make mobi'
On the Kindle readers in particular, you also have support for
PDF files, so you could download and view the PDF version as an alternative.
