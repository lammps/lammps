########################################
LAMMPS Documentation (|version| version)
########################################

LAMMPS stands for **L**\ arge-scale **A**\ tomic/**M**\ olecular
**M**\ assively **P**\ arallel **S**\ imulator.

LAMMPS is a classical molecular dynamics simulation code focusing on
materials modeling.  It was designed to run efficiently on parallel
computers and to be easy to extend and modify.  Originally developed at
Sandia National Laboratories, a US Department of Energy facility, LAMMPS
now includes contributions from many research groups and individuals
from many institutions.  Most of the funding for LAMMPS has come from
the US Department of Energy (DOE).  LAMMPS is open-source software
distributed under the terms of the GNU Public License Version 2 (GPLv2).

The `LAMMPS website <lws_>`_ has a variety of information about the
code.  It includes links to an online version of this manual, an
`online forum <https://www.lammps.org/forum.html>`_ where users can post
questions and discuss LAMMPS, and a `GitHub site
<https://github.com/lammps/lammps>`_ where all LAMMPS development is
coordinated.

----------

The content for this manual is part of the LAMMPS distribution.  The
online version always corresponds to the latest feature release version.
If needed, you can build a local copy of the manual as HTML pages or a
PDF file by following the steps on the :doc:`Build_manual` page.  If you
have difficulties viewing the pages, please :ref:`see this note
<webbrowser>`.

-----------

The manual is organized into three parts:

1. The :ref:`User Guide <user_documentation>` with information about how
   to obtain, configure, compile, install, and use LAMMPS,
2. the :ref:`Programmer Guide <programmer_documentation>` with
   information about how to use the LAMMPS library interface from
   different programming languages, how to modify and extend LAMMPS,
   the program design, internal programming interfaces, and code
   design conventions,
3. the :ref:`Command Reference <command_reference>` with detailed
   descriptions of all input script commands available in LAMMPS.

----------

.. only:: html

   After becoming familiar with LAMMPS, consider bookmarking
   :doc:`this page <Commands_all>`, since it gives quick access to
   tables with links to the documentation for all LAMMPS commands.

.. _lws: https://www.lammps.org

----------

.. _user_documentation:

************
User Guide
************

.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: User Guide
   :name: userdoc
   :includehidden:

   Intro
   Install
   Build
   Run_head
   Commands
   Packages
   Speed
   Howto
   Examples
   Tools
   Errors


.. _programmer_documentation:

******************
Programmer Guide
******************

.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: Programmer Guide
   :name: progdoc
   :includehidden:

   Library
   Python_head
   Modify
   Developer

*****************
Command Reference
*****************

.. _command_reference:
.. toctree::
   :name: reference
   :maxdepth: 1
   :caption: Command Reference

   commands_list
   fixes
   computes
   pairs
   bonds
   angles
   dihedrals
   impropers
   dumps
   fix_modify_atc_commands
   Bibliography

******************
Indices and tables
******************

.. only:: html

   * :ref:`genindex`
   * :ref:`search`

.. _webbrowser:
.. admonition:: Web Browser Compatibility
   :class: note

   The HTML version of the manual makes use of advanced features present
   in "modern" web browsers.  This can lead to incompatibilities with older
   web browsers (released more than 4 years ago) and specific vendor browsers
   (e.g. Internet Explorer on Windows; Microsoft Edge works well though)
   where parts of the pages are not rendered as expected (e.g. the layout is
   broken or mathematical expressions not typeset).  In that case we
   recommend to install/use a different/newer web browser or use
   the `PDF version of the manual <https://docs.lammps.org/Manual.pdf>`_.
