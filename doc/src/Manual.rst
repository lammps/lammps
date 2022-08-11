########################################
LAMMPS Documentation (|version| version)
########################################

LAMMPS stands for **L**\ arge-scale **A**\ tomic/**M**\ olecular
**M**\ assively **P**\ arallel **S**\ imulator.

LAMMPS is a classical molecular dynamics simulation code with a focus
on materials modeling.  It was designed to run efficiently on parallel
computers.  It was developed originally at Sandia National
Laboratories, a US Department of Energy facility.  The majority of
funding for LAMMPS has come from the US Department of Energy (DOE).
LAMMPS is an open-source code, distributed freely under the terms of
the GNU Public License Version 2 (GPLv2).

The `LAMMPS website <lws_>`_ has a variety of information about the
code.  It includes links to an on-line version of this manual, an
`online forum <https://www.lammps.org/forum.html>`_ where users can post
questions and discuss LAMMPS, and a `GitHub site
<https://github.com/lammps/lammps>`_ where all LAMMPS development is
coordinated.

----------

The content for this manual is part of the LAMMPS distribution.  The
online version always corresponds to the latest development version.
If needed, you can download or build a local copy of the manual as
HTML pages or a PDF file by following the steps on the
:doc:`Build_manual` page. If you have difficulties viewing the pages
please :ref:`see this note <webbrowser>`.

-----------

The manual is organized in three parts:
1) the :ref:`User Guide <user_documentation>` for how to install
and use LAMMPS, 2) the :ref:`Programmer Guide <programmer_documentation>`
for how to write programs using the LAMMPS library from different
programming languages and how to modify and extend LAMMPS, and 3) the
:ref:`Command Reference <command_reference>` which includes detailed
descriptions of all commands included in the LAMMPS code.

.. only:: html

   Once you are familiar with LAMMPS, you may want to bookmark
   :doc:`this page <Commands_all>` since it gives quick access
   the documentation for all LAMMPS commands.

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
