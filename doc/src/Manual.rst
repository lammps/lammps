LAMMPS Documentation
####################

|version| version
*****************

:doc:`What is a LAMMPS version? <Manual_version>`

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel
Simulator.

LAMMPS is a classical molecular dynamics simulation code with a focus
on materials modeling.  It was designed to run efficiently on parallel
computers.  It was developed originally at Sandia National
Laboratories, a US Department of Energy facility.  The majority of
funding for LAMMPS has come from the US Department of Energy (DOE).
LAMMPS is an open-source code, distributed freely under the terms of
the GNU Public License (GPL).

The `LAMMPS website <lws_>`_ has a variety of information about the code.
It includes links to an on-line version of this manual, a `mailing list <https://lammps.sandia.gov/mail.html>`_ where users can post
questions, and a `GitHub site <https://github.com/lammps/lammps>`_ where
all LAMMPS development is coordinated.

----------

The content for this manual is part of the LAMMPS distribution.  You
can build a local copy of the Manual as HTML pages or a PDF file, by
following the steps on the :doc:`Manual build <Manual_build>` doc page.
There is also a `Developer.pdf <Developer.pdf>`_ document which gives
a brief description of the basic code structure of LAMMPS.

----------

Once you are familiar with LAMMPS, you may want to bookmark :doc:`this page <Commands_all>` since it gives quick access to a doc page for
every LAMMPS command.

.. _lws: https://lammps.sandia.gov

.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: User Documentation
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
   Modify
   Python_head
   Errors
   Manual_build

.. toctree::
   :maxdepth: 2
   :numbered: 3
   :caption: Programmer Documentation
   :name: progdoc
   :includehidden:

   pg_library
   pg_python
   pg_developer
   pg_base
   pg_utils

.. toctree::
   :caption: Index
   :name: index
   :hidden:

   commands_list
   fixes
   computes
   pairs
   bonds
   angles
   dihedrals
   impropers
   fix_modify_atc_commands

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
