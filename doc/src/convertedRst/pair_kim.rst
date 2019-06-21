.. index:: pair\_style kim

pair\_style kim command
=======================

Syntax
""""""


.. parsed-literal::

   pair_style kim model

model = name of KIM model (potential)

Examples
""""""""


.. parsed-literal::

   pair_style kim ex_model_Ar_P_LJ
   pair_coeff \* \* Ar Ar

Description
"""""""""""

This pair style is a wrapper on the `Knowledge Base for Interatomic Models (OpenKIM) <https://openkim.org>`_ repository of interatomic
potentials, so that they can be used by LAMMPS scripts.

Note that in LAMMPS lingo, a KIM model driver is a pair style
(e.g. EAM or Tersoff).  A KIM model is a pair style for a particular
element or alloy and set of parameters, e.g. EAM for Cu with a
specific EAM potential file.

See the current list of `KIM model drivers <https://openkim.org/browse/model-drivers/alphabetical>`_.

See the current list of all `KIM models <https://openkim.org/browse/models/by-model-drivers>`_

To use this pair style, you must first download and install the KIM
API library from the `OpenKIM website <https://openkim.org>`_.  The KIM
section of the :doc:`Packages details <Packages_details>` doc page has
instructions on how to do this with a simple make command, when
building LAMMPS.

See the examples/kim dir for an input script that uses a KIM model
(potential) for Lennard-Jones.


----------


The argument *model* is the name of the KIM model for a specific
potential as KIM defines it.  In principle, LAMMPS can invoke any KIM
model.  You should get an error or warning message from either LAMMPS
or KIM if there is an incompatibility.

Only a single pair\_coeff command is used with the *kim* style which
specifies the mapping of LAMMPS atom types to KIM elements.  This is
done by specifying N additional arguments after the \* \* in the
pair\_coeff command, where N is the number of LAMMPS atom types:

* N element names = mapping of KIM elements to atom types

As an example, imagine the KIM model supports Si and C atoms.  If your
LAMMPS simulation has 4 atom types and you want the 1st 3 to be Si,
and the 4th to be C, you would use the following pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* Si Si Si C

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to Si as
defined within KIM.  The final C argument maps LAMMPS atom type 4 to C
as defined within KIM.


----------


In addition to the usual LAMMPS error messages, the KIM library itself
may generate errors, which should be printed to the screen.  In this
case it is also useful to check the kim.log file for additional error
information.  The file kim.log should be generated in the same
directory where LAMMPS is running.

To download, build, and install the KIM library on your system, see
the lib/kim/README file.  Once you have done this and built LAMMPS
with the KIM package installed you can run the example input scripts
in examples/kim.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the :doc:`pair\_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since KIM stores the potential parameters.
Thus, you need to re-specify the pair\_style and pair\_coeff commands in
an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This pair style is part of the KIM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This current version of pair\_style kim is compatible with the
kim-api package version 2.0.0 and higher.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
