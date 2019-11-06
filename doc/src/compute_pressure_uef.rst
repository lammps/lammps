.. index:: compute pressure/uef

compute pressure/uef command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID pressure/uef temp-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/uef = style name of this compute command
* temp-ID = ID of compute that calculates temperature, can be NULL if not needed
* zero or more keywords may be appended
* keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or *improper* or *kspace* or *fix* or *virial*

Examples
""""""""


.. parsed-literal::

   compute 1 all pressure/uef my_temp_uef
   compute 2 all pressure/uef my_temp_uef virial

Description
"""""""""""

This command is used to compute the pressure tensor in
the reference frame of the applied flow field when
:doc:`fix nvt/uef <fix_nh_uef>` or
:doc:`fix npt/uef <fix_nh_uef>` is used.
It is not necessary to use this command to compute the scalar
value of the pressure. A :doc:`compute pressure <compute_pressure>`
may be used for that purpose.

The keywords and output information are documented in
:doc:`compute\_pressure <compute_pressure>`.

Restrictions
""""""""""""


This fix is part of the USER-UEF package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

This command can only be used when :doc:`fix nvt/uef <fix_nh_uef>`
or :doc:`fix npt/uef <fix_nh_uef>` is active.

The kinetic contribution to the pressure tensor
will be accurate only when
the compute specified by *temp-ID* is a
:doc:`compute temp/uef <compute_temp_uef>`.

Related commands
""""""""""""""""

:doc:`compute pressure <compute_pressure>`,
:doc:`fix nvt/uef <fix_nh_uef>`,
:doc:`compute temp/uef <compute_temp_uef>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
