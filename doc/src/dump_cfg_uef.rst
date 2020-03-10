.. index:: dump cfg/uef

dump cfg/uef command
====================

Syntax
""""""


.. parsed-literal::

   dump ID group-ID cfg/uef N file mass type xs ys zs args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* N = dump every this many timesteps
* file = name of file to write dump info to

  .. parsed-literal::

     args = same as args for :doc:`dump custom <dump>`



Examples
""""""""


.. parsed-literal::

   dump 1 all cfg/uef 10 dump.\*.cfg mass type xs ys zs
   dump 2 all cfg/uef 100 dump.\*.cfg mass type xs ys zs id c_stress

Description
"""""""""""

This command is used to dump atomic coordinates in the
reference frame of the applied flow field when
:doc:`fix nvt/uef <fix_nh_uef>` or
:doc:`fix npt/uef <fix_nh_uef>` or is used. Only the atomic
coordinates and frame-invariant scalar quantities
will be in the flow frame. If velocities are selected
as output, for example, they will not be in the same
reference frame as the atomic positions.

Restrictions
""""""""""""


This fix is part of the USER-UEF package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

This command can only be used when :doc:`fix nvt/uef <fix_nh_uef>`
or :doc:`fix npt/uef <fix_nh_uef>` is active.

Related commands
""""""""""""""""

:doc:`dump <dump>`,
:doc:`fix nvt/uef <fix_nh_uef>`

**Default:** none
