.. index:: box

box command
===========

Syntax
""""""

.. parsed-literal::

   box keyword value ...

* one or more keyword/value pairs may be appended
* keyword = *tilt*

  .. parsed-literal::

       *tilt* value = *small* or *large*

Examples
""""""""

.. code-block:: LAMMPS

   box tilt large
   box tilt small

Description
"""""""""""

Set attributes of the simulation box.

For triclinic (non-orthogonal) simulation boxes, the *tilt* keyword
allows simulation domains to be created with arbitrary tilt factors,
e.g. via the :doc:`create_box <create_box>` or
:doc:`read_data <read_data>` commands.  Tilt factors determine how
skewed the triclinic box is; see the :doc:`Howto triclinic <Howto_triclinic>` doc page for a discussion of triclinic
boxes in LAMMPS.

LAMMPS normally requires that no tilt factor can skew the box more
than half the distance of the parallel box length, which is the 1st
dimension in the tilt factor (x for xz).  If *tilt* is set to
*small*\ , which is the default, then an error will be
generated if a box is created which exceeds this limit.  If *tilt*
is set to *large*\ , then no limit is enforced.  You can create
a box with any tilt factors you wish.

Note that if a simulation box has a large tilt factor, LAMMPS will run
less efficiently, due to the large volume of communication needed to
acquire ghost atoms around a processor's irregular-shaped sub-domain.
For extreme values of tilt, LAMMPS may also lose atoms and generate an
error.

Restrictions
""""""""""""

This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command or
:doc:`read_restart <read_restart>` command.

**Related commands:** none

Default
"""""""

The default value is tilt = small.
