.. index:: compute nbond/atom

compute nbond/atom command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID nbond/atom keyword value

* ID, group-ID are documented in :doc:`compute <compute>` command
* nbond/atom = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *bond/type*

  .. parsed-literal::

       *bond/type* value = *btype*
         *btype* = bond type included in count

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all nbond/atom
   compute 1 all nbond/atom bond/type 2

Description
"""""""""""

.. versionadded:: 4May2022

Define a computation that computes the number of bonds each atom is
part of.  Bonds which are broken are not counted in the tally.  See
the :doc:`Howto broken bonds <Howto_bpm>` page for more information.
The number of bonds will be zero for atoms not in the specified
compute group. This compute does not depend on Newton bond settings.

If the keyword *bond/type* is specified, only bonds of *btype* are
counted.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

Restrictions
""""""""""""

This compute is part of the BPM package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

Default
"""""""

none
