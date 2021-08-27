.. index:: compute nbonds/atom

compute nbonds/atom command
=======================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID nbonds/atom 

* ID, group-ID are documented in :doc:`compute <compute>` command
* nbonds/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all nbonds/atom

Description
"""""""""""

Define a computation that computes the number of bonds per-atom.
Bonds which are broken are not counted in the tally.
See :doc:`bond_style quartic <bond_quartic>` or the
:doc:`Howto bpm <Howto_bpm>` page. The number of bonds will be zero 
for atoms not in the specified compute group.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

Restrictions
""""""""""""

This fix can only be used if LAMMPS was built with the BPM
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

Default
"""""""

none

----------

