.. index:: compute event/displace

compute event/displace command
==============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID event/displace threshold

* ID, group-ID are documented in :doc:`compute <compute>` command
* event/displace = style name of this compute command
* threshold = minimum distance any particle must move to trigger an event (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all event/displace 0.5

Description
"""""""""""

Define a computation that flags an "event" if any particle in the
group has moved a distance greater than the specified threshold
distance when compared to a previously stored reference state
(i.e. the previous event).  This compute is typically used in
conjunction with the :doc:`prd <prd>` and :doc:`tad <tad>` commands,
to detect if a transition
to a new minimum energy basin has occurred.

This value calculated by the compute is equal to 0 if no particle has
moved far enough, and equal to 1 if one or more particles have moved
further than the threshold distance.

.. note::

   If the system is undergoing significant center-of-mass motion,
   due to thermal motion, an external force, or an initial net momentum,
   then this compute will not be able to distinguish that motion from
   local atom displacements and may generate "false positives."

**Output info:**

This compute calculates a global scalar (the flag).  This value can be
used by any command that uses a global scalar value from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The scalar value calculated by this compute is "intensive".  The
scalar value will be a 0 or 1 as explained above.

Restrictions
""""""""""""

This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`prd <prd>`, :doc:`tad <tad>`

**Default:** none
