.. index:: fix planeforce

fix planeforce command
======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID planeforce x y z

* ID, group-ID are documented in :doc:`fix <fix>` command
* planeforce = style name of this fix command
* x y z = 3-vector that is normal to the plane

Examples
""""""""


.. parsed-literal::

   fix hold boundary planeforce 1.0 0.0 0.0

Description
"""""""""""

Adjust the forces on each atom in the group so that only the
components of force in the plane specified by the normal vector
(x,y,z) remain.  This is done by subtracting out the component of
force perpendicular to the plane.

If the initial velocity of the atom is 0.0 (or in the plane), then it
should continue to move in the plane thereafter.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix lineforce <fix_lineforce>`

**Default:** none


