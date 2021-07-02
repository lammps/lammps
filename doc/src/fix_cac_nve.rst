.. index:: fix cac/nve

fix cac/nve command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID cac/nve

* ID, group-ID are documented in :doc:`fix <fix>` command
* cac/nve = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all cac/nve

Description
"""""""""""

Perform constant NVE integration to update position and velocity for
atoms and finite elements in the group each timestep. Note that there
is currently no formal reason to believe that this conserves energy 
for the CAC method; it is however free of external input that would supposedly
perturb the energy of the system.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

this fix requires a CAC atom style

**Default:** none
